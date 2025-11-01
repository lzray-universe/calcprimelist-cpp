#include "cpu_info.h"

#include <algorithm>
#include <climits>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <thread>
#include <vector>
#include <string>
#include <cctype>
#include <limits>

#ifdef _WIN32
#define NOMINMAX
#include <windows.h>
#else
#include <dirent.h>
#include <sched.h>
#include <sys/types.h>
#endif

namespace calcprime {
namespace {

#ifdef _WIN32

using LogicalProcessorKey=std::pair<WORD,DWORD>;

std::size_t count_shared_physical(const CACHE_RELATIONSHIP&cache,const std::map<LogicalProcessorKey,unsigned>&logical_to_core,std::size_t&logical_count) {
    std::set<unsigned>physical_cores;
    logical_count=0;
    const GROUP_AFFINITY*group_masks=reinterpret_cast<const GROUP_AFFINITY*>(&cache.GroupMask);
    for(WORD group_index=0;group_index<cache.GroupCount;++group_index) {
        const GROUP_AFFINITY&affinity=group_masks[group_index];
        KAFFINITY mask=affinity.Mask;
        for(DWORD bit=0;bit<sizeof(KAFFINITY)*CHAR_BIT;++bit) {
            if((mask&(static_cast<KAFFINITY>(1)<<bit))==0) {
                continue;
            }
            ++logical_count;
            LogicalProcessorKey key(affinity.Group,bit);
            auto it=logical_to_core.find(key);
            if(it!=logical_to_core.end()) {
                physical_cores.insert(it->second);
            }
        }
    }
    if(!physical_cores.empty()) {
        return physical_cores.size();
    }
    return logical_count;
}

CpuInfo detect_windows_cpu_info() {
    CpuInfo info;

    DWORD length=0;
    GetLogicalProcessorInformationEx(RelationProcessorCore,nullptr,&length);
    std::vector<std::uint8_t>buffer(length);
    std::map<LogicalProcessorKey,unsigned>logical_to_core;
    if(GetLogicalProcessorInformationEx(RelationProcessorCore,
                                         reinterpret_cast<PSYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX>(buffer.data()),
                                         &length)) {
        unsigned physical=0;
        unsigned logical=0;
        bool has_smt=false;
        auto*ptr=reinterpret_cast<PSYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX>(buffer.data());
        while(reinterpret_cast<std::uint8_t*>(ptr)<buffer.data()+length) {
            if(ptr->Relationship==RelationProcessorCore) {
                unsigned core_index=physical;
                ++physical;
                DWORD core_logical=0;
                for(size_t i=0;i<ptr->Processor.GroupCount;++i) {
                    const GROUP_AFFINITY&affinity=ptr->Processor.GroupMask[i];
                    KAFFINITY mask=affinity.Mask;
                    for(DWORD bit=0;bit<sizeof(KAFFINITY)*CHAR_BIT;++bit) {
                        if((mask&(static_cast<KAFFINITY>(1)<<bit))==0) {
                            continue;
                        }
                        ++core_logical;
                        LogicalProcessorKey key(affinity.Group,bit);
                        logical_to_core[key]=core_index;
                    }
                }
                logical+=core_logical;
                if(core_logical>1) {
                    has_smt=true;
                }
            }
            ptr=reinterpret_cast<PSYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX>(reinterpret_cast<std::uint8_t*>(ptr)+ptr->Size);
        }
        if(physical) {
            info.physical_cpus=physical;
        }
        if(logical) {
            info.logical_cpus=logical;
        }
        info.has_smt=has_smt;
    }

    length=0;
    GetLogicalProcessorInformationEx(RelationCache,nullptr,&length);
    buffer.assign(length,0);
    std::size_t min_l1=std::numeric_limits<std::size_t>::max();
    std::size_t min_l2=std::numeric_limits<std::size_t>::max();
    bool have_l1=false;
    bool have_l2=false;
    std::set<std::string>seen_caches;
    if(GetLogicalProcessorInformationEx(RelationCache,
                                         reinterpret_cast<PSYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX>(buffer.data()),
                                         &length)) {
        auto*ptr=reinterpret_cast<PSYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX>(buffer.data());
        while(reinterpret_cast<std::uint8_t*>(ptr)<buffer.data()+length) {
            if(ptr->Relationship==RelationCache) {
                const auto&cache=ptr->Cache;
                std::size_t logical_count=0;
                std::size_t shared_physical=count_shared_physical(cache,logical_to_core,logical_count);
                std::size_t divisor=shared_physical ? shared_physical : logical_count;
                if(divisor==0) {
                    divisor=1;
                }
                std::size_t per_core=static_cast<std::size_t>(cache.CacheSize)/divisor;
                if(per_core==0) {
                    per_core=static_cast<std::size_t>(cache.CacheSize);
                }
                if(cache.Level==1&&cache.Type==CacheData) {
                    min_l1=std::min(min_l1,per_core);
                    have_l1=true;
                } else if(cache.Level==2&&(cache.Type==CacheUnified||cache.Type==CacheData)) {
                    min_l2=std::min(min_l2,per_core);
                    have_l2=true;
                    std::ostringstream oss;
                    oss<<cache.Level<<':'<<static_cast<int>(cache.Type);
                    for(WORD group_index=0;group_index<cache.GroupCount;++group_index) {
                        const GROUP_AFFINITY&affinity=
                            reinterpret_cast<const GROUP_AFFINITY*>(&cache.GroupMask)[group_index];
                        oss<<':'<<affinity.Group<<':'<<static_cast<unsigned long long>(affinity.Mask);
                    }
                    bool inserted=seen_caches.insert(oss.str()).second;
                    if(inserted) {
                        info.l2_total_bytes+=static_cast<std::size_t>(cache.CacheSize);
                    }
                }
            }
            ptr=reinterpret_cast<PSYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX>(reinterpret_cast<std::uint8_t*>(ptr)+ptr->Size);
        }
    }

    if(have_l1) {
        info.l1_data_bytes=min_l1;
    }
    if(have_l2) {
        info.l2_bytes=min_l2;
        if(!info.l2_total_bytes) {
            std::size_t cores=info.physical_cpus ? info.physical_cpus : info.logical_cpus;
            if(cores==0) {
                cores=1;
            }
            if(info.l2_bytes>0) {
                if(info.l2_bytes>std::numeric_limits<std::size_t>::max()/cores) {
                    info.l2_total_bytes=std::numeric_limits<std::size_t>::max();
                } else {
                    info.l2_total_bytes=info.l2_bytes*cores;
                }
            }
        }
    }

    if(!info.logical_cpus) {
        info.logical_cpus=std::max(1u,std::thread::hardware_concurrency());
    }
    if(!info.physical_cpus) {
        info.physical_cpus=std::max(1u,info.logical_cpus/2);
    }
    if(info.physical_cpus==info.logical_cpus) {
        info.has_smt=false;
    }
    return info;
}

#else

unsigned count_affinity_cpus() {
    cpu_set_t set;
    CPU_ZERO(&set);
    if(sched_getaffinity(0,sizeof(set),&set)!=0) {
        return std::max(1u,std::thread::hardware_concurrency());
    }
    unsigned count=0;
    for(int i=0;i<CPU_SETSIZE;++i) {
        if(CPU_ISSET(i,&set)) {
            ++count;
        }
    }
    return count;
}

struct LinuxCpuTopologyEntry {
    int package_id=0;
    int core_id=0;
};

std::map<unsigned,LinuxCpuTopologyEntry>read_linux_cpu_topology(const std::string&base_path) {
    std::map<unsigned,LinuxCpuTopologyEntry>topology;
    DIR*dir=opendir(base_path.c_str());
    if(!dir) {
        return topology;
    }
    struct dirent*entry;
    while((entry=readdir(dir))!=nullptr) {
        if(std::strncmp(entry->d_name,"cpu",3)!=0) {
            continue;
        }
        char*end=nullptr;
        long cpu_id=std::strtol(entry->d_name+3,&end,10);
        if(!end||*end!='\0'||cpu_id<0) {
            continue;
        }
        unsigned logical_id=static_cast<unsigned>(cpu_id);
        std::string core_path=base_path+"/"+entry->d_name+"/topology/core_id";
        std::ifstream core_file(core_path);
        int core_id=static_cast<int>(cpu_id);
        if(!(core_file>>core_id)) {
            core_id=static_cast<int>(cpu_id);
        }
        std::string package_path=base_path+"/"+entry->d_name+"/topology/physical_package_id";
        std::ifstream package_file(package_path);
        int package_id=0;
        if(!(package_file>>package_id)) {
            package_id=0;
        }
        topology[logical_id]={package_id,core_id};
    }
    closedir(dir);
    return topology;
}

std::size_t parse_cache_size_string(std::string size_str) {
    if(size_str.empty()) {
        return 0;
    }
    while(!size_str.empty()&&std::isspace(static_cast<unsigned char>(size_str.back()))) {
        size_str.pop_back();
    }
    std::size_t factor=1;
    if(!size_str.empty()) {
        char suffix=size_str.back();
        if(suffix=='K'||suffix=='k') {
            factor=1024;
            size_str.pop_back();
        } else if(suffix=='M'||suffix=='m') {
            factor=1024*1024;
            size_str.pop_back();
        }
    }
    std::size_t value=std::strtoull(size_str.c_str(),nullptr,10);
    return value*factor;
}

std::vector<unsigned>parse_cpu_list(const std::string&list) {
    std::vector<unsigned>cpus;
    std::size_t pos=0;
    while(pos<list.size()) {
        while(pos<list.size()&&(std::isspace(static_cast<unsigned char>(list[pos]))||list[pos]==',')) {
            ++pos;
        }
        if(pos>=list.size()) {
            break;
        }
        std::size_t end=pos;
        while(end<list.size()&&list[end]!=','&&!std::isspace(static_cast<unsigned char>(list[end]))) {
            ++end;
        }
        if(end<=pos) {
            pos=end+1;
            continue;
        }
        std::string token=list.substr(pos,end-pos);
        auto dash=token.find('-');
        if(dash==std::string::npos) {
            char*token_end=nullptr;
            long value=std::strtol(token.c_str(),&token_end,10);
            if(token_end&&*token_end=='\0'&&value>=0) {
                cpus.push_back(static_cast<unsigned>(value));
            }
        } else {
            std::string first=token.substr(0,dash);
            std::string second=token.substr(dash+1);
            char*end_first=nullptr;
            char*end_second=nullptr;
            long start=std::strtol(first.c_str(),&end_first,10);
            long stop=std::strtol(second.c_str(),&end_second,10);
            if(end_first&&*end_first=='\0'&&end_second&&*end_second=='\0'&&start>=0&&stop>=start) {
                for(long value=start;value<=stop;++value) {
                    cpus.push_back(static_cast<unsigned>(value));
                }
            }
        }
        pos=end+1;
    }
    std::sort(cpus.begin(),cpus.end());
    cpus.erase(std::unique(cpus.begin(),cpus.end()),cpus.end());
    return cpus;
}

std::string canonical_cpu_list(const std::vector<unsigned>&cpus) {
    if(cpus.empty()) {
        return std::string();
    }
    std::ostringstream oss;
    for(std::size_t i=0;i<cpus.size();++i) {
        if(i!=0) {
            oss<<',';
        }
        oss<<cpus[i];
    }
    return oss.str();
}

std::size_t shared_physical_cores(const std::vector<unsigned>&cpus,const std::map<unsigned,LinuxCpuTopologyEntry>&topology) {
    std::set<std::pair<int,int>>unique;
    for(unsigned cpu : cpus) {
        auto it=topology.find(cpu);
        if(it!=topology.end()) {
            unique.emplace(it->second.package_id,it->second.core_id);
        } else {
            unique.emplace(-1,static_cast<int>(cpu));
        }
    }
    return unique.size();
}

struct CacheIdentifier {
    int level;
    std::string type;
    std::string id;
    std::string shared;

    bool operator<(const CacheIdentifier&other) const {
        if(level!=other.level) {
            return level<other.level;
        }
        if(type!=other.type) {
            return type<other.type;
        }
        if(id!=other.id) {
            return id<other.id;
        }
        return shared<other.shared;
    }
};

CpuInfo detect_linux_cpu_info() {
    CpuInfo info;
    unsigned affinity=count_affinity_cpus();
    unsigned logical=std::max(affinity,std::thread::hardware_concurrency());
    if(logical==0) {
        logical=1;
    }
    info.logical_cpus=logical;

    const std::string topology_base="/sys/devices/system/cpu";
    auto topology=read_linux_cpu_topology(topology_base);
    if(topology.empty()) {
        topology[0]={};
    }

    std::set<std::pair<int,int>>unique_cores;
    for(const auto&entry : topology) {
        unique_cores.emplace(entry.second.package_id,entry.second.core_id);
    }
    unsigned physical=static_cast<unsigned>(unique_cores.size());
    if(physical==0) {
        physical=logical;
    }
    info.physical_cpus=physical;
    info.has_smt=physical<logical;

    std::set<CacheIdentifier>seen;
    std::size_t min_l1=std::numeric_limits<std::size_t>::max();
    std::size_t min_l2=std::numeric_limits<std::size_t>::max();
    bool have_l1=false;
    bool have_l2=false;

    for(const auto&entry : topology) {
        unsigned cpu_id=entry.first;
        std::string cache_dir_path=topology_base+"/cpu"+std::to_string(cpu_id)+"/cache";
        DIR*cache_dir=opendir(cache_dir_path.c_str());
        if(!cache_dir) {
            continue;
        }
        struct dirent*cache_entry;
        while((cache_entry=readdir(cache_dir))!=nullptr) {
            if(std::strncmp(cache_entry->d_name,"index",5)!=0) {
                continue;
            }
            std::string entry_base=cache_dir_path+"/"+cache_entry->d_name;

            std::ifstream level_file(entry_base+"/level");
            int level=0;
            if(!(level_file>>level)) {
                continue;
            }

            std::ifstream type_file(entry_base+"/type");
            std::string type;
            if(!(type_file>>type)) {
                continue;
            }
            std::string type_lower=type;
            std::transform(type_lower.begin(),type_lower.end(),type_lower.begin(),
                           [](unsigned char c) { return static_cast<char>(std::tolower(c));});

            std::ifstream size_file(entry_base+"/size");
            std::string size_str;
            if(!(size_file>>size_str)) {
                continue;
            }
            std::size_t size_bytes=parse_cache_size_string(size_str);
            if(!size_bytes) {
                continue;
            }

            std::ifstream shared_file(entry_base+"/shared_cpu_list");
            std::string shared_list;
            if(shared_file.is_open()) {
                std::getline(shared_file,shared_list);
            }
            auto shared_cpus=parse_cpu_list(shared_list);
            if(shared_cpus.empty()) {
                shared_cpus.push_back(cpu_id);
            }

            std::ifstream id_file(entry_base+"/id");
            std::string cache_id;
            if(!(id_file>>cache_id)) {
                cache_id.clear();
            }
            CacheIdentifier identifier{level,type_lower,cache_id,canonical_cpu_list(shared_cpus)};
            auto [it,inserted]=seen.insert(identifier);
            if(!inserted) {
                continue;
            }

            std::size_t shared_physical=shared_physical_cores(shared_cpus,topology);
            if(shared_physical==0) {
                shared_physical=shared_cpus.size();
            }
            std::size_t per_core=size_bytes;
            if(shared_physical>0) {
                per_core/=shared_physical;
                if(per_core==0) {
                    per_core=size_bytes;
                }
            }

            if(level==1&&type_lower=="data") {
                min_l1=std::min(min_l1,per_core);
                have_l1=true;
            } else if(level==2&&(type_lower=="unified"||type_lower=="data")) {
                min_l2=std::min(min_l2,per_core);
                have_l2=true;
                info.l2_total_bytes+=size_bytes;
            }
        }
        closedir(cache_dir);
    }

    if(have_l1) {
        info.l1_data_bytes=min_l1;
    }
    if(have_l2) {
        info.l2_bytes=min_l2;
        if(!info.l2_total_bytes) {
            std::size_t cores=info.physical_cpus ? info.physical_cpus : info.logical_cpus;
            if(cores==0) {
                cores=1;
            }
            if(info.l2_bytes>0) {
                if(info.l2_bytes>std::numeric_limits<std::size_t>::max()/cores) {
                    info.l2_total_bytes=std::numeric_limits<std::size_t>::max();
                } else {
                    info.l2_total_bytes=info.l2_bytes*cores;
                }
            }
        }
    }
    return info;
}

#endif

}

CpuInfo detect_cpu_info() {
#ifdef _WIN32
    return detect_windows_cpu_info();
#else
    return detect_linux_cpu_info();
#endif
}

unsigned effective_thread_count(const CpuInfo&info) {
    unsigned threads=info.physical_cpus;
    if(threads==0) {
        threads=info.logical_cpus;
    }
    if(threads==0) {
        threads=1;
    }
    return threads;
}

}
