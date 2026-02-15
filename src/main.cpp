#include "base_sieve.h"
#include "cpu_info.h"
#include "marker.h"
#include "popcnt.h"
#include "prime_count.h"
#include "segmenter.h"
#include "wheel_bitmap_count.h"
#include "wheel.h"
#include "writer.h"

#include<bit>
#include<algorithm>
#include<array>
#include<atomic>
#include<chrono>
#include<cinttypes>
#include<cctype>
#include<cmath>
#include<condition_variable>
#include<cstddef>
#include<cstdint>
#include<exception>
#include<fstream>
#include<cstdio>
#include<iomanip>
#include<iostream>
#include<limits>
#include<mutex>
#include<numeric>
#include<optional>
#include<sstream>
#include<stdexcept>
#include<string>
#include<thread>
#include<vector>

#ifdef _WIN32
#define NOMINMAX
#include<windows.h>
#include<winreg.h>
#else
#include<sys/utsname.h>
#endif

namespace calcprime{

struct Options{
	std::uint64_t from=0;
	std::uint64_t to=0;
	bool has_to=false;
	bool count_only=true;
	bool print_primes=false;
	std::optional<std::uint64_t> nth;
	unsigned threads=0;
	WheelType wheel=WheelType::Mod30;
	std::size_t segment_bytes=0;
	std::size_t tile_bytes=0;
	std::string output_path;
	PrimeOutputFormat output_format=PrimeOutputFormat::Text;
	bool use_zstd=false;
	bool show_time=false;
	bool show_stats=false;
	bool use_ml=false;
	bool use_wheel_bitmap=false;
	bool self_test=false;
	bool help=false;
	std::optional<std::uint64_t> test_value;
};

std::uint64_t parse_u64(const std::string&value){
	if(value.empty()){
		throw std::invalid_argument("invalid integer: "+value);
	}

	bool has_hex_prefix=value.size()>=2&&value[0]=='0'&&
						(value[1]=='x'||value[1]=='X');
	auto exp_pos=has_hex_prefix?std::string::npos:value.find_first_of("eE");
	if(exp_pos!=std::string::npos){
		std::string mantissa_str=value.substr(0,exp_pos);
		std::string exponent_str=value.substr(exp_pos+1);
		if(mantissa_str.empty()||exponent_str.empty()){
			throw std::invalid_argument("invalid integer: "+value);
		}
		if(!std::all_of(mantissa_str.begin(),mantissa_str.end(),
						[](unsigned char ch){ return std::isdigit(ch)!=0; })){
			throw std::invalid_argument("invalid integer: "+value);
		}

		std::size_t mantissa_idx=0;
		std::uint64_t mantissa=0;
		try{
			mantissa=std::stoull(mantissa_str,&mantissa_idx,10);
		}catch(const std::exception&){
			throw std::invalid_argument("invalid integer: "+value);
		}
		if(mantissa_idx!=mantissa_str.size()){
			throw std::invalid_argument("invalid integer: "+value);
		}

		std::size_t exponent_idx=0;
		long long exponent=0;
		try{
			exponent=std::stoll(exponent_str,&exponent_idx,10);
		}catch(const std::exception&){
			throw std::invalid_argument("invalid integer: "+value);
		}
		if(exponent_idx!=exponent_str.size()){
			throw std::invalid_argument("invalid integer: "+value);
		}
		if(exponent<0){
			throw std::invalid_argument("invalid integer: "+value);
		}
		if(mantissa==0){
			return 0;
		}
		if(exponent>19){
			throw std::invalid_argument("integer too large: "+value);
		}

		std::uint64_t result=mantissa;
		for(long long i=0;i<exponent;++i){
			if(result>std::numeric_limits<std::uint64_t>::max()/10ULL){
				throw std::invalid_argument("integer too large: "+value);
			}
			result*=10ULL;
		}
		return result;
	}

	std::size_t idx=0;
	std::uint64_t result=0;
	try{
		result=std::stoull(value,&idx,0);
	}catch(const std::exception&){
		throw std::invalid_argument("invalid integer: "+value);
	}
	if(idx!=value.size()){
		throw std::invalid_argument("invalid integer: "+value);
	}
	return result;
}

std::size_t parse_size(const std::string&value){
	if(value.empty()){
		throw std::invalid_argument("invalid size");
	}
	std::size_t idx=0;
	std::uint64_t base=std::stoull(value,&idx,0);
	std::uint64_t factor=1;
	if(idx<value.size()){
		char suffix=value[idx];
		switch(suffix){
		case 'k':
		case 'K':
			factor=1024;
			++idx;
			break;
		case 'm':
		case 'M':
			factor=1024*1024;
			++idx;
			break;
		case 'g':
		case 'G':
			factor=1024ull*1024ull*1024ull;
			++idx;
			break;
		default:
			break;
		}
	}
	if(idx!=value.size()){
		throw std::invalid_argument("invalid size suffix: "+value);
	}
	if(base>std::numeric_limits<std::uint64_t>::max()/factor){
		throw std::invalid_argument("size too large: "+value);
	}
	std::uint64_t result=base*factor;
	if(result>
	   static_cast<std::uint64_t>(std::numeric_limits<std::size_t>::max())){
		throw std::invalid_argument("size too large: "+value);
	}
	return static_cast<std::size_t>(result);
}

static void extract_segment_primes(const std::vector<std::uint64_t>&bitset,
								   std::uint64_t seg_low,
								   std::size_t bit_count,
								   std::vector<std::uint64_t>&primes){
	std::size_t remaining=bit_count;
	for(std::size_t word=0;word<bitset.size()&&remaining>0;++word){
		std::uint64_t valid_mask=std::numeric_limits<std::uint64_t>::max();
		if(remaining<64){
			valid_mask=(1ULL<<remaining)-1ULL;
		}
		std::uint64_t prime_bits=(~bitset[word])&valid_mask;
		while(prime_bits){
			unsigned bit=std::countr_zero(prime_bits);
			std::uint64_t index=(static_cast<std::uint64_t>(word)<<6)+
								static_cast<std::uint64_t>(bit);
			primes.push_back(seg_low+(index<<1));
			prime_bits&=(prime_bits-1);
		}
		if(remaining>=64){
			remaining-=64;
		}else{
			remaining=0;
		}
	}
}

void parse_output_format(Options&opts,const std::string&fmt){
	if(fmt=="text"){
		opts.output_format=PrimeOutputFormat::Text;
	}else if(fmt=="binary"){
		opts.output_format=PrimeOutputFormat::Binary;
	}else if(fmt=="delta16"){
		opts.output_format=PrimeOutputFormat::Delta16;
	}else if(fmt=="zstd"||fmt=="zstd+delta"){
		opts.output_format=PrimeOutputFormat::Delta16;
		opts.use_zstd=true;
		std::fprintf(stderr,
					 "[calcprime] warning: --out-format=%s is deprecated; "
					 "use --out-format=delta16 --zstd.\n",
					 fmt.c_str());
	}else{
		throw std::invalid_argument("unsupported out-format: "+fmt);
	}
}

Options parse_options(int argc,char**argv){
	Options opts;
	for(int i=1;i<argc;++i){
		std::string arg=argv[i];
		static const std::string out_format_prefix="--out-format=";

		if(arg=="--help"||arg=="-h"){
			opts.help=true;
			return opts;
		}else if(arg.rfind(out_format_prefix,0)==0){
			parse_output_format(opts,arg.substr(out_format_prefix.size()));
		}else if(arg=="--from"){
			if(i+1>=argc){
				throw std::invalid_argument("--from requires a value");
			}
			opts.from=parse_u64(argv[++i]);
		}else if(arg=="--to"){
			if(i+1>=argc){
				throw std::invalid_argument("--to requires a value");
			}
			opts.to=parse_u64(argv[++i]);
			opts.has_to=true;
		}else if(arg=="--count"){
			opts.count_only=true;
		}else if(arg=="--print"){
			opts.print_primes=true;
			opts.count_only=false;
		}else if(arg=="--nth"){
			if(i+1>=argc){
				throw std::invalid_argument("--nth requires a value");
			}
			opts.nth=parse_u64(argv[++i]);
			opts.count_only=false;
		}else if(arg=="--threads"){
			if(i+1>=argc){
				throw std::invalid_argument("--threads requires a value");
			}
			opts.threads=static_cast<unsigned>(parse_u64(argv[++i]));
		}else if(arg=="--wheel"){
			if(i+1>=argc){
				throw std::invalid_argument("--wheel requires a value");
			}
			std::string w=argv[++i];
			if(w=="30"){
				opts.wheel=WheelType::Mod30;
			}else if(w=="210"){
				opts.wheel=WheelType::Mod210;
			}else if(w=="1155"){
				opts.wheel=WheelType::Mod1155;
			}else{
				throw std::invalid_argument("unsupported wheel: "+w);
			}
		}else if(arg=="--segment"){
			if(i+1>=argc){
				throw std::invalid_argument("--segment requires a value");
			}
			opts.segment_bytes=parse_size(argv[++i]);
		}else if(arg=="--tile"){
			if(i+1>=argc){
				throw std::invalid_argument("--tile requires a value");
			}
			opts.tile_bytes=parse_size(argv[++i]);
		}else if(arg=="--out"){
			if(i+1>=argc){
				throw std::invalid_argument("--out requires a path");
			}
			opts.output_path=argv[++i];
		}else if(arg=="--out-format"){
			if(i+1>=argc){
				throw std::invalid_argument("--out-format requires a value");
			}
			parse_output_format(opts,argv[++i]);
		}else if(arg=="--zstd"){
			opts.use_zstd=true;
		}else if(arg=="--time"){
			opts.show_time=true;
		}else if(arg=="--stats"){
			opts.show_stats=true;
		}else if(arg=="--ml"){
			opts.use_ml=true;
		}else if(arg=="--wheel-bitmap"){
			opts.use_wheel_bitmap=true;
		}else if(arg=="--stest"){
			opts.self_test=true;
		}else if(arg=="--test"){
			if(i+1>=argc){
				throw std::invalid_argument("--test requires a value");
			}
			opts.test_value=parse_u64(argv[++i]);
		}else{
			throw std::invalid_argument("unknown option: "+arg);
		}
	}
	return opts;
}

void print_usage(){
	std::cout
		<<"prime-sieve --from A --to B [options]\n"
		<<"  --count             Count primes (default)\n"
		<<"  --print             Print primes in the interval\n"
		<<"  --nth K             Find the K-th prime in the interval\n"
		<<"  --threads N         Override thread count\n"
		<<"  --wheel 30|210|1155 Select wheel factorisation (default 30)\n"
		<<"  --segment BYTES     Override segment size\n"
		<<"  --tile BYTES        Override tile size\n"
		<<"  --out PATH          Write primes to file\n"
		<<"  --out-format FMT    Output format: text (default), binary, delta16\n"
		<<"                    Deprecated aliases: zstd, zstd+delta\n"
		<<"  --zstd              Compress output stream with zstd (if supported)\n"
		<<"  --time              Print elapsed time\n"
		<<"  --stats             Print configuration statistics\n"
		<<"  --ml                Use Meissel-Lehmer for counting (only with --count)\n"
		<<"  --wheel-bitmap      Force wheel-compressed count path\n"
		<<"                       (auto-enabled for large wheel=30 counts)\n"
		<<"  --stest             Run built-in benchmark (1e6..1e11, 10 runs)\n"
		<<"  --test N           Run a Miller-Rabin primality check for N\n";
}

struct SegmentResult{
	std::uint64_t count=0;
	std::vector<std::uint64_t> primes;
	std::atomic<bool> ready{false};
};

std::string trim_copy(const std::string&value){
	std::size_t first=0;
	while(first<value.size()&&
		  std::isspace(static_cast<unsigned char>(value[first]))){
		++first;
	}
	std::size_t last=value.size();
	while(last>first&&
		  std::isspace(static_cast<unsigned char>(value[last-1]))){
		--last;
	}
	return value.substr(first,last-first);
}

std::string format_with_commas(std::uint64_t value){
	std::string digits=std::to_string(value);
	std::string out;
	out.reserve(digits.size()+digits.size()/3);
	for(std::size_t i=0;i<digits.size();++i){
		if(i>0&&((digits.size()-i)%3==0)){
			out.push_back(',');
		}
		out.push_back(digits[i]);
	}
	return out;
}

std::string format_seconds(double value){
	std::ostringstream oss;
	oss<<std::fixed<<std::setprecision(6)<<value;
	return oss.str();
}

double median_seconds(std::vector<double>samples){
	if(samples.empty()){
		return 0.0;
	}
	std::sort(samples.begin(),samples.end());
	std::size_t mid=samples.size()/2;
	if((samples.size()&1U)!=0U){
		return samples[mid];
	}
	return (samples[mid-1]+samples[mid])*0.5;
}

std::string unquote_value(const std::string&value){
	if(value.size()>=2&&
	   ((value.front()=='"'&&value.back()=='"')||
		(value.front()=='\''&&value.back()=='\''))){
		return value.substr(1,value.size()-2);
	}
	return value;
}

struct HostInfo{
	std::string cpu;
	std::string os;
};

#ifdef _WIN32
std::string utf8_from_wide(const std::wstring&value){
	if(value.empty()){
		return std::string();
	}
	int needed=WideCharToMultiByte(CP_UTF8,0,value.data(),
								   static_cast<int>(value.size()),nullptr,0,
								   nullptr,nullptr);
	if(needed<=0){
		return std::string();
	}
	std::string out(static_cast<std::size_t>(needed),'\0');
	int converted=WideCharToMultiByte(CP_UTF8,0,value.data(),
									  static_cast<int>(value.size()),
									  out.data(),needed,nullptr,nullptr);
	if(converted!=needed){
		return std::string();
	}
	return out;
}

std::wstring read_registry_string(HKEY root,const wchar_t*subkey,
								  const wchar_t*value_name){
	DWORD type=0;
	DWORD bytes=0;
	LSTATUS status=RegGetValueW(root,subkey,value_name,RRF_RT_REG_SZ,&type,
								nullptr,&bytes);
	if(status!=ERROR_SUCCESS||bytes<sizeof(wchar_t)){
		return std::wstring();
	}
	std::vector<wchar_t> buffer(bytes/sizeof(wchar_t));
	status=RegGetValueW(root,subkey,value_name,RRF_RT_REG_SZ,&type,
						buffer.data(),&bytes);
	if(status!=ERROR_SUCCESS||buffer.empty()){
		return std::wstring();
	}
	return std::wstring(buffer.data());
}

bool read_registry_dword(HKEY root,const wchar_t*subkey,
						 const wchar_t*value_name,std::uint32_t&out){
	DWORD type=0;
	DWORD bytes=sizeof(DWORD);
	DWORD value=0;
	LSTATUS status=RegGetValueW(root,subkey,value_name,RRF_RT_REG_DWORD,&type,
								&value,&bytes);
	if(status!=ERROR_SUCCESS||bytes!=sizeof(DWORD)){
		return false;
	}
	out=static_cast<std::uint32_t>(value);
	return true;
}

std::string detect_windows_cpu_name(){
	const std::wstring cpu_name_w=read_registry_string(
		HKEY_LOCAL_MACHINE,
		L"HARDWARE\\DESCRIPTION\\System\\CentralProcessor\\0",
		L"ProcessorNameString");
	return trim_copy(utf8_from_wide(cpu_name_w));
}

std::string detect_windows_os_name(){
	const wchar_t*version_key=L"SOFTWARE\\Microsoft\\Windows NT\\CurrentVersion";
	std::string product=trim_copy(
		utf8_from_wide(read_registry_string(HKEY_LOCAL_MACHINE,version_key,
											L"ProductName")));
	std::string display=trim_copy(
		utf8_from_wide(read_registry_string(HKEY_LOCAL_MACHINE,version_key,
											L"DisplayVersion")));
	if(display.empty()){
		display=trim_copy(
			utf8_from_wide(read_registry_string(HKEY_LOCAL_MACHINE,
												version_key,L"ReleaseId")));
	}
	std::string build=trim_copy(
		utf8_from_wide(read_registry_string(HKEY_LOCAL_MACHINE,version_key,
											L"CurrentBuildNumber")));
	std::uint32_t ubr=0;
	bool has_ubr=read_registry_dword(HKEY_LOCAL_MACHINE,version_key,L"UBR",
									 ubr);
	if(!product.empty()&&product.find("Windows 10")!=std::string::npos&&
	   !build.empty()){
		try{
			unsigned long build_number=std::stoul(build);
			if(build_number>=22000){
				auto pos=product.find("Windows 10");
				if(pos!=std::string::npos){
					product.replace(pos,std::string("Windows 10").size(),
									"Windows 11");
				}
			}
		}catch(const std::exception&){
			// Keep original product string on parse failures.
		}
	}

	SYSTEM_INFO sys_info{};
	GetNativeSystemInfo(&sys_info);
	std::string arch="unknown";
	switch(sys_info.wProcessorArchitecture){
	case PROCESSOR_ARCHITECTURE_AMD64:
		arch="x64";
		break;
	case PROCESSOR_ARCHITECTURE_INTEL:
		arch="x86";
		break;
	case PROCESSOR_ARCHITECTURE_ARM64:
		arch="arm64";
		break;
	case PROCESSOR_ARCHITECTURE_ARM:
		arch="arm";
		break;
	default:
		break;
	}

	std::ostringstream os;
	os<<(product.empty()?"Microsoft Windows":product);
	if(!display.empty()){
		os<<", Version "<<display;
	}
	if(!build.empty()){
		os<<" (OS Build "<<build;
		if(has_ubr){
			os<<'.'<<ubr;
		}
		os<<')';
	}
	os<<", "<<arch;
	return os.str();
}

#else

std::string detect_linux_cpu_name(){
	std::ifstream cpuinfo("/proc/cpuinfo");
	if(!cpuinfo){
		return std::string();
	}
	std::string line;
	while(std::getline(cpuinfo,line)){
		auto colon=line.find(':');
		if(colon==std::string::npos){
			continue;
		}
		std::string key=trim_copy(line.substr(0,colon));
		if(key=="model name"||key=="Hardware"||key=="Processor"){
			return trim_copy(line.substr(colon+1));
		}
	}
	return std::string();
}

std::string detect_linux_os_name(){
	std::string pretty_name;
	{
		std::ifstream os_release("/etc/os-release");
		std::string line;
		while(std::getline(os_release,line)){
			auto eq=line.find('=');
			if(eq==std::string::npos){
				continue;
			}
			std::string key=trim_copy(line.substr(0,eq));
			if(key=="PRETTY_NAME"){
				pretty_name=unquote_value(trim_copy(line.substr(eq+1)));
				break;
			}
		}
	}

	struct utsname uts{};
	bool has_uname=(uname(&uts)==0);
	std::string arch=has_uname?trim_copy(uts.machine):std::string();
	std::string kernel=has_uname?trim_copy(uts.release):std::string();

	if(pretty_name.empty()){
		if(has_uname){
			pretty_name=trim_copy(uts.sysname);
		}else{
			pretty_name="Linux";
		}
	}

	std::ostringstream os;
	os<<pretty_name;
	if(!arch.empty()){
		os<<", "<<arch;
	}
	if(!kernel.empty()){
		os<<" (kernel "<<kernel<<')';
	}
	return os.str();
}

#endif

HostInfo detect_host_info(const CpuInfo&cpu_info){
	HostInfo host;
#ifdef _WIN32
	host.cpu=detect_windows_cpu_name();
	host.os=detect_windows_os_name();
#else
	host.cpu=detect_linux_cpu_name();
	host.os=detect_linux_os_name();
#endif
	if(host.cpu.empty()){
		host.cpu="Unknown CPU";
	}
	if(cpu_info.physical_cpus!=0||cpu_info.logical_cpus!=0){
		host.cpu+=" | "+std::to_string(cpu_info.physical_cpus)+"C/"+
				  std::to_string(cpu_info.logical_cpus)+"T";
	}
	if(host.os.empty()){
		host.os="Unknown OS";
	}
	return host;
}

void validate_self_test_options(const Options&opts){
	if(opts.has_to){
		throw std::invalid_argument("--stest cannot be combined with --to");
	}
	if(opts.from!=0){
		throw std::invalid_argument("--stest cannot be combined with --from");
	}
	if(opts.print_primes||opts.nth.has_value()||!opts.output_path.empty()||
	   opts.test_value.has_value()){
		throw std::invalid_argument(
			"--stest supports count benchmarking only");
	}
	if(opts.use_zstd||opts.output_format!=PrimeOutputFormat::Text){
		throw std::invalid_argument(
			"--stest does not support output-format or zstd options");
	}
	if(opts.use_ml){
		throw std::invalid_argument("--stest does not support --ml");
	}
}

struct TimedCountResult{
	std::uint64_t prime_count=0;
	std::uint64_t elapsed_us=0;
};

TimedCountResult run_timed_count_once(const Options&opts,const CpuInfo&info,
									  std::uint64_t upper_bound){
	if(upper_bound<=opts.from||upper_bound<2){
		throw std::invalid_argument("invalid range");
	}

	unsigned threads=opts.threads?opts.threads:effective_thread_count(info);
	if(threads==0){
		threads=1;
	}

	std::uint64_t odd_begin=opts.from<=3?3:opts.from;
	if((odd_begin&1ULL)==0){
		++odd_begin;
	}
	std::uint64_t odd_end=upper_bound;
	if((odd_end&1ULL)==0){
		++odd_end;
	}
	if(odd_begin>=upper_bound||odd_end<=odd_begin){
		odd_end=odd_begin;
	}

	SieveRange range{odd_begin,odd_end};
	std::uint64_t length=(range.end>range.begin)?(range.end-range.begin):0;
	SegmentConfig config=choose_segment_config(
		info,threads,opts.segment_bytes,opts.tile_bytes,length);

	const Wheel&wheel=get_wheel(opts.wheel);
	std::uint32_t small_limit=19u;
	switch(opts.wheel){
	case WheelType::Mod30:
		small_limit=19u;
		break;
	case WheelType::Mod210:
	case WheelType::Mod1155:
		small_limit=47u;
		break;
	}

	std::uint64_t sqrt_limit=static_cast<std::uint64_t>(std::sqrt(
									 static_cast<long double>(upper_bound)))+
								 1;
	auto base_primes=simple_sieve(sqrt_limit);
	auto start_time=std::chrono::steady_clock::now();

	bool can_wheel_bitmap=supports_wheel_bitmap_count(opts.wheel);
	std::uint64_t span=
		(upper_bound>opts.from)?(upper_bound-opts.from):0ULL;
	bool auto_wheel_bitmap=
		can_wheel_bitmap&&!opts.use_wheel_bitmap&&
		opts.segment_bytes==0&&opts.wheel==WheelType::Mod30&&
		((threads<=1&&span>=1000000000ULL)||
		 (threads>1&&span>=8000000000ULL));
	bool use_wheel_bitmap_path=can_wheel_bitmap&&
							   (opts.use_wheel_bitmap||auto_wheel_bitmap);

	if(use_wheel_bitmap_path){
		std::uint64_t total=count_primes_wheel_bitmap(
			opts.from,upper_bound,threads,opts.wheel,config,base_primes,
			opts.segment_bytes!=0);
		auto end_time=std::chrono::steady_clock::now();
		auto elapsed=std::chrono::duration_cast<std::chrono::microseconds>(
						 end_time-start_time)
						 .count();
		return TimedCountResult{
			total,static_cast<std::uint64_t>(std::max<std::int64_t>(elapsed,0))};
	}

	PrimeMarker marker(wheel,config,range.begin,range.end,base_primes,
					   small_limit);
	SegmentWorkQueue queue(range,config);

	std::vector<std::thread> workers;
	workers.reserve(threads);

	std::uint64_t prefix_count=0;
	if(opts.from<=2&&upper_bound>2){
		++prefix_count;
	}
	for(std::uint16_t p16 : wheel.presieved_primes){
		std::uint64_t p=static_cast<std::uint64_t>(p16);
		if(p>=opts.from&&p<upper_bound){
			++prefix_count;
		}
	}

	std::atomic<std::uint64_t> total{prefix_count};
	for(unsigned t=0;t<threads;++t){
		workers.emplace_back([&,t](){
			auto state=marker.make_thread_state(t,threads);
			std::vector<std::uint64_t> bitset;
			std::uint64_t local_total=0;
			while(true){
				std::uint64_t segment_id=0;
				std::uint64_t seg_low=0;
				std::uint64_t seg_high=0;
				if(!queue.next(segment_id,seg_low,seg_high)){
					break;
				}
				marker.sieve_segment(state,segment_id,seg_low,seg_high,bitset);
				std::size_t bit_count=
					static_cast<std::size_t>((seg_high-seg_low)>>1);
				local_total+=count_zero_bits(bitset.data(),bit_count);
			}
			total.fetch_add(local_total,std::memory_order_relaxed);
		});
	}
	for(auto&thread : workers){
		thread.join();
	}

	auto end_time=std::chrono::steady_clock::now();
	auto elapsed=std::chrono::duration_cast<std::chrono::microseconds>(
					 end_time-start_time)
					 .count();
	return TimedCountResult{
		total.load(std::memory_order_relaxed),
		static_cast<std::uint64_t>(std::max<std::int64_t>(elapsed,0))};
}

struct SelfTestRow{
	std::uint64_t upper_bound=0;
	int runs=0;
	int ok=0;
	double best_seconds=0.0;
	double median_seconds=0.0;
	double mean_seconds=0.0;
	std::uint64_t median_throughput=0;
};

int run_self_test(const Options&opts){
	validate_self_test_options(opts);

	CpuInfo info=detect_cpu_info();
	HostInfo host=detect_host_info(info);

	constexpr int runs_per_point=10;
	const std::array<std::uint64_t,6> points={
		1000000ULL,10000000ULL,100000000ULL,
		1000000000ULL,10000000000ULL,100000000000ULL};
	std::vector<SelfTestRow> rows;
	rows.reserve(points.size());

	for(std::uint64_t upper_bound : points){
		SelfTestRow row;
		row.upper_bound=upper_bound;
		row.runs=runs_per_point;
		std::vector<double> samples;
		samples.reserve(runs_per_point);
		for(int run=0;run<runs_per_point;++run){
			try{
				TimedCountResult result=run_timed_count_once(opts,info,
															 upper_bound);
				(void)result.prime_count;
				double seconds=
					static_cast<double>(result.elapsed_us)/1000000.0;
				samples.push_back(seconds);
				++row.ok;
			}catch(const std::exception&ex){
				std::cerr<<"[stest] --to "<<upper_bound<<" run "<<(run+1)
						 <<" failed: "<<ex.what()<<"\n";
			}
		}
		if(!samples.empty()){
			row.best_seconds=
				*std::min_element(samples.begin(),samples.end());
			row.median_seconds=median_seconds(samples);
			double sum=std::accumulate(samples.begin(),samples.end(),0.0);
			row.mean_seconds=sum/static_cast<double>(samples.size());
			if(row.median_seconds>0.0){
				long double throughput=
					static_cast<long double>(upper_bound)/
					static_cast<long double>(row.median_seconds);
				row.median_throughput=
					static_cast<std::uint64_t>(throughput+0.5L);
			}
		}
		rows.push_back(row);
	}

	std::cout<<"**Test data (10 runs/point):** (table from `"<<host.cpu
			 <<"` `"<<host.os<<"`)\n\n";
	std::cout<<"|     N (`--to`) | runs | ok | best (s) | median (s) | mean (s)"
			 <<" | throughput median (range/s) |\n";
	std::cout<<"| -------------: | ---: | -: | -------: | ---------: | -------:"
			 <<" | --------------------------: |\n";
	for(const SelfTestRow&row : rows){
		std::cout<<"| "<<std::setw(13)<<format_with_commas(row.upper_bound)
				 <<" | "<<std::setw(4)<<row.runs<<" | "<<std::setw(2)
				 <<row.ok<<" | "<<std::setw(10)
				 <<format_seconds(row.best_seconds)<<" | "<<std::setw(11)
				 <<format_seconds(row.median_seconds)<<" | "<<std::setw(10)
				 <<format_seconds(row.mean_seconds)<<" | "<<std::setw(26)
				 <<format_with_commas(row.median_throughput)<<" |\n";
	}

	return 0;
}

int run_cli(int argc,char**argv){
	try{
		Options opts=parse_options(argc,argv);
		if(opts.help){
			print_usage();
			return 0;
		}
		if(opts.self_test){
			return run_self_test(opts);
		}
#if !defined(CALCPRIME_HAS_ZSTD)
		if(opts.use_zstd){
			throw std::invalid_argument("zstd not supported in this build");
		}
#endif
		if(opts.test_value.has_value()&&!opts.has_to){
			bool is_prime=miller_rabin_is_prime(opts.test_value.value());
			std::cout<<(is_prime?"prime":"composite")<<"\n";
			return 0;
		}
		if(!opts.has_to){
			print_usage();
			return 1;
		}
		if(opts.test_value.has_value()){
			bool is_prime=miller_rabin_is_prime(opts.test_value.value());
			std::cout<<(is_prime?"prime":"composite")<<"\n";
		}
		if(opts.to<=opts.from||opts.to<2){
			throw std::invalid_argument("invalid range");
		}

		CpuInfo info=detect_cpu_info();
		unsigned threads=opts.threads?opts.threads:effective_thread_count(info);
		if(opts.nth.has_value()){
			threads=1; // ensure deterministic ordering
		}
		if(threads==0){
			threads=1;
		}

		std::uint64_t odd_begin=opts.from<=3?3:opts.from;
		if((odd_begin&1ULL)==0){
			++odd_begin;
		}
		std::uint64_t odd_end=opts.to;
		if((odd_end&1ULL)==0){
			++odd_end;
		}
		if(odd_begin>=opts.to||odd_end<=odd_begin){
			odd_end=odd_begin;
		}

		SieveRange range{odd_begin,odd_end};
		std::uint64_t length=(range.end>range.begin)?(range.end-range.begin):0;

		SegmentConfig config=choose_segment_config(
			info,threads,opts.segment_bytes,opts.tile_bytes,length);
		const Wheel&wheel=get_wheel(opts.wheel);
		std::uint32_t small_limit=19u;
		switch(opts.wheel){
		case WheelType::Mod30:
			small_limit=19u;
			break;
		case WheelType::Mod210:
			small_limit=47u;
			break;
		case WheelType::Mod1155:
			small_limit=47u;
			break;
		}
		std::size_t num_segments=
			length?static_cast<std::size_t>((length+config.segment_span-1)/
											config.segment_span)
				  :0;

		std::uint64_t sqrt_limit=static_cast<std::uint64_t>(std::sqrt(
									 static_cast<long double>(opts.to)))+
								 1;
		auto base_primes=simple_sieve(sqrt_limit);

		bool is_count_mode=
			opts.count_only||(!opts.print_primes&&!opts.nth.has_value());
		auto start_time=std::chrono::steady_clock::now();
		if(opts.use_ml){
			if(opts.print_primes||opts.nth.has_value()){
				throw std::invalid_argument("--ml supports count mode only");
			}
			std::uint64_t total=
				meissel_count(opts.from,opts.to,base_primes,threads);
			auto end_time=std::chrono::steady_clock::now();
			std::cout<<total<<"\n";

			if(opts.show_stats){
				std::cout<<"Threads: "<<threads<<"\n";
				std::cout<<"Segment bytes: "<<config.segment_bytes<<"\n";
				std::cout<<"Tile bytes: "<<config.tile_bytes<<"\n";
				std::cout<<"L1d: "<<info.l1_data_bytes<<"  L2: "<<info.l2_bytes
						 <<"\n";
			}

			if(opts.show_time){
				auto elapsed=
					std::chrono::duration_cast<std::chrono::microseconds>(
						end_time-start_time)
						.count();
				std::cout<<"Elapsed: "<<elapsed<<" us\n";
			}
			return 0;
		}

		bool can_wheel_bitmap=is_count_mode&&!opts.print_primes&&
							 !opts.nth.has_value()&&
							 supports_wheel_bitmap_count(opts.wheel);
		std::uint64_t span=
			(opts.to>opts.from)?(opts.to-opts.from):0ULL;
		bool auto_wheel_bitmap=
			can_wheel_bitmap&&!opts.use_wheel_bitmap&&
			opts.segment_bytes==0&&opts.wheel==WheelType::Mod30&&
			((threads<=1&&span>=1000000000ULL)||
			 (threads>1&&span>=8000000000ULL));
		bool use_wheel_bitmap_path=can_wheel_bitmap&&
								   (opts.use_wheel_bitmap||auto_wheel_bitmap);

		if(use_wheel_bitmap_path){
			std::uint64_t total=count_primes_wheel_bitmap(
				opts.from,opts.to,threads,opts.wheel,config,base_primes,
				opts.segment_bytes!=0);
			auto end_time=std::chrono::steady_clock::now();
			std::cout<<total<<"\n";

			if(opts.show_stats){
				std::cout<<"Wheel bitmap: "
						 <<(opts.use_wheel_bitmap?"forced":"auto")<<"\n";
				std::cout<<"Threads: "<<threads<<"\n";
				std::cout<<"Segment bytes: "<<config.segment_bytes<<"\n";
				std::cout<<"Tile bytes: "<<config.tile_bytes<<"\n";
				std::cout<<"L1d: "<<info.l1_data_bytes<<"  L2: "<<info.l2_bytes
						 <<"\n";
			}

			if(opts.show_time){
				auto elapsed=
					std::chrono::duration_cast<std::chrono::microseconds>(
						end_time-start_time)
						.count();
				std::cout<<"Elapsed: "<<elapsed<<" us\n";
			}
			return 0;
		}

		PrimeMarker marker(wheel,config,range.begin,range.end,base_primes,
						   small_limit);
		SegmentWorkQueue queue(range,config);

		std::vector<SegmentResult> segment_results(num_segments);
		std::mutex segment_ready_mutex;
		std::condition_variable segment_ready_cv;
		std::atomic<bool> stop{false};
		std::atomic<bool> nth_found{false};
		std::uint64_t nth_value=0;
		std::uint64_t nth_target=opts.nth.value_or(0);

		std::vector<std::thread> workers;
		workers.reserve(threads);

		bool include_two=opts.from<=2&&opts.to>2;
		std::vector<std::uint64_t> prefix_primes;
		if(include_two){
			prefix_primes.push_back(2);
		}
		for(std::uint16_t p16 : wheel.presieved_primes){
			std::uint64_t p=static_cast<std::uint64_t>(p16);
			if(p>=opts.from&&p<opts.to){
				prefix_primes.push_back(p);
			}
		}
		std::uint64_t prefix_count=prefix_primes.size();

		if(opts.nth.has_value()&&opts.nth.value()<=prefix_count){
			std::cout<<prefix_primes[opts.nth.value()-1]<<"\n";
			return 0;
		}

		if(is_count_mode&&!opts.print_primes&&!opts.nth.has_value()){
			std::atomic<std::uint64_t> total{prefix_count};
			for(unsigned t=0;t<threads;++t){
				workers.emplace_back([&,t](){
					auto state=marker.make_thread_state(t,threads);
					std::vector<std::uint64_t> bitset;
					std::uint64_t local_total=0;
					while(true){
						std::uint64_t segment_id=0;
						std::uint64_t seg_low=0;
						std::uint64_t seg_high=0;
						if(!queue.next(segment_id,seg_low,seg_high)){
							break;
						}
						marker.sieve_segment(state,segment_id,seg_low,seg_high,
											 bitset);
						std::size_t bit_count=
							static_cast<std::size_t>((seg_high-seg_low)>>1);
						local_total+=count_zero_bits(bitset.data(),bit_count);
					}
					total.fetch_add(local_total,std::memory_order_relaxed);
				});
			}
			for(auto&th : workers){
				th.join();
			}

			auto end_time=std::chrono::steady_clock::now();
			std::cout<<total.load(std::memory_order_relaxed)<<"\n";

			if(opts.show_stats){
				std::cout<<"Threads: "<<threads<<"\n";
				std::cout<<"Segment bytes: "<<config.segment_bytes<<"\n";
				std::cout<<"Tile bytes: "<<config.tile_bytes<<"\n";
				std::cout<<"L1d: "<<info.l1_data_bytes<<"  L2: "<<info.l2_bytes
						 <<"\n";
			}

			if(opts.show_time){
				auto elapsed=
					std::chrono::duration_cast<std::chrono::microseconds>(
						end_time-start_time)
						.count();
				std::cout<<"Elapsed: "<<elapsed<<" us\n";
			}

			return 0;
		}

		PrimeWriter writer(opts.print_primes,opts.output_path,
						   opts.output_format,opts.use_zstd);
		std::mutex writer_exception_mutex;
		std::exception_ptr writer_exception;
		std::thread writer_feeder;

		for(unsigned t=0;t<threads;++t){
			workers.emplace_back([&,t](){
				auto state=marker.make_thread_state(t,threads);
				std::vector<std::uint64_t> bitset;
				std::uint64_t cumulative=prefix_count;
				while(!stop.load(std::memory_order_relaxed)){
					std::uint64_t segment_id=0;
					std::uint64_t seg_low=0;
					std::uint64_t seg_high=0;
					if(!queue.next(segment_id,seg_low,seg_high)){
						break;
					}
					marker.sieve_segment(state,segment_id,seg_low,seg_high,
										 bitset);
					std::size_t bit_count=
						static_cast<std::size_t>((seg_high-seg_low)>>1);
					std::uint64_t local_count=
						count_zero_bits(bitset.data(),bit_count);
					if(segment_id<segment_results.size()){
						segment_results[segment_id].count=local_count;
					}
					bool need_primes=
						opts.print_primes||(opts.nth.has_value()&&threads==1);
					if(need_primes&&segment_id<segment_results.size()){
						std::vector<std::uint64_t> primes;
						primes.reserve(static_cast<std::size_t>(local_count));
						extract_segment_primes(bitset,seg_low,bit_count,primes);
						segment_results[segment_id].primes=std::move(primes);
						{
							std::lock_guard<std::mutex> lock(
								segment_ready_mutex);
							segment_results[segment_id].ready.store(
								true,std::memory_order_release);
						}
						segment_ready_cv.notify_one();
						if(opts.nth.has_value()&&threads==1&&
						   !nth_found.load(std::memory_order_relaxed)){
							std::uint64_t base=cumulative;
							std::uint64_t new_total=base+local_count;
							if(nth_target>base&&nth_target<=new_total){
								std::size_t index=
									static_cast<std::size_t>(nth_target-base-1);
								if(index<
								   segment_results[segment_id].primes.size()){
									nth_value=segment_results[segment_id]
												  .primes[index];
									nth_found.store(true,
													std::memory_order_relaxed);
									stop.store(true,std::memory_order_relaxed);
								}
							}
							cumulative=new_total;
						}
					}else{
						if(opts.nth.has_value()&&threads==1&&
						   !nth_found.load(std::memory_order_relaxed)){
							std::uint64_t base=cumulative;
							std::uint64_t new_total=base+local_count;
							if(nth_target>base&&nth_target<=new_total){
								// Need to extract primes to find nth
								std::vector<std::uint64_t> primes;
								primes.reserve(
									static_cast<std::size_t>(local_count));
								extract_segment_primes(bitset,seg_low,bit_count,
													   primes);
								std::size_t index=
									static_cast<std::size_t>(nth_target-base-1);
								if(index<primes.size()){
									nth_value=primes[index];
									nth_found.store(true,
													std::memory_order_relaxed);
									stop.store(true,std::memory_order_relaxed);
								}
								if(segment_id<segment_results.size()){
									segment_results[segment_id].primes=
										std::move(primes);
									{
										std::lock_guard<std::mutex> lock(
											segment_ready_mutex);
										segment_results[segment_id].ready.store(
											true,std::memory_order_release);
									}
									segment_ready_cv.notify_one();
								}
							}
							cumulative=new_total;
						}
					}
				}
			});
		}

		if(opts.print_primes){
			std::vector<std::uint64_t> prefix_copy=prefix_primes;
			writer_feeder=std::thread([&,prefix_copy]() mutable{
				try{
					if(!prefix_copy.empty()){
						writer.write_segment(prefix_copy);
					}
					for(std::size_t next=0;next<segment_results.size();++next){
						SegmentResult&res=segment_results[next];
						std::unique_lock<std::mutex> lock(segment_ready_mutex);
						segment_ready_cv.wait(lock,[&]{
							return res.ready.load(std::memory_order_acquire)||
								   stop.load(std::memory_order_relaxed);
						});
						bool ready=res.ready.load(std::memory_order_acquire);
						if(!ready){
							break;
						}
						res.ready.store(false,std::memory_order_relaxed);
						std::vector<std::uint64_t> primes=std::move(res.primes);
						lock.unlock();
						writer.write_segment(primes);
					}
					writer.flush();
				}catch(...){
					std::lock_guard<std::mutex> err_lock(
						writer_exception_mutex);
					writer_exception=std::current_exception();
					stop.store(true,std::memory_order_relaxed);
					segment_ready_cv.notify_all();
				}
			});
		}

		for(auto&th : workers){
			th.join();
		}

		segment_ready_cv.notify_all();

		auto end_time=std::chrono::steady_clock::now();

		std::uint64_t total=prefix_count;
		for(const auto&res : segment_results){
			total+=res.count;
		}

		if(is_count_mode){
			std::cout<<total<<"\n";
		}

		if(writer_feeder.joinable()){
			writer_feeder.join();
		}

		std::exception_ptr pending_writer_exception;
		{
			std::lock_guard<std::mutex> lock(writer_exception_mutex);
			pending_writer_exception=writer_exception;
		}

		try{
			writer.finish();
		}catch(...){
			if(!pending_writer_exception){
				pending_writer_exception=std::current_exception();
			}
		}

		if(pending_writer_exception){
			std::rethrow_exception(pending_writer_exception);
		}

		if(opts.nth.has_value()){
			if(!nth_found.load()){
				std::cerr<<"nth prime not found within range\n";
				return 1;
			}
			std::cout<<nth_value<<"\n";
		}

		if(opts.show_stats){
			std::cout<<"Threads: "<<threads<<"\n";
			std::cout<<"Segment bytes: "<<config.segment_bytes<<"\n";
			std::cout<<"Tile bytes: "<<config.tile_bytes<<"\n";
			std::cout<<"L1d: "<<info.l1_data_bytes<<"  L2: "<<info.l2_bytes
					 <<"\n";
		}

		if(opts.show_time){
			auto elapsed=std::chrono::duration_cast<std::chrono::microseconds>(
							 end_time-start_time)
							 .count();
			std::cout<<"Elapsed: "<<elapsed<<" us\n";
		}

		return 0;
	}catch(const std::exception&ex){
		std::cerr<<"Error: "<<ex.what()<<"\n";
		return 1;
	}
}

} // namespace calcprime

#ifndef CALCPRIME_DLL_BUILD
int main(int argc,char**argv){ return calcprime::run_cli(argc,argv); }
#endif










