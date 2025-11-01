#include "base_sieve.h"
#include "cpu_info.h"
#include "marker.h"
#include "popcnt.h"
#include "prime_count.h"
#include "segmenter.h"
#include "wheel.h"
#include "writer.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cinttypes>
#include <cstddef>
#include <cstdint>
#include <cmath>
#include <condition_variable>
#include <exception>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <optional>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>
#include <limits>

namespace calcprime {

struct Options {
    std::uint64_t from=0;
    std::uint64_t to=0;
    bool has_to=false;
    bool count_only=true;
    bool print_primes=false;
    std::optional<std::uint64_t>nth;
    unsigned threads=0;
    WheelType wheel=WheelType::Mod30;
    std::size_t segment_bytes=0;
    std::size_t tile_bytes=0;
    std::string output_path;
    PrimeOutputFormat output_format=PrimeOutputFormat::Text;
    bool show_time=false;
    bool show_stats=false;
    bool use_ml=false;
    bool help=false;
    std::optional<std::uint64_t>test_value;
};

std::uint64_t parse_u64(const std::string&value) {
    if(value.empty()) {
        throw std::invalid_argument("invalid integer: "+value);
    }

    auto exp_pos=value.find_first_of("eE");
    if(exp_pos!=std::string::npos) {
        std::string mantissa_str=value.substr(0,exp_pos);
        std::string exponent_str=value.substr(exp_pos+1);
        if(mantissa_str.empty()||exponent_str.empty()) {
            throw std::invalid_argument("invalid integer: "+value);
        }

        std::size_t mantissa_idx=0;
        std::uint64_t mantissa=0;
        try {
            mantissa=std::stoull(mantissa_str,&mantissa_idx,0);
        } catch(const std::exception&) {
            throw std::invalid_argument("invalid integer: "+value);
        }
        if(mantissa_idx!=mantissa_str.size()) {
            throw std::invalid_argument("invalid integer: "+value);
        }

        std::size_t exponent_idx=0;
        long long exponent=0;
        try {
            exponent=std::stoll(exponent_str,&exponent_idx,10);
        } catch(const std::exception&) {
            throw std::invalid_argument("invalid integer: "+value);
        }
        if(exponent_idx!=exponent_str.size()) {
            throw std::invalid_argument("invalid integer: "+value);
        }
        if(exponent<0) {
            throw std::invalid_argument("invalid integer: "+value);
        }

        std::uint64_t result=mantissa;
        for(long long i=0;i<exponent;++i) {
            if(result>std::numeric_limits<std::uint64_t>::max()/10ULL) {
                throw std::invalid_argument("integer too large: "+value);
            }
            result*=10ULL;
        }
        return result;
    }

    std::size_t idx=0;
    std::uint64_t result=0;
    try {
        result=std::stoull(value,&idx,0);
    } catch(const std::exception&) {
        throw std::invalid_argument("invalid integer: "+value);
    }
    if(idx!=value.size()) {
        throw std::invalid_argument("invalid integer: "+value);
    }
    return result;
}

std::size_t parse_size(const std::string&value) {
    if(value.empty()) {
        throw std::invalid_argument("invalid size");
    }
    std::size_t idx=0;
    std::uint64_t base=std::stoull(value,&idx,0);
    std::size_t factor=1;
    if(idx<value.size()) {
        char suffix=value[idx];
        switch(suffix) {
        case'k':
        case'K':
            factor=1024;
            ++idx;
            break;
        case'm':
        case'M':
            factor=1024*1024;
            ++idx;
            break;
        case'g':
        case'G':
            factor=1024ull*1024ull*1024ull;
            ++idx;
            break;
        default:
            break;
        }
    }
    if(idx!=value.size()) {
        throw std::invalid_argument("invalid size suffix: "+value);
    }
    std::uint64_t result=base*factor;
    if(result>static_cast<std::uint64_t>(std::numeric_limits<std::size_t>::max())) {
        throw std::invalid_argument("size too large: "+value);
    }
    return static_cast<std::size_t>(result);
}

Options parse_options(int argc,char**argv) {
    Options opts;
    for(int i=1;i<argc;++i) {
        std::string arg=argv[i];
        static const std::string out_format_prefix="--out-format=";

        if(arg=="--help"||arg=="-h") {
            opts.help=true;
            return opts;
        } else if(arg.rfind(out_format_prefix,0)==0) {
            std::string fmt=arg.substr(out_format_prefix.size());
            if(fmt=="text") {
                opts.output_format=PrimeOutputFormat::Text;
            } else if(fmt=="binary") {
                opts.output_format=PrimeOutputFormat::Binary;
            } else if(fmt=="zstd"||fmt=="zstd+delta") {
                opts.output_format=PrimeOutputFormat::ZstdDelta;
            } else {
                throw std::invalid_argument("unsupported out-format: "+fmt);
            }
        } else if(arg=="--from") {
            if(i+1>=argc) {
                throw std::invalid_argument("--from requires a value");
            }
            opts.from=parse_u64(argv[++i]);
        } else if(arg=="--to") {
            if(i+1>=argc) {
                throw std::invalid_argument("--to requires a value");
            }
            opts.to=parse_u64(argv[++i]);
            opts.has_to=true;
        } else if(arg=="--count") {
            opts.count_only=true;
        } else if(arg=="--print") {
            opts.print_primes=true;
            opts.count_only=false;
        } else if(arg=="--nth") {
            if(i+1>=argc) {
                throw std::invalid_argument("--nth requires a value");
            }
            opts.nth=parse_u64(argv[++i]);
            opts.count_only=false;
        } else if(arg=="--threads") {
            if(i+1>=argc) {
                throw std::invalid_argument("--threads requires a value");
            }
            opts.threads=static_cast<unsigned>(parse_u64(argv[++i]));
        } else if(arg=="--wheel") {
            if(i+1>=argc) {
                throw std::invalid_argument("--wheel requires a value");
            }
            std::string w=argv[++i];
            if(w=="30") {
                opts.wheel=WheelType::Mod30;
            } else if(w=="210") {
                opts.wheel=WheelType::Mod210;
            } else if(w=="1155") {
                opts.wheel=WheelType::Mod1155;
            } else {
                throw std::invalid_argument("unsupported wheel: "+w);
            }
        } else if(arg=="--segment") {
            if(i+1>=argc) {
                throw std::invalid_argument("--segment requires a value");
            }
            opts.segment_bytes=parse_size(argv[++i]);
        } else if(arg=="--tile") {
            if(i+1>=argc) {
                throw std::invalid_argument("--tile requires a value");
            }
            opts.tile_bytes=parse_size(argv[++i]);
        } else if(arg=="--out") {
            if(i+1>=argc) {
                throw std::invalid_argument("--out requires a path");
            }
            opts.output_path=argv[++i];
        } else if(arg=="--out-format") {
            if(i+1>=argc) {
                throw std::invalid_argument("--out-format requires a value");
            }
            std::string fmt=argv[++i];
            if(fmt=="text") {
                opts.output_format=PrimeOutputFormat::Text;
            } else if(fmt=="binary") {
                opts.output_format=PrimeOutputFormat::Binary;
            } else if(fmt=="zstd"||fmt=="zstd+delta") {
                opts.output_format=PrimeOutputFormat::ZstdDelta;
            } else {
                throw std::invalid_argument("unsupported out-format: "+fmt);
            }
        } else if(arg=="--time") {
            opts.show_time=true;
        } else if(arg=="--stats") {
            opts.show_stats=true;
        } else if(arg=="--ml") {
            opts.use_ml=true;
        } else if(arg=="--test") {
            if(i+1>=argc) {
                throw std::invalid_argument("--test requires a value");
            }
            opts.test_value=parse_u64(argv[++i]);
        } else {
            throw std::invalid_argument("unknown option: "+arg);
        }
    }
    return opts;
}

void print_usage() {
    std::cout<<"prime-sieve --from A --to B [options]\n"
              <<"  --count             Count primes (default)\n"
              <<"  --print             Print primes in the interval\n"
              <<"  --nth K             Find the K-th prime in the interval\n"
              <<"  --threads N         Override thread count\n"
              <<"  --wheel 30|210|1155 Select wheel factorisation (default 30)\n"
              <<"  --segment BYTES     Override segment size\n"
              <<"  --tile BYTES        Override tile size\n"
              <<"  --out PATH          Write primes to file\n"
              <<"  --out-format FMT    Output format: text (default), binary, zstd\n"
              <<"  --time              Print elapsed time\n"
              <<"  --stats             Print configuration statistics\n"
              <<"  --ml                Use Meissel-Lehmer counting for --count\n"
              <<"  --test N           Run a Miller-Rabin primality check for N\n";
}

struct SegmentResult {
    std::uint64_t count=0;
    std::vector<std::uint64_t>primes;
    std::atomic<bool>ready{false};
};

int run_cli(int argc,char**argv) {
    try {
        Options opts=parse_options(argc,argv);
        if(opts.help) {
            print_usage();
            return 0;
        }
        if(opts.test_value.has_value()&&!opts.has_to) {
            bool is_prime=miller_rabin_is_prime(opts.test_value.value());
            std::cout<<(is_prime ?"prime" :"composite")<<"\n";
            return 0;
        }
        if(!opts.has_to) {
            print_usage();
            return 1;
        }
        if(opts.test_value.has_value()) {
            bool is_prime=miller_rabin_is_prime(opts.test_value.value());
            std::cout<<(is_prime ?"prime" :"composite")<<"\n";
        }
        if(opts.to<=opts.from||opts.to<2) {
            throw std::invalid_argument("invalid range");
        }

        CpuInfo info=detect_cpu_info();
        unsigned threads=opts.threads?opts.threads:effective_thread_count(info);
        if(opts.nth.has_value()) {
            threads=1;
        }
        if(threads==0) {
            threads=1;
        }

        std::uint64_t odd_begin=opts.from<=3?3:opts.from;
        if((odd_begin&1ULL)==0) {
            ++odd_begin;
        }
        if(odd_begin>=opts.to) {
            odd_begin=3;
        }
        std::uint64_t odd_end=opts.to;
        if((odd_end&1ULL)==0) {
            ++odd_end;
        }
        if(odd_end<=odd_begin) {
            odd_end=odd_begin;
        }

        SieveRange range{odd_begin,odd_end};
        std::uint64_t length=(range.end>range.begin)?(range.end-range.begin):0;

        SegmentConfig config=
            choose_segment_config(info,threads,opts.segment_bytes,opts.tile_bytes,length);
        const Wheel&wheel=get_wheel(opts.wheel);
        std::uint32_t small_limit=29u;
        switch(opts.wheel) {
        case WheelType::Mod30:
            small_limit=29u;
            break;
        case WheelType::Mod210:
            small_limit=47u;
            break;
        case WheelType::Mod1155:
            small_limit=47u;
            break;
        }
        std::size_t num_segments=length?static_cast<std::size_t>((length+config.segment_span-1)/config.segment_span):0;

        std::uint64_t sqrt_limit=static_cast<std::uint64_t>(std::sqrt(static_cast<long double>(opts.to)))+1;
        auto base_primes=simple_sieve(sqrt_limit);

        bool is_count_mode=opts.count_only||(!opts.print_primes&&!opts.nth.has_value());
        auto start_time=std::chrono::steady_clock::now();

        if(opts.use_ml&&is_count_mode) {
            std::uint64_t result=meissel_count(opts.from,opts.to,base_primes,threads);
            auto end_time=std::chrono::steady_clock::now();

            std::cout<<result<<"\n";

            if(opts.show_stats) {
                unsigned ml_threads=threads==0 ? 1u : threads;
                std::cout<<"Threads: "<<ml_threads<<"\n";
                std::cout<<"Segment bytes: 0\n";
                std::cout<<"Tile bytes: 0\n";
                std::cout<<"L1d: "<<info.l1_data_bytes<<"  L2: "<<info.l2_bytes<<"\n";
            }

            if(opts.show_time) {
                auto elapsed=std::chrono::duration_cast<std::chrono::microseconds>(end_time-start_time).count();
                std::cout<<"Elapsed: "<<elapsed<<" us\n";
            }

            return 0;
        }

        PrimeMarker marker(wheel,config,range.begin,range.end,base_primes,small_limit);
        SegmentWorkQueue queue(range,config);

        std::vector<SegmentResult>segment_results(num_segments);
        std::mutex segment_ready_mutex;
        std::condition_variable segment_ready_cv;
        std::atomic<bool>stop{false};
        std::atomic<bool>nth_found{false};
        std::uint64_t nth_value=0;
        std::uint64_t nth_target=opts.nth.value_or(0);

        std::vector<std::thread>workers;
        workers.reserve(threads);

        bool include_two=opts.from<=2&&opts.to>2;
        std::vector<std::uint64_t>prefix_primes;
        if(include_two) {
            prefix_primes.push_back(2);
        }
        std::vector<std::uint64_t>wheel_primes;
        switch(opts.wheel) {
        case WheelType::Mod30:
            wheel_primes={3,5};
            break;
        case WheelType::Mod210:
            wheel_primes={3,5,7};
            break;
        case WheelType::Mod1155:
            wheel_primes={3,5,7,11};
            break;
        }
        for(std::uint64_t p : wheel_primes) {
            if(p>=opts.from&&p<opts.to) {
                prefix_primes.push_back(p);
            }
        }
        std::uint64_t prefix_count=prefix_primes.size();

        if(opts.nth.has_value()&&opts.nth.value()<=prefix_count) {
            std::cout<<prefix_primes[opts.nth.value()-1]<<"\n";
            return 0;
        }

        PrimeWriter writer(opts.print_primes,opts.output_path,opts.output_format);
        std::mutex writer_exception_mutex;
        std::exception_ptr writer_exception;
        std::thread writer_feeder;

        for(unsigned t=0;t<threads;++t) {
            workers.emplace_back([&,t]() {
                auto state=marker.make_thread_state(t,threads);
                std::vector<std::uint64_t>bitset;
                std::uint64_t cumulative=prefix_count;
                while(!stop.load(std::memory_order_relaxed)) {
                    std::uint64_t segment_id=0;
                    std::uint64_t seg_low=0;
                    std::uint64_t seg_high=0;
                    if(!queue.next(segment_id,seg_low,seg_high)) {
                        break;
                    }
                    marker.sieve_segment(state,segment_id,seg_low,seg_high,bitset);
                    std::size_t bit_count=static_cast<std::size_t>((seg_high-seg_low)>>1);
                    std::uint64_t local_count=count_zero_bits(bitset.data(),bit_count);
                    if(segment_id<segment_results.size()) {
                        segment_results[segment_id].count=local_count;
                    }
                    bool need_primes=opts.print_primes||(opts.nth.has_value()&&threads==1);
                    if(need_primes&&segment_id<segment_results.size()) {
                        std::vector<std::uint64_t>primes;
                        primes.reserve(static_cast<std::size_t>(local_count));
                        std::uint64_t value=seg_low;
                        std::size_t produced=0;
                        for(std::size_t word=0;word<bitset.size()&&produced<bit_count;++word) {
                            std::uint64_t composite=bitset[word];
                            for(std::size_t bit=0;bit<64&&produced<bit_count;++bit,++produced,value+=2) {
                                if(composite&(1ULL<<bit)) {
                                    continue;
                                }
                                primes.push_back(value);
                            }
                        }
                        segment_results[segment_id].primes=std::move(primes);
                        {
                            std::lock_guard<std::mutex>lock(segment_ready_mutex);
                            segment_results[segment_id].ready.store(true,std::memory_order_release);
                        }
                        segment_ready_cv.notify_all();
                        if(opts.nth.has_value()&&threads==1&&!nth_found.load(std::memory_order_relaxed)) {
                            std::uint64_t base=cumulative;
                            std::uint64_t new_total=base+local_count;
                            if(nth_target>base&&nth_target<=new_total) {
                                std::size_t index=static_cast<std::size_t>(nth_target-base-1);
                                if(index<segment_results[segment_id].primes.size()) {
                                    nth_value=segment_results[segment_id].primes[index];
                                    nth_found.store(true,std::memory_order_relaxed);
                                    stop.store(true,std::memory_order_relaxed);
                                }
                            }
                            cumulative=new_total;
                        }
                    } else {
                        if(opts.nth.has_value()&&threads==1&&!nth_found.load(std::memory_order_relaxed)) {
                            std::uint64_t base=cumulative;
                            std::uint64_t new_total=base+local_count;
                            if(nth_target>base&&nth_target<=new_total) {

                                std::vector<std::uint64_t>primes;
                                primes.reserve(static_cast<std::size_t>(local_count));
                                std::uint64_t value=seg_low;
                                std::size_t produced=0;
                                for(std::size_t word=0;word<bitset.size()&&produced<bit_count;++word) {
                                    std::uint64_t composite=bitset[word];
                                    for(std::size_t bit=0;bit<64&&produced<bit_count;++bit,++produced,value+=2) {
                                        if(composite&(1ULL<<bit)) {
                                            continue;
                                        }
                                        primes.push_back(value);
                                    }
                                }
                                std::size_t index=static_cast<std::size_t>(nth_target-base-1);
                                if(index<primes.size()) {
                                    nth_value=primes[index];
                                    nth_found.store(true,std::memory_order_relaxed);
                                    stop.store(true,std::memory_order_relaxed);
                                }
                                if(segment_id<segment_results.size()) {
                                    segment_results[segment_id].primes=std::move(primes);
                                    {
                                        std::lock_guard<std::mutex>lock(segment_ready_mutex);
                                        segment_results[segment_id].ready.store(true,std::memory_order_release);
                                    }
                                    segment_ready_cv.notify_all();
                                }
                            }
                            cumulative=new_total;
                        }
                    }
                }
            });
        }

        if(opts.print_primes) {
            std::vector<std::uint64_t>prefix_copy=prefix_primes;
            writer_feeder=std::thread([&,prefix_copy]() mutable {
                try {
                    if(!prefix_copy.empty()) {
                        writer.write_segment(prefix_copy);
                    }
                    for(std::size_t next=0;next<segment_results.size();++next) {
                        SegmentResult&res=segment_results[next];
                        std::unique_lock<std::mutex>lock(segment_ready_mutex);
                        segment_ready_cv.wait(lock,[&] {
                            return res.ready.load(std::memory_order_acquire)||
                                   stop.load(std::memory_order_relaxed);
                        });
                        bool ready=res.ready.load(std::memory_order_acquire);
                        if(!ready) {
                            break;
                        }
                        res.ready.store(false,std::memory_order_relaxed);
                        std::vector<std::uint64_t>primes=std::move(res.primes);
                        lock.unlock();
                        writer.write_segment(primes);
                    }
                    writer.flush();
                } catch(...) {
                    std::lock_guard<std::mutex>err_lock(writer_exception_mutex);
                    writer_exception=std::current_exception();
                }
            });
        }

        for(auto&th : workers) {
            th.join();
        }

        segment_ready_cv.notify_all();

        auto end_time=std::chrono::steady_clock::now();

        std::uint64_t total=prefix_count;
        for(const auto&res : segment_results) {
            total+=res.count;
        }

        if(is_count_mode) {
            std::cout<<total<<"\n";
        }

        if(writer_feeder.joinable()) {
            writer_feeder.join();
        }

        std::exception_ptr pending_writer_exception;
        {
            std::lock_guard<std::mutex>lock(writer_exception_mutex);
            pending_writer_exception=writer_exception;
        }

        try {
            writer.finish();
        } catch(...) {
            if(!pending_writer_exception) {
                pending_writer_exception=std::current_exception();
            }
        }

        if(pending_writer_exception) {
            std::rethrow_exception(pending_writer_exception);
        }

        if(opts.nth.has_value()) {
            if(!nth_found.load()) {
                std::cerr<<"nth prime not found within range\n";
                return 1;
            }
            std::cout<<nth_value<<"\n";
        }

        if(opts.show_stats) {
            std::cout<<"Threads: "<<threads<<"\n";
            std::cout<<"Segment bytes: "<<config.segment_bytes<<"\n";
            std::cout<<"Tile bytes: "<<config.tile_bytes<<"\n";
            std::cout<<"L1d: "<<info.l1_data_bytes<<"  L2: "<<info.l2_bytes<<"\n";
        }

        if(opts.show_time) {
            auto elapsed=std::chrono::duration_cast<std::chrono::microseconds>(end_time-start_time).count();
            std::cout<<"Elapsed: "<<elapsed<<" us\n";
        }

        return 0;
    } catch(const std::exception&ex) {
        std::cerr<<"Error: "<<ex.what()<<"\n";
        return 1;
    }
}

}

#ifndef CALCPRIME_DLL_BUILD
int main(int argc,char**argv) {
    return calcprime::run_cli(argc,argv);
}
#endif
