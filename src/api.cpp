#include "calcprime/api.h"

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
#include <cmath>
#include <condition_variable>
#include <exception>
#include <memory>
#include <mutex>
#include <new>
#include <string>
#include <thread>
#include <vector>

namespace {

bool is_valid_wheel(calcprime_wheel_type type) {
    switch(type) {
    case CALCPRIME_WHEEL_MOD30:
    case CALCPRIME_WHEEL_MOD210:
    case CALCPRIME_WHEEL_MOD1155:
        return true;
    }
    return false;
}

bool is_valid_output_format(calcprime_output_format format) {
    switch(format) {
    case CALCPRIME_OUTPUT_TEXT:
    case CALCPRIME_OUTPUT_BINARY:
    case CALCPRIME_OUTPUT_ZSTD_DELTA:
        return true;
    }
    return false;
}

calcprime::WheelType to_cpp_wheel(calcprime_wheel_type type) {
    switch(type) {
    case CALCPRIME_WHEEL_MOD30:
        return calcprime::WheelType::Mod30;
    case CALCPRIME_WHEEL_MOD210:
        return calcprime::WheelType::Mod210;
    case CALCPRIME_WHEEL_MOD1155:
        return calcprime::WheelType::Mod1155;
    }
    return calcprime::WheelType::Mod30;
}

calcprime::PrimeOutputFormat to_cpp_output(calcprime_output_format format) {
    switch(format) {
    case CALCPRIME_OUTPUT_TEXT:
        return calcprime::PrimeOutputFormat::Text;
    case CALCPRIME_OUTPUT_BINARY:
        return calcprime::PrimeOutputFormat::Binary;
    case CALCPRIME_OUTPUT_ZSTD_DELTA:
        return calcprime::PrimeOutputFormat::ZstdDelta;
    }
    return calcprime::PrimeOutputFormat::Text;
}

calcprime_wheel_type to_c_wheel(calcprime::WheelType type) {
    switch(type) {
    case calcprime::WheelType::Mod30:
        return CALCPRIME_WHEEL_MOD30;
    case calcprime::WheelType::Mod210:
        return CALCPRIME_WHEEL_MOD210;
    case calcprime::WheelType::Mod1155:
        return CALCPRIME_WHEEL_MOD1155;
    }
    return CALCPRIME_WHEEL_MOD30;
}

calcprime_output_format to_c_output(calcprime::PrimeOutputFormat format) {
    switch(format) {
    case calcprime::PrimeOutputFormat::Text:
        return CALCPRIME_OUTPUT_TEXT;
    case calcprime::PrimeOutputFormat::Binary:
        return CALCPRIME_OUTPUT_BINARY;
    case calcprime::PrimeOutputFormat::ZstdDelta:
        return CALCPRIME_OUTPUT_ZSTD_DELTA;
    }
    return CALCPRIME_OUTPUT_TEXT;
}

calcprime_segment_config to_c_segment_config(const calcprime::SegmentConfig&config) {
    calcprime_segment_config out{};
    out.segment_bytes=config.segment_bytes;
    out.tile_bytes=config.tile_bytes;
    out.segment_bits=config.segment_bits;
    out.tile_bits=config.tile_bits;
    out.segment_span=config.segment_span;
    out.tile_span=config.tile_span;
    return out;
}

calcprime_cpu_info to_c_cpu_info(const calcprime::CpuInfo&info) {
    calcprime_cpu_info out{};
    out.logical_cpus=info.logical_cpus;
    out.physical_cpus=info.physical_cpus;
    out.l1_data_bytes=info.l1_data_bytes;
    out.l2_bytes=info.l2_bytes;
    out.l2_total_bytes=info.l2_total_bytes;
    out.has_smt=info.has_smt ? 1 : 0;
    return out;
}

struct SegmentResult {
    std::uint64_t count=0;
    std::vector<std::uint64_t>primes;
    std::atomic<bool>ready{false};
};

struct RangeOptions {
    std::uint64_t from=0;
    std::uint64_t to=0;
    unsigned threads=0;
    calcprime::WheelType wheel=calcprime::WheelType::Mod30;
    std::size_t segment_bytes=0;
    std::size_t tile_bytes=0;
    std::uint64_t nth_index=0;
    bool collect_primes=false;
    bool use_meissel=false;
    bool write_to_file=false;
    calcprime::PrimeOutputFormat output_format=calcprime::PrimeOutputFormat::Text;
    std::string output_path;
    calcprime_prime_chunk_callback prime_callback=nullptr;
    void*prime_user_data=nullptr;
    calcprime_progress_callback progress_callback=nullptr;
    void*progress_user_data=nullptr;
    calcprime_cancel_token*cancel_token=nullptr;
};

RangeOptions make_range_options(const calcprime_range_options&opts) {
    RangeOptions result;
    result.from=opts.from;
    result.to=opts.to;
    result.threads=opts.threads;
    result.wheel=to_cpp_wheel(opts.wheel);
    result.segment_bytes=opts.segment_bytes;
    result.tile_bytes=opts.tile_bytes;
    result.nth_index=opts.nth_index;
    result.collect_primes=opts.collect_primes!=0;
    result.use_meissel=opts.use_meissel!=0;
    result.write_to_file=opts.write_to_file!=0;
    result.output_format=to_cpp_output(opts.output_format);
    if(opts.output_path) {
        result.output_path=opts.output_path;
    }
    result.prime_callback=opts.prime_callback;
    result.prime_user_data=opts.prime_callback_user_data;
    result.progress_callback=opts.progress_callback;
    result.progress_user_data=opts.progress_user_data;
    result.cancel_token=opts.cancel_token;
    return result;
}

}

struct calcprime_cancel_token {
    std::atomic<bool>cancelled{false};
};

struct calcprime_range_run_result {
    calcprime_status status=CALCPRIME_STATUS_SUCCESS;
    calcprime_range_stats stats{};
    std::uint64_t total_count=0;
    int nth_found=0;
    std::uint64_t nth_value=0;
    bool primes_collected=false;
    std::vector<std::vector<std::uint64_t>>prime_chunks;
    std::uint64_t stored_prime_total=0;
    std::string error_message;
};

namespace calcprime {
int run_cli(int argc,char**argv);
}

extern"C" int calcprime_run_cli(int argc,char**argv) {
    return calcprime::run_cli(argc,argv);
}

extern"C" std::uint64_t calcprime_meissel_count(std::uint64_t from,std::uint64_t to,unsigned threads) {
    std::uint64_t sqrt_limit=0;
    if(to>1) {
        sqrt_limit=static_cast<std::uint64_t>(std::sqrt(static_cast<long double>(to)))+
                     1;
    }
    auto primes=calcprime::simple_sieve(sqrt_limit);
    return calcprime::meissel_count(from,to,primes,threads);
}

extern"C" int calcprime_miller_rabin_is_prime(std::uint64_t n) {
    return calcprime::miller_rabin_is_prime(n) ? 1 : 0;
}

extern"C" int calcprime_simple_sieve(std::uint64_t limit,std::uint32_t**out_primes,std::size_t*out_count) {
    if(!out_primes||!out_count) {
        return-1;
    }
    *out_primes=nullptr;
    *out_count=0;
    auto primes=calcprime::simple_sieve(limit);
    if(primes.empty()) {
        return 0;
    }
    std::unique_ptr<std::uint32_t[]>buffer;
    try {
        buffer=std::make_unique<std::uint32_t[]>(primes.size());
    } catch(const std::bad_alloc&) {
        return-1;
    }
    std::copy(primes.begin(),primes.end(),buffer.get());
    *out_count=primes.size();
    *out_primes=buffer.release();
    return 0;
}

extern"C" void calcprime_release_u32_buffer(std::uint32_t*buffer) {
    delete[] buffer;
}

extern"C" std::uint64_t calcprime_meissel_count_with_primes(std::uint64_t from,std::uint64_t to,const std::uint32_t*primes,std::size_t prime_count,unsigned threads) {
    std::vector<std::uint32_t>prime_vec;
    if(primes&&prime_count>0) {
        prime_vec.assign(primes,primes+prime_count);
    }
    return calcprime::meissel_count(from,to,prime_vec,threads);
}

extern"C" calcprime_cpu_info calcprime_detect_cpu_info(void) {
    auto cpp_info=calcprime::detect_cpu_info();
    calcprime_cpu_info info{};
    info.logical_cpus=cpp_info.logical_cpus;
    info.physical_cpus=cpp_info.physical_cpus;
    info.l1_data_bytes=cpp_info.l1_data_bytes;
    info.l2_bytes=cpp_info.l2_bytes;
    info.l2_total_bytes=cpp_info.l2_total_bytes;
    info.has_smt=cpp_info.has_smt ? 1 : 0;
    return info;
}

extern"C" unsigned calcprime_effective_thread_count(const calcprime_cpu_info*info) {
    if(!info) {
        return 0;
    }
    calcprime::CpuInfo cpp_info{};
    cpp_info.logical_cpus=info->logical_cpus;
    cpp_info.physical_cpus=info->physical_cpus;
    cpp_info.l1_data_bytes=info->l1_data_bytes;
    cpp_info.l2_bytes=info->l2_bytes;
    cpp_info.l2_total_bytes=info->l2_total_bytes;
    cpp_info.has_smt=info->has_smt!=0;
    return calcprime::effective_thread_count(cpp_info);
}

extern"C" calcprime_segment_config calcprime_choose_segment_config(const calcprime_cpu_info*info,unsigned threads,std::size_t requested_segment_bytes,std::size_t requested_tile_bytes,std::uint64_t range_length) {
    calcprime_segment_config config{};
    calcprime::CpuInfo cpp_info{};
    if(info) {
        cpp_info.logical_cpus=info->logical_cpus;
        cpp_info.physical_cpus=info->physical_cpus;
        cpp_info.l1_data_bytes=info->l1_data_bytes;
        cpp_info.l2_bytes=info->l2_bytes;
        cpp_info.l2_total_bytes=info->l2_total_bytes;
        cpp_info.has_smt=info->has_smt!=0;
    }
    auto cpp_config=calcprime::choose_segment_config(cpp_info,threads,requested_segment_bytes,requested_tile_bytes,range_length);
    config.segment_bytes=cpp_config.segment_bytes;
    config.tile_bytes=cpp_config.tile_bytes;
    config.segment_bits=cpp_config.segment_bits;
    config.tile_bits=cpp_config.tile_bits;
    config.segment_span=cpp_config.segment_span;
    config.tile_span=cpp_config.tile_span;
    return config;
}

extern"C" int calcprime_range_options_init(calcprime_range_options*options) {
    if(!options) {
        return-1;
    }
    options->from=0;
    options->to=0;
    options->threads=0;
    options->wheel=CALCPRIME_WHEEL_MOD30;
    options->segment_bytes=0;
    options->tile_bytes=0;
    options->nth_index=0;
    options->collect_primes=0;
    options->use_meissel=0;
    options->write_to_file=0;
    options->output_format=CALCPRIME_OUTPUT_TEXT;
    options->output_path=nullptr;
    options->prime_callback=nullptr;
    options->prime_callback_user_data=nullptr;
    options->progress_callback=nullptr;
    options->progress_user_data=nullptr;
    options->cancel_token=nullptr;
    return 0;
}

extern"C" calcprime_cancel_token*calcprime_cancel_token_create(void) {
    return new (std::nothrow) calcprime_cancel_token();
}

extern"C" void calcprime_cancel_token_destroy(calcprime_cancel_token*token) {
    delete token;
}

extern"C" void calcprime_cancel_token_request(calcprime_cancel_token*token) {
    if(!token) {
        return;
    }
    token->cancelled.store(true,std::memory_order_release);
}

namespace {

enum class FailureKind {
    None,
    Writer,
    PrimeCallback,
    Progress
};

calcprime_status assign_exception(calcprime_range_run_result&result,FailureKind kind,const std::exception_ptr&ex_ptr) {
    if(ex_ptr) {
        try {
            std::rethrow_exception(ex_ptr);
        } catch(const std::exception&ex) {
            result.error_message=ex.what();
        } catch(...) {
            result.error_message="unknown error";
        }
    }
    switch(kind) {
    case FailureKind::Writer:
        return CALCPRIME_STATUS_IO_ERROR;
    case FailureKind::PrimeCallback:
    case FailureKind::Progress:
        return CALCPRIME_STATUS_INTERNAL_ERROR;
    case FailureKind::None:
        break;
    }
    return CALCPRIME_STATUS_INTERNAL_ERROR;
}

}

extern"C" calcprime_status calcprime_run_range(const calcprime_range_options*options,calcprime_range_run_result**out_result) {
    if(!out_result) {
        return CALCPRIME_STATUS_INVALID_ARGUMENT;
    }

    auto result=std::make_unique<calcprime_range_run_result>();
    result->status=CALCPRIME_STATUS_INVALID_ARGUMENT;

    if(!options) {
        result->error_message="options pointer is null";
        *out_result=result.release();
        return CALCPRIME_STATUS_INVALID_ARGUMENT;
    }

    if(!is_valid_wheel(options->wheel)) {
        result->error_message="invalid wheel selection";
        *out_result=result.release();
        return CALCPRIME_STATUS_INVALID_ARGUMENT;
    }
    if(!is_valid_output_format(options->output_format)) {
        result->error_message="invalid output format";
        *out_result=result.release();
        return CALCPRIME_STATUS_INVALID_ARGUMENT;
    }

    RangeOptions opts=make_range_options(*options);

    result->status=CALCPRIME_STATUS_SUCCESS;
    result->stats.from=opts.from;
    result->stats.to=opts.to;
    result->stats.elapsed_us=0;
    result->stats.threads=0;
    result->stats.cpu=calcprime_cpu_info{};
    result->stats.segment=calcprime_segment_config{};
    result->stats.wheel=options->wheel;
    result->stats.output_format=options->output_format;
    result->stats.segments_total=0;
    result->stats.segments_processed=0;
    result->stats.prime_count=0;
    result->stats.nth_index=opts.nth_index;
    result->stats.nth_found=0;
    result->stats.use_meissel=opts.use_meissel ? 1 : 0;
    result->stats.completed=0;
    result->stats.cancelled=0;
    result->primes_collected=opts.collect_primes;
    result->prime_chunks.clear();
    result->stored_prime_total=0;
    result->total_count=0;
    result->nth_found=0;
    result->nth_value=0;
    result->error_message.clear();

    if(opts.to<=opts.from||opts.to<2) {
        result->status=CALCPRIME_STATUS_INVALID_ARGUMENT;
        result->error_message="invalid range";
        *out_result=result.release();
        return CALCPRIME_STATUS_INVALID_ARGUMENT;
    }

    bool need_prime_delivery=opts.collect_primes||opts.write_to_file||(opts.prime_callback!=nullptr);
    if(opts.use_meissel&&(need_prime_delivery||opts.nth_index!=0)) {
        result->status=CALCPRIME_STATUS_INVALID_ARGUMENT;
        result->error_message="Meissel counting cannot emit primes";
        *out_result=result.release();
        return CALCPRIME_STATUS_INVALID_ARGUMENT;
    }

    calcprime::CpuInfo cpu_info=calcprime::detect_cpu_info();
    result->stats.cpu=to_c_cpu_info(cpu_info);

    unsigned threads=opts.threads ? opts.threads : calcprime::effective_thread_count(cpu_info);
    if(opts.nth_index!=0) {
        threads=1;
    }
    if(threads==0) {
        threads=1;
    }
    result->stats.threads=threads;

    auto start_time=std::chrono::steady_clock::now();

    if(opts.use_meissel) {
        try {
            std::uint64_t sqrt_limit=0;
            if(opts.to>1) {
                sqrt_limit=static_cast<std::uint64_t>(std::sqrt(static_cast<long double>(opts.to)))+
                             1;
            }
            auto primes=calcprime::simple_sieve(sqrt_limit);
            std::uint64_t count=calcprime::meissel_count(opts.from,opts.to,primes,threads);
            auto end_time=std::chrono::steady_clock::now();
            auto elapsed=std::chrono::duration_cast<std::chrono::microseconds>(end_time-start_time);
            result->total_count=count;
            result->stats.prime_count=count;
            result->stats.elapsed_us=static_cast<std::uint64_t>(elapsed.count());
            result->stats.segment=calcprime_segment_config{};
            result->stats.segments_total=0;
            result->stats.segments_processed=0;
            result->stats.completed=1;
        } catch(const std::exception&ex) {
            result->status=CALCPRIME_STATUS_INTERNAL_ERROR;
            result->error_message=ex.what();
        } catch(...) {
            result->status=CALCPRIME_STATUS_INTERNAL_ERROR;
            result->error_message="unknown error";
        }

        if(opts.progress_callback&&result->status==CALCPRIME_STATUS_SUCCESS) {
            try {
                opts.progress_callback(1.0,opts.progress_user_data);
            } catch(const std::exception&ex) {
                result->status=CALCPRIME_STATUS_INTERNAL_ERROR;
                result->error_message=ex.what();
            } catch(...) {
                result->status=CALCPRIME_STATUS_INTERNAL_ERROR;
                result->error_message="unknown error";
            }
        }

        *out_result=result.release();
        return (*out_result)->status;
    }

    calcprime::WheelType wheel_type=opts.wheel;
    const calcprime::Wheel&wheel=calcprime::get_wheel(wheel_type);

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

    calcprime::SieveRange range{odd_begin,odd_end};
    std::uint64_t length=(range.end>range.begin)?(range.end-range.begin):0;

    calcprime::SegmentConfig config=calcprime::choose_segment_config(cpu_info,threads,opts.segment_bytes,opts.tile_bytes,length);
    result->stats.segment=to_c_segment_config(config);

    std::size_t num_segments=length?static_cast<std::size_t>((length+config.segment_span-1)/config.segment_span):0;
    result->stats.segments_total=num_segments;

    std::uint64_t sqrt_limit=static_cast<std::uint64_t>(std::sqrt(static_cast<long double>(opts.to)))+
                            1;
    auto base_primes=calcprime::simple_sieve(sqrt_limit);

    std::uint32_t small_limit=29u;
    switch(wheel_type) {
    case calcprime::WheelType::Mod30:
        small_limit=29u;
        break;
    case calcprime::WheelType::Mod210:
        small_limit=47u;
        break;
    case calcprime::WheelType::Mod1155:
        small_limit=47u;
        break;
    }

    bool need_segment_storage=need_prime_delivery;
    bool need_primes_for_nth=opts.nth_index!=0;

    calcprime::PrimeMarker marker(wheel,config,range.begin,range.end,base_primes,small_limit);
    calcprime::SegmentWorkQueue queue(range,config);

    std::vector<SegmentResult>segment_results(num_segments);
    std::mutex segment_ready_mutex;
    std::condition_variable segment_ready_cv;
    std::atomic<bool>stop{false};
    std::atomic<bool>nth_found_flag{false};
    std::uint64_t nth_value=0;
    std::atomic<std::size_t>segments_processed{0};
    std::mutex progress_mutex;
    bool progress_cancelled=false;
    bool callback_cancelled=false;
    bool external_cancelled=false;
    FailureKind failure_kind=FailureKind::None;
    std::exception_ptr stored_exception;

    std::unique_ptr<calcprime::PrimeWriter>writer;
    if(opts.write_to_file) {
        try {
            writer=std::make_unique<calcprime::PrimeWriter>(true,opts.output_path,opts.output_format);
        } catch(const std::exception&ex) {
            result->status=CALCPRIME_STATUS_IO_ERROR;
            result->error_message=ex.what();
            *out_result=result.release();
            return (*out_result)->status;
        } catch(...) {
            result->status=CALCPRIME_STATUS_IO_ERROR;
            result->error_message="failed to initialize writer";
            *out_result=result.release();
            return (*out_result)->status;
        }
    }
    calcprime::PrimeWriter*writer_ptr=writer.get();

    auto deliver_chunk=[&](std::vector<std::uint64_t>&&chunk)->bool {
        if(chunk.empty()) {
            return true;
        }
        std::size_t chunk_size=chunk.size();
        if(writer_ptr) {
            try {
                writer_ptr->write_segment(chunk);
            } catch(...) {
                if(failure_kind==FailureKind::None) {
                    failure_kind=FailureKind::Writer;
                    stored_exception=std::current_exception();
                }
                return false;
            }
        }
        if(opts.prime_callback) {
            int callback_result=0;
            try {
                callback_result=opts.prime_callback(chunk.data(),chunk.size(),opts.prime_user_data);
            } catch(...) {
                if(failure_kind==FailureKind::None) {
                    failure_kind=FailureKind::PrimeCallback;
                    stored_exception=std::current_exception();
                }
                return false;
            }
            if(callback_result!=0) {
                callback_cancelled=true;
                return false;
            }
        }
        if(opts.collect_primes) {
            result->prime_chunks.emplace_back(std::move(chunk));
            result->stored_prime_total+=static_cast<std::uint64_t>(chunk_size);
        }
        return true;
    };

    std::vector<std::uint64_t>prefix_primes;
    if(opts.from<=2&&opts.to>2) {
        prefix_primes.push_back(2);
    }
    std::vector<std::uint64_t>wheel_primes;
    switch(wheel_type) {
    case calcprime::WheelType::Mod30:
        wheel_primes={3,5};
        break;
    case calcprime::WheelType::Mod210:
        wheel_primes={3,5,7};
        break;
    case calcprime::WheelType::Mod1155:
        wheel_primes={3,5,7,11};
        break;
    }
    for(std::uint64_t p : wheel_primes) {
        if(p>=opts.from&&p<opts.to) {
            prefix_primes.push_back(p);
        }
    }
    std::uint64_t prefix_count=static_cast<std::uint64_t>(prefix_primes.size());

    if(opts.nth_index!=0&&opts.nth_index<=prefix_count) {
        nth_value=prefix_primes[static_cast<std::size_t>(opts.nth_index-1)];
        nth_found_flag.store(true,std::memory_order_release);
        stop.store(true,std::memory_order_release);
    }

    if(!prefix_primes.empty()) {
        if(!deliver_chunk(std::move(prefix_primes))) {
            stop.store(true,std::memory_order_release);
        }
    }

    if(opts.progress_callback) {
        try {
            opts.progress_callback(0.0,opts.progress_user_data);
        } catch(...) {
            if(failure_kind==FailureKind::None) {
                failure_kind=FailureKind::Progress;
                stored_exception=std::current_exception();
            }
            stop.store(true,std::memory_order_release);
        }
    }

    std::vector<std::thread>workers;
    workers.reserve(threads);

    const std::uint64_t nth_target=opts.nth_index;
    const std::uint64_t prefix_total=prefix_count;

    for(unsigned t=0;t<threads;++t) {
        workers.emplace_back([&,t]() {
            auto state=marker.make_thread_state(t,threads);
            std::vector<std::uint64_t>bitset;
            std::uint64_t cumulative=prefix_total;
            while(!stop.load(std::memory_order_acquire)) {
                if(opts.cancel_token&&opts.cancel_token->cancelled.load(std::memory_order_acquire)) {
                    external_cancelled=true;
                    stop.store(true,std::memory_order_release);
                    break;
                }
                std::uint64_t segment_id=0;
                std::uint64_t seg_low=0;
                std::uint64_t seg_high=0;
                if(!queue.next(segment_id,seg_low,seg_high)) {
                    break;
                }
                marker.sieve_segment(state,segment_id,seg_low,seg_high,bitset);
                std::size_t bit_count=static_cast<std::size_t>((seg_high-seg_low)>>1);
                std::uint64_t local_count=calcprime::count_zero_bits(bitset.data(),bit_count);
                if(segment_id<segment_results.size()) {
                    segment_results[segment_id].count=local_count;
                }

                std::vector<std::uint64_t>primes;
                bool need_primes=need_segment_storage||(need_primes_for_nth&&threads==1);
                if(need_primes&&local_count>0) {
                    primes.reserve(static_cast<std::size_t>(local_count));
                    std::uint64_t value=seg_low;
                    std::size_t produced=0;
                    for(std::size_t word=0;word<bitset.size()&&produced<bit_count;++word) {
                        std::uint64_t composite=bitset[word];
                        for(std::size_t bit=0;bit<64&&produced<bit_count;
                             ++bit,++produced,value+=2) {
                            if(composite&(1ULL<<bit)) {
                                continue;
                            }
                            primes.push_back(value);
                        }
                    }
                }

                if(need_primes_for_nth&&threads==1&&!nth_found_flag.load(std::memory_order_acquire)) {
                    std::uint64_t base=cumulative;
                    std::uint64_t new_total=base+local_count;
                    if(nth_target>base&&nth_target<=new_total) {
                        if(primes.empty()&&local_count>0) {
                            primes.reserve(static_cast<std::size_t>(local_count));
                            std::uint64_t value=seg_low;
                            std::size_t produced=0;
                            for(std::size_t word=0;word<bitset.size()&&produced<bit_count;++word) {
                                std::uint64_t composite=bitset[word];
                                for(std::size_t bit=0;bit<64&&produced<bit_count;
                                     ++bit,++produced,value+=2) {
                                    if(composite&(1ULL<<bit)) {
                                        continue;
                                    }
                                    primes.push_back(value);
                                }
                            }
                        }
                        std::size_t index=static_cast<std::size_t>(nth_target-base-1);
                        if(index<primes.size()) {
                            nth_value=primes[index];
                            nth_found_flag.store(true,std::memory_order_release);
                            stop.store(true,std::memory_order_release);
                        }
                    }
                    cumulative=new_total;
                }

                if(need_segment_storage&&segment_id<segment_results.size()) {
                    segment_results[segment_id].primes=std::move(primes);
                    segment_results[segment_id].ready.store(true,std::memory_order_release);
                    segment_ready_cv.notify_all();
                }

                std::size_t completed=segments_processed.fetch_add(1,std::memory_order_acq_rel)+1;
                if(opts.progress_callback&&!progress_cancelled) {
                    std::lock_guard<std::mutex>lock(progress_mutex);
                    if(!progress_cancelled) {
                        double progress_value=(num_segments==0)
                                                    ? 1.0
                                                    : static_cast<double>(completed)/
                                                          static_cast<double>(num_segments);
                        if(progress_value>1.0) {
                            progress_value=1.0;
                        }
                        int progress_result=0;
                        try {
                            progress_result=opts.progress_callback(progress_value,opts.progress_user_data);
                        } catch(...) {
                            if(failure_kind==FailureKind::None) {
                                failure_kind=FailureKind::Progress;
                                stored_exception=std::current_exception();
                            }
                            progress_cancelled=true;
                            stop.store(true,std::memory_order_release);
                            break;
                        }
                        if(progress_result!=0) {
                            progress_cancelled=true;
                            stop.store(true,std::memory_order_release);
                        }
                    }
                }
            }
        });
    }

    std::thread delivery_thread;
    if(need_segment_storage&&num_segments>0) {
        delivery_thread=std::thread([&]() {
            for(std::size_t idx=0;idx<num_segments;++idx) {
                SegmentResult&seg=segment_results[idx];
                std::vector<std::uint64_t>primes;
                {
                    std::unique_lock<std::mutex>lock(segment_ready_mutex);
                    segment_ready_cv.wait(lock,[&]() {
                        return seg.ready.load(std::memory_order_acquire)||
                               stop.load(std::memory_order_acquire);
                    });
                    bool ready=seg.ready.load(std::memory_order_acquire);
                    if(!ready) {
                        if(stop.load(std::memory_order_acquire)) {
                            break;
                        }
                        continue;
                    }
                    seg.ready.store(false,std::memory_order_release);
                    primes=std::move(seg.primes);
                }
                if(!deliver_chunk(std::move(primes))) {
                    stop.store(true,std::memory_order_release);
                    break;
                }
            }
            if(writer_ptr) {
                try {
                    writer_ptr->flush();
                } catch(...) {
                    if(failure_kind==FailureKind::None) {
                        failure_kind=FailureKind::Writer;
                        stored_exception=std::current_exception();
                    }
                    stop.store(true,std::memory_order_release);
                }
            }
        });
    }

    for(auto&worker : workers) {
        if(worker.joinable()) {
            worker.join();
        }
    }

    segment_ready_cv.notify_all();

    if(delivery_thread.joinable()) {
        delivery_thread.join();
    }

    if(writer_ptr) {
        try {
            writer_ptr->finish();
        } catch(...) {
            if(failure_kind==FailureKind::None) {
                failure_kind=FailureKind::Writer;
                stored_exception=std::current_exception();
            }
        }
    }

    std::size_t processed=segments_processed.load(std::memory_order_acquire);
    result->stats.segments_processed=processed;

    std::uint64_t total=prefix_count;
    for(const auto&seg : segment_results) {
        total+=seg.count;
    }
    result->total_count=total;
    result->stats.prime_count=total;

    bool nth_found=nth_found_flag.load(std::memory_order_acquire);
    if(nth_found) {
        result->nth_found=1;
        result->stats.nth_found=1;
        result->nth_value=nth_value;
    }

    bool cancelled=external_cancelled||progress_cancelled||callback_cancelled;
    result->stats.cancelled=cancelled?1:0;

    if(failure_kind!=FailureKind::None) {
        result->status=assign_exception(*result,failure_kind,stored_exception);
    } else if(cancelled) {
        result->status=CALCPRIME_STATUS_CANCELLED;
        if(callback_cancelled) {
            result->error_message="prime callback requested cancellation";
        } else if(progress_cancelled) {
            result->error_message="progress callback requested cancellation";
        } else {
            result->error_message="operation cancelled";
        }
    } else {
        result->status=CALCPRIME_STATUS_SUCCESS;
    }

    if(opts.nth_index!=0&&!nth_found&&result->status==CALCPRIME_STATUS_SUCCESS) {
        result->status=CALCPRIME_STATUS_INTERNAL_ERROR;
        result->error_message="nth prime not found within range";
    }

    if(opts.progress_callback&&!progress_cancelled&&failure_kind==FailureKind::None&&!external_cancelled) {
        try {
            opts.progress_callback(1.0,opts.progress_user_data);
        } catch(...) {
            if(result->status==CALCPRIME_STATUS_SUCCESS) {
                try {
                    throw;
                } catch(const std::exception&ex) {
                    result->status=CALCPRIME_STATUS_INTERNAL_ERROR;
                    result->error_message=ex.what();
                } catch(...) {
                    result->status=CALCPRIME_STATUS_INTERNAL_ERROR;
                    result->error_message="unknown error";
                }
            }
        }
    }

    auto end_time=std::chrono::steady_clock::now();
    auto elapsed=std::chrono::duration_cast<std::chrono::microseconds>(end_time-start_time);
    result->stats.elapsed_us=static_cast<std::uint64_t>(elapsed.count());

    bool completed=(processed==num_segments)&&!cancelled&&failure_kind==FailureKind::None&&(!nth_found||num_segments==0);
    result->stats.completed=completed?1:0;

    *out_result=result.release();
    return (*out_result)->status;
}

extern"C" calcprime_status calcprime_range_result_status(const calcprime_range_run_result*result) {
    if(!result) {
        return CALCPRIME_STATUS_INVALID_ARGUMENT;
    }
    return result->status;
}

extern"C" const char* calcprime_range_result_error_message(const calcprime_range_run_result*result) {
    if(!result||result->error_message.empty()) {
        return nullptr;
    }
    return result->error_message.c_str();
}

extern"C" std::uint64_t calcprime_range_result_count(const calcprime_range_run_result*result) {
    if(!result) {
        return 0;
    }
    return result->total_count;
}

extern"C" int calcprime_range_result_nth_prime(const calcprime_range_run_result*result,std::uint64_t*out_value) {
    if(!result||!out_value) {
        return-1;
    }
    if(!result->nth_found) {
        *out_value=0;
        return-1;
    }
    *out_value=result->nth_value;
    return 0;
}

extern"C" int calcprime_range_result_stats(const calcprime_range_run_result*result,calcprime_range_stats*out_stats) {
    if(!result||!out_stats) {
        return-1;
    }
    *out_stats=result->stats;
    return 0;
}

extern"C" std::size_t calcprime_range_result_segment_count(const calcprime_range_run_result*result) {
    if(!result||!result->primes_collected) {
        return 0;
    }
    return result->prime_chunks.size();
}

extern"C" int calcprime_range_result_segment(const calcprime_range_run_result*result,std::size_t index,const std::uint64_t**out_primes,std::size_t*out_count) {
    if(!result||!result->primes_collected||!out_primes||!out_count) {
        return-1;
    }
    if(index>=result->prime_chunks.size()) {
        return-1;
    }
    const auto&chunk=result->prime_chunks[index];
    *out_primes=chunk.data();
    *out_count=chunk.size();
    return 0;
}

extern"C" int calcprime_range_result_copy_primes(const calcprime_range_run_result*result,std::uint64_t*buffer,std::size_t capacity,std::size_t*out_written) {
    if(!result||!out_written) {
        return-1;
    }
    *out_written=static_cast<std::size_t>(result->stored_prime_total);
    if(!result->primes_collected) {
        return-1;
    }
    if(!buffer) {
        return 0;
    }
    if(capacity<result->stored_prime_total) {
        return-1;
    }
    std::size_t offset=0;
    for(const auto&chunk : result->prime_chunks) {
        if(chunk.empty()) {
            continue;
        }
        std::copy(chunk.begin(),chunk.end(),buffer+offset);
        offset+=chunk.size();
    }
    return 0;
}

extern"C" void calcprime_range_result_release(calcprime_range_run_result*result) {
    delete result;
}

extern"C" std::uint64_t calcprime_popcount_u64(std::uint64_t value) {
    return calcprime::popcount_u64(value);
}

extern"C" std::uint64_t calcprime_count_zero_bits(const std::uint64_t*bits,std::size_t bit_count) {
    if(!bits) {
        return 0;
    }
    return calcprime::count_zero_bits(bits,bit_count);
}
