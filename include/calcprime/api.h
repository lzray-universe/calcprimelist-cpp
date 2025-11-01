#pragma once

#include <cstddef>
#include <cstdint>

#if defined(_WIN32) && defined(CALCPRIME_DLL_EXPORT)
#  define CALCPRIME_API __declspec(dllexport)
#else
#  define CALCPRIME_API
#endif

#ifdef __cplusplus
extern"C" {
#endif

typedef struct calcprime_cpu_info {
    unsigned logical_cpus;
    unsigned physical_cpus;
    std::size_t l1_data_bytes;
    std::size_t l2_bytes;
    std::size_t l2_total_bytes;
    int has_smt;
} calcprime_cpu_info;

typedef struct calcprime_segment_config {
    std::size_t segment_bytes;
    std::size_t tile_bytes;
    std::size_t segment_bits;
    std::size_t tile_bits;
    std::uint64_t segment_span;
    std::uint64_t tile_span;
} calcprime_segment_config;

typedef enum calcprime_status {
    CALCPRIME_STATUS_SUCCESS=0,
    CALCPRIME_STATUS_INVALID_ARGUMENT=1,
    CALCPRIME_STATUS_CANCELLED=2,
    CALCPRIME_STATUS_IO_ERROR=3,
    CALCPRIME_STATUS_INTERNAL_ERROR=4
} calcprime_status;

typedef enum calcprime_wheel_type {
    CALCPRIME_WHEEL_MOD30=0,
    CALCPRIME_WHEEL_MOD210=1,
    CALCPRIME_WHEEL_MOD1155=2
} calcprime_wheel_type;

typedef enum calcprime_output_format {
    CALCPRIME_OUTPUT_TEXT=0,
    CALCPRIME_OUTPUT_BINARY=1,
    CALCPRIME_OUTPUT_ZSTD_DELTA=2
} calcprime_output_format;

struct calcprime_cancel_token;
typedef struct calcprime_cancel_token calcprime_cancel_token;

typedef int (*calcprime_prime_chunk_callback)(const std::uint64_t*primes,std::size_t count,void*user_data);

typedef int (*calcprime_progress_callback)(double progress,void*user_data);

typedef struct calcprime_range_options {
    std::uint64_t from;
    std::uint64_t to;
    unsigned threads;
    calcprime_wheel_type wheel;
    std::size_t segment_bytes;
    std::size_t tile_bytes;
    std::uint64_t nth_index;
    int collect_primes;
    int use_meissel;
    int write_to_file;
    calcprime_output_format output_format;
    const char*output_path;
    calcprime_prime_chunk_callback prime_callback;
    void*prime_callback_user_data;
    calcprime_progress_callback progress_callback;
    void*progress_user_data;
    calcprime_cancel_token*cancel_token;
} calcprime_range_options;

typedef struct calcprime_range_stats {
    std::uint64_t from;
    std::uint64_t to;
    std::uint64_t elapsed_us;
    unsigned threads;
    calcprime_cpu_info cpu;
    calcprime_segment_config segment;
    calcprime_wheel_type wheel;
    calcprime_output_format output_format;
    std::size_t segments_total;
    std::size_t segments_processed;
    std::uint64_t prime_count;
    std::uint64_t nth_index;
    int nth_found;
    int use_meissel;
    int completed;
    int cancelled;
} calcprime_range_stats;

struct calcprime_range_run_result;
typedef struct calcprime_range_run_result calcprime_range_run_result;

CALCPRIME_API int calcprime_run_cli(int argc,char**argv);

CALCPRIME_API std::uint64_t calcprime_meissel_count(std::uint64_t from,std::uint64_t to,unsigned threads);

CALCPRIME_API int calcprime_miller_rabin_is_prime(std::uint64_t n);

CALCPRIME_API int calcprime_simple_sieve(std::uint64_t limit,std::uint32_t**out_primes,std::size_t*out_count);

CALCPRIME_API void calcprime_release_u32_buffer(std::uint32_t*buffer);

CALCPRIME_API std::uint64_t calcprime_meissel_count_with_primes(std::uint64_t from,std::uint64_t to,const std::uint32_t*primes,std::size_t prime_count,unsigned threads);

CALCPRIME_API calcprime_cpu_info calcprime_detect_cpu_info(void);

CALCPRIME_API unsigned calcprime_effective_thread_count(const calcprime_cpu_info*info);

CALCPRIME_API calcprime_segment_config calcprime_choose_segment_config(const calcprime_cpu_info*info,unsigned threads,std::size_t requested_segment_bytes,std::size_t requested_tile_bytes,std::uint64_t range_length);

CALCPRIME_API int calcprime_range_options_init(calcprime_range_options*options);

CALCPRIME_API calcprime_cancel_token*calcprime_cancel_token_create(void);
CALCPRIME_API void calcprime_cancel_token_destroy(calcprime_cancel_token*token);
CALCPRIME_API void calcprime_cancel_token_request(calcprime_cancel_token*token);

CALCPRIME_API calcprime_status calcprime_run_range(const calcprime_range_options*options,calcprime_range_run_result**out_result);

CALCPRIME_API calcprime_status calcprime_range_result_status(const calcprime_range_run_result*result);

CALCPRIME_API const char* calcprime_range_result_error_message(const calcprime_range_run_result*result);

CALCPRIME_API std::uint64_t calcprime_range_result_count(const calcprime_range_run_result*result);

CALCPRIME_API int calcprime_range_result_nth_prime(const calcprime_range_run_result*result,std::uint64_t*out_value);

CALCPRIME_API int calcprime_range_result_stats(const calcprime_range_run_result*result,calcprime_range_stats*out_stats);

CALCPRIME_API std::size_t calcprime_range_result_segment_count(const calcprime_range_run_result*result);

CALCPRIME_API int calcprime_range_result_segment(const calcprime_range_run_result*result,std::size_t index,const std::uint64_t**out_primes,std::size_t*out_count);

CALCPRIME_API int calcprime_range_result_copy_primes(const calcprime_range_run_result*result,std::uint64_t*buffer,std::size_t capacity,std::size_t*out_written);

CALCPRIME_API void calcprime_range_result_release(calcprime_range_run_result*result);

CALCPRIME_API std::uint64_t calcprime_popcount_u64(std::uint64_t value);

CALCPRIME_API std::uint64_t calcprime_count_zero_bits(const std::uint64_t*bits,std::size_t bit_count);

#ifdef __cplusplus
}
#endif

#undef CALCPRIME_API
