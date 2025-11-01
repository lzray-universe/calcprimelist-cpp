#pragma once

#include "bucket.h"
#include "segmenter.h"
#include "wheel.h"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace calcprime {

struct LargePrimeState {
    std::uint32_t prime;
    std::uint64_t next_value;
    std::uint64_t stride;
};

struct TileView {
    std::uint64_t start_value;
    std::size_t bit_offset;
    std::size_t bit_count;
    std::uint64_t*word_ptr;
    std::size_t word_count;
};

class PrimeMarker {
public:
    PrimeMarker(const Wheel&wheel,SegmentConfig config,std::uint64_t range_begin,std::uint64_t range_end,const std::vector<std::uint32_t>&primes,std::uint32_t small_prime_limit=29);

    struct ThreadState {
        BucketRing bucket;
        std::vector<std::uint64_t>small_positions;
        std::vector<LargePrimeState>large_states;
        std::vector<std::uint64_t>medium_positions;
    };

    ThreadState make_thread_state(std::size_t thread_index,std::size_t thread_count) const;

    void sieve_segment(ThreadState&state,std::uint64_t segment_id,std::uint64_t segment_low,std::uint64_t segment_high,std::vector<std::uint64_t>&bitset) const;

    const SegmentConfig&config() const { return config_;}

private:
    const Wheel&wheel_;
    SegmentConfig config_;
    std::uint64_t range_begin_;
    std::uint64_t range_end_;
    std::vector<std::uint32_t>small_primes_;
    std::vector<std::uint64_t>small_initial_;
    std::vector<const SmallPrimePattern*>small_prime_patterns_;
    std::vector<std::uint32_t>medium_primes_;
    std::vector<std::uint64_t>medium_initial_;
    std::vector<LargePrimeState>large_primes_template_;

    static std::uint64_t first_hit(std::uint32_t prime,std::uint64_t start);
    void apply_small_primes(ThreadState&state,const TileView&tile) const;
    void apply_medium_primes(ThreadState&state,const TileView&tile,std::size_t segment_index) const;
    void apply_large_primes(ThreadState&state,std::uint64_t segment_id,std::uint64_t segment_low,std::uint64_t segment_high,std::vector<std::uint64_t>&bitset) const;
};

}
