#pragma once

#include "cpu_info.h"

#include <atomic>
#include <cstddef>
#include <cstdint>

namespace calcprime {

struct SieveRange {
    std::uint64_t begin;
    std::uint64_t end;
};

struct SegmentConfig {
    std::size_t segment_bytes;
    std::size_t tile_bytes;
    std::size_t segment_bits;
    std::size_t tile_bits;
    std::uint64_t segment_span;
    std::uint64_t tile_span;
};

SegmentConfig choose_segment_config(const CpuInfo&info,unsigned threads,std::size_t requested_segment_bytes,std::size_t requested_tile_bytes,std::uint64_t range_length);

class SegmentWorkQueue {
public:
    SegmentWorkQueue(SieveRange range,const SegmentConfig&config);

    bool next(std::uint64_t&segment_id,std::uint64_t&segment_low,std::uint64_t&segment_high);

private:
    SieveRange range_;
    SegmentConfig config_;
    std::atomic<std::uint64_t>next_segment_;
    std::uint64_t length_;
};

}
