#include "segmenter.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace calcprime {
namespace {

std::size_t align_to(std::size_t value,std::size_t alignment) {
    if(alignment==0) {
        return value;
    }
    std::size_t remainder=value%alignment;
    if(remainder==0) {
        return value;
    }
    std::size_t add=alignment-remainder;
    if(value>std::numeric_limits<std::size_t>::max()-add) {
        std::size_t max_aligned=std::numeric_limits<std::size_t>::max()-
                                  (std::numeric_limits<std::size_t>::max()%alignment);
        return max_aligned;
    }
    return value+add;
}

std::size_t align_down(std::size_t value,std::size_t alignment) {
    if(alignment==0||value==0) {
        return value;
    }
    return value-(value%alignment);
}

std::size_t clamp_floor_to_size_t(long double value) {
    if(!std::isfinite(value)||value<=0.0L) {
        return 0;
    }
    long double max_size=static_cast<long double>(std::numeric_limits<std::size_t>::max());
    if(value>=max_size) {
        return std::numeric_limits<std::size_t>::max();
    }
    return static_cast<std::size_t>(std::floor(value));
}

}

SegmentConfig choose_segment_config(const CpuInfo&info,unsigned threads,std::size_t requested_segment_bytes,std::size_t requested_tile_bytes,std::uint64_t range_length) {
    std::size_t l1=info.l1_data_bytes ? info.l1_data_bytes : 32*1024;
    std::size_t l2=info.l2_bytes ? info.l2_bytes : 1024*1024;

    std::size_t thread_count=threads ? static_cast<std::size_t>(threads) : 1;

    std::size_t total_l2=info.l2_total_bytes;
    if(!total_l2) {
        std::size_t cores=info.physical_cpus ? info.physical_cpus : info.logical_cpus;
        if(cores==0) {
            cores=thread_count ? thread_count : 1;
        }
        if(l2>0) {
            if(l2>std::numeric_limits<std::size_t>::max()/cores) {
                total_l2=std::numeric_limits<std::size_t>::max();
            } else {
                total_l2=l2*cores;
            }
        }
    }

    std::size_t segment_bytes=requested_segment_bytes;
    std::size_t cap_limit_bytes=0;
    if(!segment_bytes) {
        constexpr long double k0=1562.5L;
        constexpr long double beta=0.0625L;
        constexpr long double alpha_g=0.833333L;
        constexpr long double min_segment=8.0L*1024.0L;

        long double R=static_cast<long double>(range_length);
        long double s_fixed=0.0L;
        if(R>0.0L) {
            long double scaled_R=R/1.0e10L;
            long double k_r=k0;
            if(scaled_R>0.0L) {
                k_r*=std::pow(scaled_R,beta);
            }
            if(k_r>0.0L) {
                s_fixed=R/(16.0L*k_r);
            }
        }

        long double s_min=0.0L;
        if(R>0.0L) {
            if(R<=1.0e9L) {
                long double ratio=R/1.0e8L;
                s_min=8.0L*1024.0L*std::pow(ratio,1.05L);
            } else {
                long double ratio=R/1.0e9L;
                s_min=90.0L*1024.0L*std::pow(ratio,-0.5L);
            }
        }

        long double base=std::max({min_segment,s_fixed,s_min});
        if(total_l2) {
            long double s_max=static_cast<long double>(total_l2)*alpha_g;
            base=std::min(base,s_max);
            cap_limit_bytes=clamp_floor_to_size_t(s_max);
        }

        if(!std::isfinite(base)||base<=0.0L) {
            base=min_segment;
        }

        long double max_size=static_cast<long double>(std::numeric_limits<std::size_t>::max());
        if(base>=max_size) {
            segment_bytes=std::numeric_limits<std::size_t>::max();
        } else {
            long double rounded=std::floor(base+0.5L);
            if(rounded<=0.0L) {
                rounded=min_segment;
            }
            if(rounded>=max_size) {
                segment_bytes=std::numeric_limits<std::size_t>::max();
            } else {
                segment_bytes=align_to(static_cast<std::size_t>(rounded),128);
            }
        }

        if(segment_bytes==0) {
            segment_bytes=8*1024;
        }
    } else {
        segment_bytes=align_to(requested_segment_bytes,128);
    }

    segment_bytes=align_to(segment_bytes,128);
    if(cap_limit_bytes) {
        std::size_t cap_aligned=align_down(cap_limit_bytes,128);
        if(cap_aligned==0&&cap_limit_bytes) {
            cap_aligned=cap_limit_bytes;
        }
        if(cap_aligned&&segment_bytes>cap_aligned) {
            segment_bytes=cap_aligned;
        }
    }
    if(segment_bytes<8*1024) {
        segment_bytes=8*1024;
    }

    std::size_t tile_bytes=requested_tile_bytes;
    if(!tile_bytes) {
        std::size_t target=std::max<std::size_t>(l1,8*1024);
        tile_bytes=align_to(target,128);
    } else {
        tile_bytes=align_to(requested_tile_bytes,128);
    }
    tile_bytes=std::min(tile_bytes,segment_bytes);

    SegmentConfig config{};
    config.segment_bytes=segment_bytes;
    config.tile_bytes=tile_bytes;
    config.segment_bits=segment_bytes*8;
    config.tile_bits=tile_bytes*8;
    config.segment_span=static_cast<std::uint64_t>(config.segment_bits)*2ULL;
    config.tile_span=static_cast<std::uint64_t>(config.tile_bits)*2ULL;
    return config;
}

SegmentWorkQueue::SegmentWorkQueue(SieveRange range,const SegmentConfig&config)
    : range_(range),config_(config),next_segment_(0) {
    length_=(range_.end>range_.begin)?(range_.end-range_.begin):0;
}

bool SegmentWorkQueue::next(std::uint64_t&segment_id,std::uint64_t&segment_low,std::uint64_t&segment_high) {
    std::uint64_t idx=next_segment_.fetch_add(1,std::memory_order_relaxed);
    std::uint64_t span=config_.segment_span;
    std::uint64_t offset=idx*span;
    if(span!=0&&offset/span!=idx) {
        return false;
    }
    if(offset>=length_) {
        return false;
    }
    segment_id=idx;
    segment_low=range_.begin+offset;
    std::uint64_t remaining=length_-offset;
    std::uint64_t span_length=span;
    if(span_length>remaining) {
        span_length=remaining;
    }
    segment_high=range_.begin+offset+span_length;
    return segment_low<segment_high;
}

}
