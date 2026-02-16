#include "segmenter.h"

#include<algorithm>
#include<cmath>
#include<limits>

namespace calcprime{
namespace{

std::size_t align_to(std::size_t value,std::size_t alignment){
	if(alignment==0){
		return value;
	}
	std::size_t remainder=value%alignment;
	if(remainder==0){
		return value;
	}
	std::size_t add=alignment-remainder;
	if(value>std::numeric_limits<std::size_t>::max()-add){
		std::size_t max_aligned=
			std::numeric_limits<std::size_t>::max()-
			(std::numeric_limits<std::size_t>::max()%alignment);
		return max_aligned;
	}
	return value+add;
}

std::size_t align_down(std::size_t value,std::size_t alignment){
	if(alignment==0||value==0){
		return value;
	}
	return value-(value%alignment);
}

std::size_t clamp_floor_to_size_t(long double value){
	if(!std::isfinite(value)||value<=0.0L){
		return 0;
	}
	long double max_size=
		static_cast<long double>(std::numeric_limits<std::size_t>::max());
	if(value>=max_size){
		return std::numeric_limits<std::size_t>::max();
	}
	return static_cast<std::size_t>(std::floor(value));
}

} // namespace

SegmentConfig choose_segment_config(const CpuInfo&info,unsigned threads,
									std::size_t requested_segment_bytes,
									std::size_t requested_tile_bytes,
									std::uint64_t range_length){
	std::size_t l1=info.l1_data_bytes?info.l1_data_bytes:32*1024;
	std::size_t l2=info.l2_bytes?info.l2_bytes:1024*1024;

	std::size_t thread_count=threads?static_cast<std::size_t>(threads):1;

	std::size_t total_l2=info.l2_total_bytes;
	if(!total_l2){
		std::size_t cores=
			info.physical_cpus?info.physical_cpus:info.logical_cpus;
		if(cores==0){
			cores=thread_count?thread_count:1;
		}
		if(l2>0){
			if(l2>std::numeric_limits<std::size_t>::max()/cores){
				total_l2=std::numeric_limits<std::size_t>::max();
			}else{
				total_l2=l2*cores;
			}
		}
	}

	std::size_t segment_bytes=requested_segment_bytes;
	std::size_t cap_limit_bytes=0;
	if(!segment_bytes){
		constexpr long double k0=1562.5L;
		constexpr long double beta=0.0625L;
		constexpr long double alpha_g=0.833333L;
		constexpr long double min_segment=8.0L*1024.0L;

		long double R=static_cast<long double>(range_length);
		long double s_fixed=0.0L;
		if(R>0.0L){
			long double scaled_R=R/1.0e10L;
			long double k_r=k0;
			if(scaled_R>0.0L){
				k_r*=std::pow(scaled_R,beta);
			}
			if(k_r>0.0L){
				s_fixed=R/(16.0L*k_r);
			}
		}

		long double s_min=0.0L;
		if(R>0.0L){
			if(R<=1.0e9L){
				long double ratio=R/1.0e8L;
				s_min=8.0L*1024.0L*std::pow(ratio,1.05L);
			}else{
				long double ratio=R/1.0e9L;
				s_min=90.0L*1024.0L*std::pow(ratio,-0.5L);
			}
		}

		long double base=std::max({min_segment,s_fixed,s_min});
		if(total_l2){
			long double s_max=static_cast<long double>(total_l2)*alpha_g;
			base=std::min(base,s_max);
			cap_limit_bytes=clamp_floor_to_size_t(s_max);
		}

		if(!std::isfinite(base)||base<=0.0L){
			base=min_segment;
		}

		long double max_size=
			static_cast<long double>(std::numeric_limits<std::size_t>::max());
		if(base>=max_size){
			segment_bytes=std::numeric_limits<std::size_t>::max();
		}else{
			long double rounded=std::floor(base+0.5L);
			if(rounded<=0.0L){
				rounded=min_segment;
			}
			if(rounded>=max_size){
				segment_bytes=std::numeric_limits<std::size_t>::max();
			}else{
				segment_bytes=align_to(static_cast<std::size_t>(rounded),128);
			}
		}

		if(segment_bytes==0){
			segment_bytes=8*1024;
		}
	}else{
		segment_bytes=align_to(requested_segment_bytes,128);
	}

	segment_bytes=align_to(segment_bytes,128);
	if(!requested_segment_bytes){
		if(thread_count<=1){
			segment_bytes=std::max<std::size_t>(segment_bytes,1024*1024);
		}else{
			segment_bytes=std::max<std::size_t>(segment_bytes,768*1024);
		}
	}
	if(cap_limit_bytes){
		std::size_t cap_aligned=align_down(cap_limit_bytes,128);
		if(cap_aligned==0&&cap_limit_bytes){
			cap_aligned=cap_limit_bytes;
		}
		if(cap_aligned&&segment_bytes>cap_aligned){
			segment_bytes=cap_aligned;
		}
	}
	if(segment_bytes<8*1024){
		segment_bytes=8*1024;
	}

	std::size_t tile_bytes=requested_tile_bytes;
	if(!tile_bytes){
		std::size_t target=std::max<std::size_t>(l1,8*1024);
		if(thread_count<=1&&target<64*1024){
			target=64*1024;
		}
		tile_bytes=align_to(target,128);
	}else{
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

SegmentConfig choose_worker_segment_config(const CpuInfo&info,
										   const SegmentConfig&base_config,
										   unsigned worker_index,
										   unsigned thread_count,
										   std::size_t requested_tile_bytes,
										   std::uint64_t range_span,
										   CoreSchedulingMode mode){
	(void)range_span;
	SegmentConfig config=base_config;
	if(requested_tile_bytes!=0||thread_count<=1||
	   mode==CoreSchedulingMode::Legacy){
		return config;
	}
	if(config.segment_bytes==0){
		return config;
	}
	if(!info.has_hybrid){
		return config;
	}
	bool performance_worker=
		is_performance_worker(info,worker_index,thread_count,mode);

	std::size_t class_l1=
		performance_worker?info.performance_l1_data_bytes
						  :info.efficiency_l1_data_bytes;
	if(!class_l1){
		class_l1=info.l1_data_bytes;
	}
	if(!class_l1){
		class_l1=32*1024;
	}
	std::size_t class_l2=
		performance_worker?info.performance_l2_bytes:info.efficiency_l2_bytes;
	if(!class_l2){
		class_l2=info.l2_bytes;
	}

	std::size_t tile_bytes=align_to(std::max<std::size_t>(class_l1,8*1024),128);
	if(class_l2){
		std::size_t l2_cap=align_down(
			std::max<std::size_t>(class_l2/4,8*1024),128);
		if(l2_cap){
			tile_bytes=std::min(tile_bytes,l2_cap);
		}
	}
	tile_bytes=std::min(tile_bytes,config.segment_bytes);
	if(tile_bytes<8*1024){
		tile_bytes=std::min(config.segment_bytes,static_cast<std::size_t>(8*1024));
	}
	tile_bytes=align_to(tile_bytes,128);
	if(tile_bytes==0){
		tile_bytes=std::min(config.segment_bytes,static_cast<std::size_t>(8*1024));
	}
	if(tile_bytes>config.segment_bytes){
		tile_bytes=config.segment_bytes;
	}
	if(tile_bytes==config.tile_bytes){
		return config;
	}
	config.tile_bytes=tile_bytes;
	config.tile_bits=tile_bytes*8;
	config.tile_span=static_cast<std::uint64_t>(config.tile_bits)*2ULL;
	return config;
}

SegmentWorkQueue::SegmentWorkQueue(SieveRange range,const SegmentConfig&config)
	: range_(range),config_(config),next_segment_(0){
	length_=(range_.end>range_.begin)?(range_.end-range_.begin):0;
	if(length_==0||config_.segment_span==0){
		total_segments_=0;
	}else{
		total_segments_=
			length_/config_.segment_span+
			((length_%config_.segment_span)!=0ULL?1ULL:0ULL);
	}
}

bool SegmentWorkQueue::next(std::uint64_t&segment_id,std::uint64_t&segment_low,
							std::uint64_t&segment_high){
	std::uint64_t begin_id=0;
	std::uint64_t end_id=0;
	if(!next_chunk(1,begin_id,end_id)||begin_id>=end_id){
		return false;
	}
	segment_id=begin_id;
	return segment_bounds(segment_id,segment_low,segment_high);
}

bool SegmentWorkQueue::next_chunk(std::uint64_t requested_segments,
								  std::uint64_t&segment_begin_id,
								  std::uint64_t&segment_end_id){
	segment_begin_id=0;
	segment_end_id=0;
	if(total_segments_==0){
		return false;
	}
	if(requested_segments==0){
		requested_segments=1;
	}
	std::uint64_t begin=
		next_segment_.fetch_add(requested_segments,std::memory_order_relaxed);
	if(begin>=total_segments_){
		return false;
	}
	std::uint64_t end=begin+requested_segments;
	if(end<begin||end>total_segments_){
		end=total_segments_;
	}
	segment_begin_id=begin;
	segment_end_id=end;
	return true;
}

bool SegmentWorkQueue::segment_bounds(std::uint64_t segment_id,
									  std::uint64_t&segment_low,
									  std::uint64_t&segment_high) const{
	if(config_.segment_span==0||segment_id>=total_segments_){
		return false;
	}
	if(segment_id>std::numeric_limits<std::uint64_t>::max()/
					  config_.segment_span){
		return false;
	}
	std::uint64_t offset=segment_id*config_.segment_span;
	if(offset>=length_){
		return false;
	}
	if(range_.begin>std::numeric_limits<std::uint64_t>::max()-offset){
		return false;
	}
	segment_low=range_.begin+offset;
	std::uint64_t remaining=length_-offset;
	std::uint64_t span_length=config_.segment_span;
	if(span_length>remaining){
		span_length=remaining;
	}
	if(segment_low>
	   std::numeric_limits<std::uint64_t>::max()-span_length){
		return false;
	}
	segment_high=segment_low+span_length;
	return segment_low<segment_high;
}

} // namespace calcprime
