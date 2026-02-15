#include "marker.h"

#include<algorithm>

namespace calcprime{
namespace{

std::size_t words_for_bits(std::size_t bits){ return (bits+63)/64; }
std::uint64_t ceil_div_u64(std::uint64_t value,std::uint64_t divisor){
	return value/divisor+((value%divisor)!=0ULL?1ULL:0ULL);
}

const SmallPrimePattern*find_small_pattern(const Wheel&wheel,
										   std::uint32_t prime){
	for(const auto&pattern : wheel.small_patterns){
		if(pattern.prime==prime){
			return &pattern;
		}
	}
	return nullptr;
}

} // namespace

std::uint64_t PrimeMarker::first_hit(std::uint32_t prime,std::uint64_t start){
	std::uint64_t begin=
		static_cast<std::uint64_t>(prime)*static_cast<std::uint64_t>(prime);
	if(begin<start){
		begin=start;
	}
	std::uint64_t remainder=begin%prime;
	if(remainder){
		begin+=prime-remainder;
	}
	if((begin&1ULL)==0){
		begin+=prime;
	}
	return begin;
}

PrimeMarker::PrimeMarker(const Wheel&wheel,SegmentConfig config,
						 std::uint64_t range_begin,std::uint64_t range_end,
						 const std::vector<std::uint32_t>&primes,
						 std::uint32_t small_prime_limit)
	: wheel_(wheel),config_(config),range_begin_(range_begin),
	  range_end_(range_end){
	std::uint64_t large_threshold=config_.segment_span/2ULL;
	for(std::uint32_t prime : primes){
		if(prime<2){
			continue;
		}
		if(prime==2){
			continue; // even numbers already excluded
		}
		if(wheel_.presieve_modulus%prime==0){
			continue; // already removed by wheel presieve
		}
		if(prime<=small_prime_limit){
			small_primes_.push_back(prime);
			small_initial_.push_back(first_hit(prime,range_begin_));
			small_prime_patterns_.push_back(find_small_pattern(wheel_,prime));
		}else if(static_cast<std::uint64_t>(prime)<=large_threshold){
			medium_primes_.push_back(prime);
			medium_initial_.push_back(first_hit(prime,range_begin_));
		}else{
			LargePrimeState state;
			state.prime=prime;
			state.stride=static_cast<std::uint64_t>(prime)*2ULL;
			state.next_value=first_hit(prime,range_begin_);
			large_primes_template_.push_back(state);
		}
	}
}

PrimeMarker::ThreadState
PrimeMarker::make_thread_state(std::size_t thread_index,
							   std::size_t thread_count) const{
	if(thread_count==0){
		thread_count=1;
	}
	ThreadState state;
	state.bucket.reset(0);
	state.small_positions=small_initial_;
	state.medium_positions=medium_initial_;
	state.medium_next.assign(medium_primes_.size(),-1);
	std::size_t count=0;
	for(std::size_t i=0;i<large_primes_template_.size();++i){
		if(i%thread_count==thread_index){
			++count;
		}
	}
	state.large_states.reserve(count);
	for(std::size_t i=0;i<large_primes_template_.size();++i){
		if(i%thread_count!=thread_index){
			continue;
		}
		state.large_states.push_back(large_primes_template_[i]);
		auto&lp=state.large_states.back();
		if(lp.next_value>=range_end_){
			continue;
		}
		std::uint64_t segment=(lp.next_value-range_begin_)/config_.segment_span;
		std::uint64_t base=range_begin_+segment*config_.segment_span;
		if((base&1ULL)==0){
			++base;
		}
		std::uint64_t offset=(lp.next_value-base)>>1;
		state.bucket.push(
			segment,BucketEntry{lp.prime,segment,offset,lp.next_value,&lp});
	}
	return state;
}

void PrimeMarker::apply_small_primes(ThreadState&state,
									 const TileView&tile) const{
	if(tile.bit_count==0){
		return;
	}
	const std::uint32_t tile_bits=static_cast<std::uint32_t>(tile.bit_count);
	std::uint64_t tile_end=tile.start_value+tile.bit_count*2ULL;
	for(std::size_t i=0;i<small_primes_.size();++i){
		std::uint32_t prime=small_primes_[i];
		std::uint64_t step=static_cast<std::uint64_t>(prime)*2ULL;
		std::uint64_t pos=state.small_positions[i];
		if(pos<tile.start_value){
			std::uint64_t delta=tile.start_value-pos;
			std::uint64_t skip=(delta+step-1)/step;
			pos+=skip*step;
		}
		if(pos>=tile_end){
			state.small_positions[i]=pos;
			continue;
		}
		const SmallPrimePattern*pattern=small_prime_patterns_[i];
		if(pattern){
			std::size_t bit_index=
				static_cast<std::size_t>((pos-tile.start_value)>>1);
			std::size_t word_index=bit_index/64;
			if(word_index<tile.word_count){
				std::size_t bit_in_word=bit_index%64;
				std::uint32_t phase=pattern->start_phase[bit_in_word];
				std::uint64_t*word_ptr=tile.word_ptr+word_index;
				std::uint64_t*word_end=tile.word_ptr+tile.word_count;
				std::uint64_t mask=pattern->masks[phase];
				if(bit_in_word!=0){
					mask&=(~0ULL)<<bit_in_word;
				}
				word_ptr[0]|=mask;
				phase=pattern->next_phase[phase];
				for(++word_ptr;word_ptr<word_end;++word_ptr){
					word_ptr[0]|=pattern->masks[phase];
					phase=pattern->next_phase[phase];
				}
			}
			std::uint64_t delta=tile_end-pos;
			std::uint64_t skip=(delta+step-1)/step;
			pos+=skip*step;
			state.small_positions[i]=pos;
		}else{
			std::uint32_t bit_index=
				static_cast<std::uint32_t>((pos-tile.start_value)>>1);
			std::uint32_t bit_step=prime;
			while(bit_index<tile_bits){
				tile.word_ptr[bit_index>>6]|=(1ULL<<(bit_index&63ULL));
				bit_index+=bit_step;
			}
			state.small_positions[i]=
				tile.start_value+(static_cast<std::uint64_t>(bit_index)<<1);
		}
	}
}

void PrimeMarker::apply_medium_primes(ThreadState&state,
									  const TileView&tile,
									  std::uint64_t segment_low,
									  std::uint64_t segment_high,
									  std::size_t tile_index,
									  std::size_t tile_count) const{
	if(tile.bit_count==0){
		return;
	}
	if(tile_index>=state.medium_tile_heads.size()){
		return;
	}
	std::int32_t head=state.medium_tile_heads[tile_index];
	state.medium_tile_heads[tile_index]=-1;
	if(head<0){
		return;
	}
	const std::uint32_t tile_bits=static_cast<std::uint32_t>(tile.bit_count);
	std::uint64_t tile_end=tile.start_value+tile.bit_count*2ULL;
	while(head>=0){
		std::size_t i=static_cast<std::size_t>(head);
		head=state.medium_next[i];
		std::uint32_t prime=medium_primes_[i];
		std::uint64_t step=static_cast<std::uint64_t>(prime)*2ULL;
		std::uint64_t pos=state.medium_positions[i];
		if(pos<tile.start_value){
			std::uint64_t delta=tile.start_value-pos;
			std::uint64_t skip=(delta+step-1)/step;
			pos+=skip*step;
		}
		if(pos<tile_end){
			std::uint32_t bit_index=
				static_cast<std::uint32_t>((pos-tile.start_value)>>1);
			std::uint32_t bit_step=prime;
			while(bit_index<tile_bits){
				tile.word_ptr[bit_index>>6]|=(1ULL<<(bit_index&63ULL));
				bit_index+=bit_step;
			}
			pos=tile.start_value+(static_cast<std::uint64_t>(bit_index)<<1);
		}
		state.medium_positions[i]=pos;
		if(pos<segment_high){
			std::size_t next_tile=static_cast<std::size_t>(
				(pos-segment_low)/config_.tile_span);
			if(next_tile>=tile_count){
				next_tile=tile_count-1;
			}
			state.medium_next[i]=state.medium_tile_heads[next_tile];
			state.medium_tile_heads[next_tile]=static_cast<std::int32_t>(i);
		}else{
			state.medium_next[i]=-1;
		}
	}
}

void PrimeMarker::apply_large_primes(ThreadState&state,std::uint64_t segment_id,
									 std::uint64_t segment_low,
									 std::uint64_t segment_high,
									 std::vector<std::uint64_t>&bitset) const{
	auto hits=state.bucket.take(segment_id);
	for(auto&entry : hits){
		if(entry.value>=segment_low&&entry.value<segment_high){
			std::size_t bit_index=(entry.value-segment_low)>>1;
			bitset[bit_index/64]|=(1ULL<<(bit_index%64));
		}
		if(!entry.owner){
			continue;
		}
		std::uint64_t next=entry.value+entry.owner->stride;
		entry.owner->next_value=next;
		if(next>=range_end_){
			continue;
		}
		std::uint64_t seg=(next-range_begin_)/config_.segment_span;
		std::uint64_t base=range_begin_+seg*config_.segment_span;
		if((base&1ULL)==0){
			++base;
		}
		std::uint64_t offset=(next-base)>>1;
		state.bucket.push(
			seg,BucketEntry{entry.owner->prime,seg,offset,next,entry.owner});
	}
}

void PrimeMarker::sieve_segment(ThreadState&state,std::uint64_t segment_id,
								std::uint64_t segment_low,
								std::uint64_t segment_high,
								std::vector<std::uint64_t>&bitset) const{
	if(segment_high<=segment_low){
		bitset.clear();
		return;
	}
	std::size_t bit_count=
		static_cast<std::size_t>((segment_high-segment_low)>>1);
	if(bit_count==0){
		bitset.clear();
		return;
	}
	std::size_t word_count=words_for_bits(bit_count);
	bitset.resize(word_count);

	wheel_.fill_presieve(segment_low,bit_count,bitset.data());
	apply_large_primes(state,segment_id,segment_low,segment_high,bitset);

	std::size_t tile_count=
		static_cast<std::size_t>(ceil_div_u64(
			segment_high-segment_low,config_.tile_span));
	if(tile_count==0){
		tile_count=1;
	}
	state.medium_tile_heads.assign(tile_count,-1);
	if(state.medium_next.size()!=medium_primes_.size()){
		state.medium_next.assign(medium_primes_.size(),-1);
	}
	for(std::size_t i=0;i<medium_primes_.size();++i){
		std::uint32_t prime=medium_primes_[i];
		std::uint64_t step=static_cast<std::uint64_t>(prime)*2ULL;
		std::uint64_t pos=state.medium_positions[i];
		if(pos<segment_low){
			std::uint64_t delta=segment_low-pos;
			std::uint64_t skip=(delta+step-1)/step;
			pos+=skip*step;
			state.medium_positions[i]=pos;
		}
		if(pos>=segment_high){
			state.medium_next[i]=-1;
			continue;
		}
		std::size_t next_tile=static_cast<std::size_t>(
			(pos-segment_low)/config_.tile_span);
		if(next_tile>=tile_count){
			next_tile=tile_count-1;
		}
		state.medium_next[i]=state.medium_tile_heads[next_tile];
		state.medium_tile_heads[next_tile]=static_cast<std::int32_t>(i);
	}

	std::uint64_t tile_low=segment_low;
	std::size_t bit_offset=0;
	std::size_t tile_index=0;
	while(tile_low<segment_high){
		std::uint64_t tile_high=
			std::min<std::uint64_t>(segment_high,tile_low+config_.tile_span);
		std::size_t tile_bits=static_cast<std::size_t>((tile_high-tile_low)>>1);
		std::size_t tile_words=words_for_bits(tile_bits);
		TileView tile{tile_low,bit_offset,tile_bits,
					  bitset.data()+(bit_offset/64),tile_words};
		apply_small_primes(state,tile);
		apply_medium_primes(state,tile,segment_low,segment_high,tile_index,
							tile_count);
		if(tile_bits%64!=0&&tile_words>0){
			std::uint64_t mask=(1ULL<<(tile_bits%64))-1;
			tile.word_ptr[tile_words-1]&=mask;
		}
		tile_low=tile_high;
		bit_offset+=tile_bits;
		++tile_index;
	}
}

} // namespace calcprime
