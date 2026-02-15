#include "wheel_bitmap_count.h"

#include "popcnt.h"

#include<algorithm>
#include<array>
#include<atomic>
#include<cstdint>
#include<limits>
#include<numeric>
#include<thread>
#include<vector>

namespace calcprime{
namespace{

struct WheelBitmapKernel{
	std::uint32_t modulus=0;
	std::uint32_t wheel_factor_max=0;
	std::uint32_t min_value=0;
	std::uint16_t phase_count=0;
	std::uint64_t full_mask_u64=0;
	std::vector<std::uint32_t> prefix_primes;
	std::vector<std::uint16_t> residues;
	std::vector<std::uint16_t> steps;
	std::vector<std::int16_t> residue_to_bit;
	std::vector<std::int16_t> residue_to_phase;
};

struct PrimeState{
	std::uint64_t block=0;
	std::uint8_t phase=0;
	std::array<std::uint32_t,48> block_inc{};
	std::array<std::uint8_t,48> bit_for_phase{};
};

struct PrimeState30{
	std::uint64_t block=0;
	std::uint8_t phase=0;
	std::array<std::uint32_t,8> block_inc{};
	std::array<std::uint8_t,8> bit_for_phase{};
};

struct DensePattern30{
	std::uint32_t prime=0;
	std::uint16_t phase_count=0;
	std::uint64_t start_block=0;
	std::vector<std::uint8_t> block_masks;
	std::vector<std::uint16_t> next_phase_block;
	std::vector<std::uint64_t> word_masks;
	std::vector<std::uint16_t> next_phase_word;
};

struct DensePattern210{
	std::uint32_t prime=0;
	std::uint64_t start_block=0;
	std::vector<std::uint64_t> block_masks;
	std::vector<std::uint16_t> next_phase_block;
};

constexpr std::array<std::uint8_t,256> build_popcnt8_table(){
	std::array<std::uint8_t,256> table{};
	for(std::size_t i=0;i<table.size();++i){
		std::uint8_t value=static_cast<std::uint8_t>(i);
		std::uint8_t count=0;
		while(value){
			value&=static_cast<std::uint8_t>(value-1);
			++count;
		}
		table[i]=count;
	}
	return table;
}

constexpr auto kPopcnt8Table=build_popcnt8_table();

constexpr std::array<std::array<std::uint64_t,8>,8> build_packed_bit_masks30(){
	std::array<std::array<std::uint64_t,8>,8> table{};
	for(std::size_t slot=0;slot<8;++slot){
		for(std::size_t bit=0;bit<8;++bit){
			table[slot][bit]=
				(1ULL<<static_cast<unsigned>(bit+slot*8));
		}
	}
	return table;
}

constexpr auto kPackedBitMask30=build_packed_bit_masks30();

std::uint64_t ceil_div_u64(std::uint64_t value,std::uint64_t divisor){
	return value/divisor+((value%divisor)!=0ULL?1ULL:0ULL);
}

WheelBitmapKernel build_kernel(std::uint32_t modulus,
							   std::vector<std::uint32_t> prefix_primes,
							   std::uint32_t min_value){
	WheelBitmapKernel kernel;
	kernel.modulus=modulus;
	kernel.prefix_primes=std::move(prefix_primes);
	kernel.min_value=min_value;
	kernel.wheel_factor_max=
		kernel.prefix_primes.empty()?0u:kernel.prefix_primes.back();
	kernel.residue_to_bit.assign(modulus,-1);
	kernel.residue_to_phase.assign(modulus,-1);

	for(std::uint32_t r=0;r<modulus;++r){
		if(std::gcd(r,modulus)==1){
			std::int16_t idx=
				static_cast<std::int16_t>(kernel.residues.size());
			kernel.residues.push_back(static_cast<std::uint16_t>(r));
			kernel.residue_to_bit[r]=idx;
			kernel.residue_to_phase[r]=idx;
		}
	}

	kernel.steps.reserve(kernel.residues.size());
	for(std::size_t i=0;i<kernel.residues.size();++i){
		std::uint32_t a=kernel.residues[i];
		std::uint32_t b=kernel.residues[(i+1)%kernel.residues.size()];
		std::uint32_t step=(b+modulus-a)%modulus;
		if(step==0){
			step=modulus;
		}
		kernel.steps.push_back(static_cast<std::uint16_t>(step));
	}

	kernel.phase_count=static_cast<std::uint16_t>(kernel.steps.size());
	kernel.full_mask_u64=0;
	for(std::size_t i=0;i<kernel.residues.size();++i){
		kernel.full_mask_u64|=(1ULL<<i);
	}
	return kernel;
}

const WheelBitmapKernel&kernel30(){
	static const WheelBitmapKernel kernel=
		build_kernel(30,{2,3,5},7);
	return kernel;
}

const WheelBitmapKernel&kernel210(){
	static const WheelBitmapKernel kernel=
		build_kernel(210,{2,3,5,7},11);
	return kernel;
}

template<typename MarkT> std::uint32_t popcount_mark(MarkT value);

template<> std::uint32_t popcount_mark<std::uint8_t>(std::uint8_t value){
	return kPopcnt8Table[value];
}

template<> std::uint32_t popcount_mark<std::uint64_t>(std::uint64_t value){
	return static_cast<std::uint32_t>(popcount_u64(value));
}

template<typename MarkT>
MarkT build_valid_mask(const WheelBitmapKernel&kernel,std::uint64_t block_base,
					   std::uint64_t low,std::uint64_t high){
	MarkT mask=0;
	for(std::size_t i=0;i<kernel.residues.size();++i){
		std::uint64_t value=
			block_base+static_cast<std::uint64_t>(kernel.residues[i]);
		if(value>=low&&value<high){
			mask|=static_cast<MarkT>(static_cast<MarkT>(1)<<i);
		}
	}
	return mask;
}

void advance_state_once(const WheelBitmapKernel&kernel,PrimeState&state){
	std::uint8_t idx=state.phase;
	state.block+=state.block_inc[idx];
	++idx;
	if(idx==kernel.phase_count){
		idx=0;
	}
	state.phase=idx;
}

inline void advance_state_once30(PrimeState30&state){
	std::uint8_t idx=state.phase;
	state.block+=state.block_inc[idx];
	state.phase=static_cast<std::uint8_t>((idx+1U)&7U);
}

inline void mark_sparse_state30(PrimeState30&state,
								std::uint64_t seg_block_begin,
								std::uint64_t seg_block_end,
								std::uint64_t* marks_words){
	std::uint64_t block=state.block;
	std::uint8_t phase=state.phase;
	const auto&inc=state.block_inc;
	const auto&bits=state.bit_for_phase;

	auto mark_phase=[&](std::uint8_t p){
		std::uint64_t rel=block-seg_block_begin;
		marks_words[static_cast<std::size_t>(rel>>3)]|=
			kPackedBitMask30[static_cast<std::size_t>(rel&7ULL)]
							[static_cast<std::size_t>(bits[p])];
		block+=inc[p];
	};

	while(phase!=0U&&block<seg_block_end){
		mark_phase(phase);
		phase=static_cast<std::uint8_t>((phase+1U)&7U);
	}

	while(block<seg_block_end){
		mark_phase(0); if(block>=seg_block_end){ phase=1; break; }
		mark_phase(1); if(block>=seg_block_end){ phase=2; break; }
		mark_phase(2); if(block>=seg_block_end){ phase=3; break; }
		mark_phase(3); if(block>=seg_block_end){ phase=4; break; }
		mark_phase(4); if(block>=seg_block_end){ phase=5; break; }
		mark_phase(5); if(block>=seg_block_end){ phase=6; break; }
		mark_phase(6); if(block>=seg_block_end){ phase=7; break; }
		mark_phase(7); if(block>=seg_block_end){ phase=0; break; }
		phase=0;
	}

	state.block=block;
	state.phase=phase;
}

bool initialize_state(const WheelBitmapKernel&kernel,std::uint32_t prime,
					  std::uint64_t chunk_low,std::uint64_t chunk_high,
					  PrimeState&state){
	if(chunk_high<=chunk_low){
		return false;
	}
	std::uint64_t p64=static_cast<std::uint64_t>(prime);
	if(p64==0||p64>(chunk_high-1ULL)/p64){
		return false; // no hit from p*p in this chunk
	}

	std::uint64_t value=p64*p64;
	if(value<chunk_low){
		std::uint64_t remainder=chunk_low%p64;
		value=remainder?chunk_low+(p64-remainder):chunk_low;
	}
	if((value&1ULL)==0){
		if(value>std::numeric_limits<std::uint64_t>::max()-p64){
			return false;
		}
		value+=p64;
	}
	if(value>=chunk_high){
		return false;
	}

	std::uint64_t odd_stride=p64*2ULL;
	while(value<chunk_high){
		std::uint64_t q=value/p64;
		std::int16_t phase=
			kernel.residue_to_phase[static_cast<std::size_t>(q%
														kernel.modulus)];
		if(phase>=0){
			std::uint64_t block=value/kernel.modulus;
			state.block=block;
			state.phase=static_cast<std::uint8_t>(phase);
			std::uint64_t pmod=static_cast<std::uint64_t>(prime)%kernel.modulus;
			for(std::size_t i=0;i<kernel.phase_count;++i){
				std::uint64_t nres=(pmod*kernel.residues[i])%kernel.modulus;
				std::int16_t bit=
					kernel.residue_to_bit[static_cast<std::size_t>(nres)];
				state.bit_for_phase[i]=static_cast<std::uint8_t>(bit);
				std::uint64_t delta=
					static_cast<std::uint64_t>(prime)*
					static_cast<std::uint64_t>(kernel.steps[i]);
				state.block_inc[i]=static_cast<std::uint32_t>(
					(nres+delta)/kernel.modulus);
			}
			return true;
		}
		if(value>std::numeric_limits<std::uint64_t>::max()-odd_stride){
			return false;
		}
		value+=odd_stride;
	}
	return false;
}

bool initialize_state30(const WheelBitmapKernel&kernel,std::uint32_t prime,
						std::uint64_t chunk_low,std::uint64_t chunk_high,
						PrimeState30&state){
	if(chunk_high<=chunk_low){
		return false;
	}
	std::uint64_t p64=static_cast<std::uint64_t>(prime);
	if(p64==0||p64>(chunk_high-1ULL)/p64){
		return false;
	}

	std::uint64_t value=p64*p64;
	if(value<chunk_low){
		std::uint64_t remainder=chunk_low%p64;
		value=remainder?chunk_low+(p64-remainder):chunk_low;
	}
	if((value&1ULL)==0){
		if(value>std::numeric_limits<std::uint64_t>::max()-p64){
			return false;
		}
		value+=p64;
	}
	if(value>=chunk_high){
		return false;
	}

	const std::uint64_t odd_stride=p64*2ULL;
	while(value<chunk_high){
		std::uint64_t q=value/p64;
		std::int16_t phase=
			kernel.residue_to_phase[static_cast<std::size_t>(q%30ULL)];
		if(phase>=0){
			std::uint64_t block=value/30ULL;
			state.block=block;
			state.phase=static_cast<std::uint8_t>(phase);
			std::uint64_t pmod=static_cast<std::uint64_t>(prime)%30ULL;
			for(std::size_t i=0;i<8;++i){
				std::uint64_t nres=(pmod*kernel.residues[i])%30ULL;
				std::int16_t bit=
					kernel.residue_to_bit[static_cast<std::size_t>(nres)];
				state.bit_for_phase[i]=static_cast<std::uint8_t>(bit);
				std::uint64_t delta=
					static_cast<std::uint64_t>(prime)*
					static_cast<std::uint64_t>(kernel.steps[i]);
				state.block_inc[i]=
					static_cast<std::uint32_t>((nres+delta)/30ULL);
			}
			return true;
		}
		if(value>std::numeric_limits<std::uint64_t>::max()-odd_stride){
			return false;
		}
		value+=odd_stride;
	}
	return false;
}

std::vector<DensePattern30>
build_dense_patterns30(const WheelBitmapKernel&kernel,
					   const std::vector<std::uint32_t>&base_primes,
					   std::uint32_t dense_limit){
	std::vector<DensePattern30> patterns;
	patterns.reserve(32);
	for(std::uint32_t prime : base_primes){
		if(prime<=kernel.wheel_factor_max){
			continue;
		}
		if(prime>dense_limit){
			break;
		}
		if(prime>std::numeric_limits<std::uint64_t>::max()/prime){
			break;
		}
		DensePattern30 pattern{};
		pattern.prime=prime;
		pattern.phase_count=static_cast<std::uint16_t>(prime);
		pattern.start_block=
			(static_cast<std::uint64_t>(prime)*
			 static_cast<std::uint64_t>(prime))/
			kernel.modulus;
		pattern.block_masks.assign(prime,0);
		pattern.next_phase_block.assign(prime,0);
		pattern.word_masks.assign(prime,0);
		pattern.next_phase_word.assign(prime,0);

		std::uint32_t mod_step=kernel.modulus%prime;
		for(std::uint32_t phase=0;phase<prime;++phase){
			std::uint32_t rem=phase;
			std::uint8_t block_mask=0;
			for(std::size_t i=0;i<kernel.residues.size();++i){
				std::uint32_t residue=kernel.residues[i]%prime;
				if((rem+residue)%prime==0){
					block_mask|=
						static_cast<std::uint8_t>(1u<<i);
				}
			}
			pattern.block_masks[phase]=block_mask;
			pattern.next_phase_block[phase]=
				static_cast<std::uint16_t>((phase+mod_step)%prime);
		}

		for(std::uint32_t phase=0;phase<prime;++phase){
			std::uint32_t rem=phase;
			std::uint64_t packed=0;
			for(std::size_t b=0;b<8;++b){
				packed|=(static_cast<std::uint64_t>(
							pattern.block_masks[rem])
						<<(b*8));
				rem=(rem+mod_step)%prime;
			}
			pattern.word_masks[phase]=packed;
			pattern.next_phase_word[phase]=
				static_cast<std::uint16_t>(rem);
		}

		patterns.push_back(std::move(pattern));
	}
	return patterns;
}

std::vector<DensePattern210>
build_dense_patterns210(const WheelBitmapKernel&kernel,
						const std::vector<std::uint32_t>&base_primes,
						std::uint32_t dense_limit){
	std::vector<DensePattern210> patterns;
	patterns.reserve(24);
	for(std::uint32_t prime : base_primes){
		if(prime<=kernel.wheel_factor_max){
			continue;
		}
		if(prime>dense_limit){
			break;
		}
		if(prime>std::numeric_limits<std::uint64_t>::max()/prime){
			break;
		}
		if(static_cast<std::uint64_t>(prime)*
			   static_cast<std::uint64_t>(prime)<
		   kernel.modulus){
			continue;
		}

		DensePattern210 pattern{};
		pattern.prime=prime;
		pattern.start_block=
			(static_cast<std::uint64_t>(prime)*
			 static_cast<std::uint64_t>(prime))/
			kernel.modulus;
		pattern.block_masks.assign(prime,0ULL);
		pattern.next_phase_block.assign(prime,0);

		const std::uint32_t mod_step=kernel.modulus%prime;
		for(std::uint32_t phase=0;phase<prime;++phase){
			std::uint64_t mask=0ULL;
			std::uint32_t rem=phase;
			for(std::size_t i=0;i<kernel.residues.size();++i){
				std::uint32_t residue=kernel.residues[i]%prime;
				if((rem+residue)%prime==0){
					mask|=(1ULL<<i);
				}
			}
			pattern.block_masks[phase]=mask;
			pattern.next_phase_block[phase]=
				static_cast<std::uint16_t>((phase+mod_step)%prime);
		}

		patterns.push_back(std::move(pattern));
	}
	return patterns;
}

std::uint64_t count_with_kernel30_fast(const WheelBitmapKernel&kernel,
									   std::uint64_t from,std::uint64_t to,
									   unsigned threads,
									   const SegmentConfig&config,
									   bool user_segment_override,
									   const std::vector<std::uint32_t>&
										   base_primes){
	std::uint64_t total=0;
	for(std::uint32_t prime : kernel.prefix_primes){
		std::uint64_t p64=static_cast<std::uint64_t>(prime);
		if(p64>=from&&p64<to){
			++total;
		}
	}

	std::uint64_t low=std::max<std::uint64_t>(from,kernel.min_value);
	if(to<=low){
		return total;
	}

	std::uint64_t block_begin=low/kernel.modulus;
	std::uint64_t block_end=ceil_div_u64(to,kernel.modulus);
	if(block_end<=block_begin){
		return total;
	}

	std::uint64_t segment_blocks=config.segment_bytes;
	if(segment_blocks==0){
		segment_blocks=1;
	}

	std::uint64_t total_blocks=block_end-block_begin;
	std::uint64_t total_segments=ceil_div_u64(total_blocks,segment_blocks);
	if(total_segments==0){
		return total;
	}

	unsigned worker_count=threads?threads:1;
	if(static_cast<std::uint64_t>(worker_count)>total_segments){
		worker_count=static_cast<unsigned>(total_segments);
	}
	if(worker_count==0){
		worker_count=1;
	}
	if(!user_segment_override&&worker_count>1){
		const std::uint64_t tuned_segment_blocks=512ULL*1024ULL;
		if(segment_blocks>tuned_segment_blocks){
			segment_blocks=tuned_segment_blocks;
			total_segments=ceil_div_u64(total_blocks,segment_blocks);
		}
	}

	constexpr std::uint32_t kDenseLimit=97;
	const auto dense_patterns=
		build_dense_patterns30(kernel,base_primes,kDenseLimit);
	const bool use_dynamic_schedule=(worker_count>1);
	std::uint64_t task_segments=0;
	std::uint64_t task_count=0;
	std::atomic<std::uint64_t> next_task{0};
	if(use_dynamic_schedule){
		task_segments=
			total_segments/(static_cast<std::uint64_t>(worker_count)*8ULL);
		if(task_segments<1ULL){
			task_segments=1ULL;
		}
		if(task_segments>32ULL){
			task_segments=32ULL;
		}
		task_count=ceil_div_u64(total_segments,task_segments);
	}

	std::atomic<std::uint64_t> total_acc{total};
	std::vector<std::thread> workers;
	workers.reserve(worker_count);

	for(unsigned worker=0;worker<worker_count;++worker){
		workers.emplace_back([&,worker](){
			std::vector<PrimeState30> sparse_states;
			sparse_states.reserve(base_primes.size());

			std::vector<std::uint64_t> marks_words(
				static_cast<std::size_t>((segment_blocks+7ULL)/8ULL),0);
			std::uint64_t local_total=0;
			const bool low_aligned=(low%kernel.modulus)==0ULL;
			const bool high_aligned=(to%kernel.modulus)==0ULL;
			struct DenseRun30{
				const DensePattern30* pattern=nullptr;
				std::size_t word_start=0;
				std::uint16_t phase=0;
			};
			struct DenseTail30{
				const DensePattern30* pattern=nullptr;
				std::size_t idx=0;
				std::uint16_t phase=0;
			};
			std::vector<DenseRun30> dense_runs;
			std::vector<DenseTail30> dense_tails;
			dense_runs.reserve(dense_patterns.size());
			dense_tails.reserve(dense_patterns.size());

			auto process_segment_range=
				[&](std::uint64_t seg_first,std::uint64_t seg_last){
					if(seg_first>=seg_last){
						return;
					}

					std::uint64_t chunk_block_begin=
						block_begin+seg_first*segment_blocks;
					std::uint64_t chunk_block_end=std::min<std::uint64_t>(
						block_end,block_begin+seg_last*segment_blocks);
					std::uint64_t chunk_low=
						chunk_block_begin*kernel.modulus;
					std::uint64_t chunk_high=
						(chunk_block_end==block_end)?to:
													 chunk_block_end*kernel.modulus;

					sparse_states.clear();
					for(std::uint32_t prime : base_primes){
						if(prime<=std::max(kernel.wheel_factor_max,kDenseLimit)){
							continue;
						}
						if(chunk_high>0&&
						   static_cast<std::uint64_t>(prime)>
							   (chunk_high-1ULL)/static_cast<std::uint64_t>(prime)){
							break;
						}
						PrimeState30 state{};
						if(initialize_state30(kernel,prime,chunk_low,chunk_high,
											  state)){
							sparse_states.push_back(state);
						}
					}

					for(std::uint64_t seg=seg_first;seg<seg_last;++seg){
						std::uint64_t seg_block_begin=
							block_begin+seg*segment_blocks;
						std::uint64_t seg_block_end=std::min<std::uint64_t>(
							block_end,seg_block_begin+segment_blocks);
						std::uint64_t block_count=seg_block_end-seg_block_begin;
						if(block_count==0){
							continue;
						}
						std::size_t word_count=
							static_cast<std::size_t>((block_count+7ULL)/8ULL);
						std::fill_n(marks_words.begin(),word_count,0ULL);

						dense_runs.clear();
						dense_tails.clear();
						std::size_t dense_full_words=
							static_cast<std::size_t>(block_count>>3);
						std::size_t dense_tail_start=dense_full_words<<3;

						for(const auto&pattern : dense_patterns){
							if(pattern.start_block>=seg_block_end){
								continue;
							}
							std::uint64_t begin_block=std::max<std::uint64_t>(
								seg_block_begin,pattern.start_block);
							if(begin_block>=seg_block_end){
								continue;
							}
							std::uint64_t rel=begin_block-seg_block_begin;
							std::size_t idx=static_cast<std::size_t>(rel);
							std::uint32_t phase=
								static_cast<std::uint32_t>(
									(begin_block*kernel.modulus)%pattern.prime);

							while(idx<block_count&&((idx&7u)!=0u)){
								marks_words[idx>>3]|=
									(static_cast<std::uint64_t>(
										 pattern.block_masks[phase])
									 <<((idx&7u)*8u));
								phase=pattern.next_phase_block[phase];
								++idx;
							}

							if(idx<dense_tail_start){
								dense_runs.push_back(DenseRun30{
									&pattern,idx>>3,
									static_cast<std::uint16_t>(phase)});
							}else if(idx<block_count){
								dense_tails.push_back(DenseTail30{
									&pattern,idx,
									static_cast<std::uint16_t>(phase)});
							}
						}

						for(std::size_t word_idx=0;word_idx<dense_full_words;
							++word_idx){
							std::uint64_t dense_word=0ULL;
							for(auto&run : dense_runs){
								if(word_idx<run.word_start){
									continue;
								}
								dense_word|=run.pattern->word_masks[run.phase];
								run.phase=run.pattern->next_phase_word[run.phase];
							}
							marks_words[word_idx]|=dense_word;
						}

						if(dense_tail_start<block_count){
							for(const auto&tail : dense_tails){
								std::size_t idx=tail.idx;
								std::uint32_t phase=tail.phase;
								while(idx<block_count){
									marks_words[idx>>3]|=
										(static_cast<std::uint64_t>(
											 tail.pattern->block_masks[phase])
										 <<((idx&7u)*8u));
									phase=tail.pattern->next_phase_block[phase];
									++idx;
								}
							}
							for(const auto&run : dense_runs){
								std::size_t idx=dense_tail_start;
								std::uint32_t phase=run.phase;
								while(idx<block_count){
									marks_words[idx>>3]|=
										(static_cast<std::uint64_t>(
											 run.pattern->block_masks[phase])
										 <<((idx&7u)*8u));
									phase=run.pattern->next_phase_block[phase];
									++idx;
								}
							}
						}

						for(auto&state : sparse_states){
							while(state.block<seg_block_begin){
								advance_state_once30(state);
							}
							if(state.block<seg_block_end){
								mark_sparse_state30(
									state,seg_block_begin,seg_block_end,
									marks_words.data());
							}
						}

						bool lower_full=(seg_block_begin>block_begin)||
										(low_aligned&&seg_block_begin==block_begin);
						bool upper_full=(seg_block_end<block_end)||
										(high_aligned&&seg_block_end==block_end);
						bool segment_full=lower_full&&upper_full;

						std::size_t full_words=
							static_cast<std::size_t>(block_count>>3);
						if(segment_full){
							for(std::size_t word_index=0;word_index<full_words;
								++word_index){
								local_total+=
									64ULL-popcount_u64(marks_words[word_index]);
							}

							std::size_t tail_index=full_words<<3;
							while(tail_index<block_count){
								std::uint64_t word=marks_words[tail_index>>3];
								std::uint8_t composite=static_cast<std::uint8_t>(
									(word>>((tail_index&7u)*8u))&0xFFu);
								local_total+=static_cast<std::uint64_t>(
									8U-kPopcnt8Table[composite]);
								++tail_index;
							}
						}else{
							std::size_t word_index=0;
							for(;word_index<full_words;++word_index){
								std::uint64_t abs_block=
									seg_block_begin+
									static_cast<std::uint64_t>(word_index)*8ULL;
								std::uint64_t word=marks_words[word_index];
								for(std::size_t b=0;b<8;++b){
									std::uint64_t block_idx=
										abs_block+static_cast<std::uint64_t>(b);
									std::uint64_t block_base=
										block_idx*kernel.modulus;
									std::uint8_t valid=
										build_valid_mask<std::uint8_t>(
											kernel,block_base,low,to);
									if(valid==0){
										continue;
									}
									std::uint8_t composite=
										static_cast<std::uint8_t>(
											(word>>(b*8u))&0xFFu);
									composite=
										static_cast<std::uint8_t>(composite&valid);
									local_total+=static_cast<std::uint64_t>(
										kPopcnt8Table[valid]-
										kPopcnt8Table[composite]);
								}
							}

							std::size_t tail_index=full_words<<3;
							while(tail_index<block_count){
								std::uint64_t abs_block=
									seg_block_begin+
									static_cast<std::uint64_t>(tail_index);
								std::uint64_t block_base=abs_block*kernel.modulus;
								std::uint8_t valid=build_valid_mask<std::uint8_t>(
									kernel,block_base,low,to);
								if(valid){
									std::uint64_t word=marks_words[tail_index>>3];
									std::uint8_t composite=
										static_cast<std::uint8_t>(
											(word>>((tail_index&7u)*8u))&0xFFu);
									composite=
										static_cast<std::uint8_t>(composite&valid);
									local_total+=static_cast<std::uint64_t>(
										kPopcnt8Table[valid]-
										kPopcnt8Table[composite]);
								}
								++tail_index;
							}
						}
					}
				};

			if(use_dynamic_schedule){
				while(true){
					std::uint64_t task=
						next_task.fetch_add(1ULL,std::memory_order_relaxed);
					if(task>=task_count){
						break;
					}
					std::uint64_t seg_first=task*task_segments;
					std::uint64_t seg_last=std::min<std::uint64_t>(
						total_segments,seg_first+task_segments);
					process_segment_range(seg_first,seg_last);
				}
			}else{
				std::uint64_t seg_first=
					(total_segments*worker)/worker_count;
				std::uint64_t seg_last=
					(total_segments*(worker+1ULL))/worker_count;
				process_segment_range(seg_first,seg_last);
			}

			total_acc.fetch_add(local_total,std::memory_order_relaxed);
		});
	}

	for(auto&worker : workers){
		worker.join();
	}
	return total_acc.load(std::memory_order_relaxed);
}

std::uint64_t count_with_kernel210_fast(const WheelBitmapKernel&kernel,
										std::uint64_t from,std::uint64_t to,
										unsigned threads,
										const SegmentConfig&config,
										bool user_segment_override,
										const std::vector<std::uint32_t>&
											base_primes){
	std::uint64_t total=0;
	for(std::uint32_t prime : kernel.prefix_primes){
		std::uint64_t p64=static_cast<std::uint64_t>(prime);
		if(p64>=from&&p64<to){
			++total;
		}
	}

	std::uint64_t low=std::max<std::uint64_t>(from,kernel.min_value);
	if(to<=low){
		return total;
	}

	std::uint64_t block_begin=low/kernel.modulus;
	std::uint64_t block_end=ceil_div_u64(to,kernel.modulus);
	if(block_end<=block_begin){
		return total;
	}

	std::uint64_t segment_blocks=config.segment_bytes/sizeof(std::uint64_t);
	if(segment_blocks==0){
		segment_blocks=1;
	}

	std::uint64_t total_blocks=block_end-block_begin;
	std::uint64_t total_segments=ceil_div_u64(total_blocks,segment_blocks);
	if(total_segments==0){
		return total;
	}

	unsigned worker_count=threads?threads:1;
	if(static_cast<std::uint64_t>(worker_count)>total_segments){
		worker_count=static_cast<unsigned>(total_segments);
	}
	if(worker_count==0){
		worker_count=1;
	}
	if(!user_segment_override&&worker_count>1){
		const std::uint64_t tuned_segment_blocks=
			(1ULL*1024ULL*1024ULL)/sizeof(std::uint64_t);
		if(tuned_segment_blocks>0&&segment_blocks>tuned_segment_blocks){
			segment_blocks=tuned_segment_blocks;
			total_segments=ceil_div_u64(total_blocks,segment_blocks);
		}
	}

	constexpr std::uint32_t kDenseLimit210=127;
	const auto dense_patterns=
		build_dense_patterns210(kernel,base_primes,kDenseLimit210);

	std::atomic<std::uint64_t> total_acc{total};
	std::vector<std::thread> workers;
	workers.reserve(worker_count);

	for(unsigned worker=0;worker<worker_count;++worker){
		workers.emplace_back([&,worker](){
			std::uint64_t seg_first=
				(total_segments*worker)/worker_count;
			std::uint64_t seg_last=
				(total_segments*(worker+1ULL))/worker_count;
			if(seg_first>=seg_last){
				return;
			}

			std::uint64_t chunk_block_begin=
				block_begin+seg_first*segment_blocks;
			std::uint64_t chunk_block_end=std::min<std::uint64_t>(
				block_end,block_begin+seg_last*segment_blocks);
			std::uint64_t chunk_low=
				chunk_block_begin*kernel.modulus;
			std::uint64_t chunk_high=
				(chunk_block_end==block_end)?to:
											 chunk_block_end*kernel.modulus;

			std::vector<PrimeState> states;
			states.reserve(base_primes.size());
			for(std::uint32_t prime : base_primes){
				if(prime<=kernel.wheel_factor_max){
					continue;
				}
				if(prime<=kDenseLimit210&&
				   static_cast<std::uint64_t>(prime)*
						   static_cast<std::uint64_t>(prime)>=
					   kernel.modulus){
					continue;
				}
				if(chunk_high>0&&
				   static_cast<std::uint64_t>(prime)>
					   (chunk_high-1ULL)/static_cast<std::uint64_t>(prime)){
					break;
				}
				PrimeState state{};
				if(initialize_state(kernel,prime,chunk_low,chunk_high,state)){
					states.push_back(state);
				}
			}

			std::vector<std::uint64_t> marks(
				static_cast<std::size_t>(segment_blocks),0ULL);
			std::uint64_t local_total=0;
			const bool low_aligned=(low%kernel.modulus)==0ULL;
			const bool high_aligned=(to%kernel.modulus)==0ULL;
			const std::uint64_t full_mask=kernel.full_mask_u64;
			struct DenseRun210{
				const DensePattern210* pattern=nullptr;
				std::size_t start_idx=0;
				std::uint16_t phase=0;
			};
			std::vector<DenseRun210> dense_runs;
			dense_runs.reserve(dense_patterns.size());

			for(std::uint64_t seg=seg_first;seg<seg_last;++seg){
				std::uint64_t seg_block_begin=
					block_begin+seg*segment_blocks;
				std::uint64_t seg_block_end=std::min<std::uint64_t>(
					block_end,seg_block_begin+segment_blocks);
				std::uint64_t block_count=seg_block_end-seg_block_begin;
				if(block_count==0){
					continue;
				}
				std::fill_n(marks.begin(),static_cast<std::size_t>(block_count),
							0ULL);

				dense_runs.clear();
				for(const auto&pattern : dense_patterns){
					if(pattern.start_block>=seg_block_end){
						continue;
					}
					std::uint64_t begin_block=std::max<std::uint64_t>(
						seg_block_begin,pattern.start_block);
					if(begin_block>=seg_block_end){
						continue;
					}
					std::size_t idx=
						static_cast<std::size_t>(begin_block-seg_block_begin);
					std::uint32_t phase=static_cast<std::uint32_t>(
						(begin_block*kernel.modulus)%pattern.prime);
					dense_runs.push_back(DenseRun210{
						&pattern,idx,static_cast<std::uint16_t>(phase)});
				}

				for(std::size_t idx=0;idx<block_count;++idx){
					std::uint64_t dense_mask=0ULL;
					for(auto&run : dense_runs){
						if(idx<run.start_idx){
							continue;
						}
						dense_mask|=run.pattern->block_masks[run.phase];
						run.phase=run.pattern->next_phase_block[run.phase];
					}
					marks[idx]|=dense_mask;
				}

				for(auto&state : states){
					while(state.block<seg_block_begin){
						advance_state_once(kernel,state);
					}
					while(state.block<seg_block_end){
						std::uint8_t phase=state.phase;
						std::uint8_t bit=state.bit_for_phase[phase];
						marks[static_cast<std::size_t>(state.block-
													 seg_block_begin)]|=
							(1ULL<<bit);
						advance_state_once(kernel,state);
					}
				}

				bool lower_full=(seg_block_begin>block_begin)||
								(low_aligned&&seg_block_begin==block_begin);
				bool upper_full=(seg_block_end<block_end)||
								(high_aligned&&seg_block_end==block_end);
				bool segment_full=lower_full&&upper_full;

				if(segment_full){
					for(std::uint64_t i=0;i<block_count;++i){
						std::uint64_t composite=
							marks[static_cast<std::size_t>(i)]&full_mask;
						local_total+=static_cast<std::uint64_t>(
							kernel.phase_count-popcount_u64(composite));
					}
				}else{
					for(std::uint64_t i=0;i<block_count;++i){
						std::uint64_t block_index=seg_block_begin+i;
						std::uint64_t block_base=block_index*kernel.modulus;
						std::uint64_t valid=build_valid_mask<std::uint64_t>(
							kernel,block_base,low,to);
						if(valid==0ULL){
							continue;
						}

						std::uint64_t composite=
							marks[static_cast<std::size_t>(i)]&valid;
						local_total+=static_cast<std::uint64_t>(
							popcount_u64(valid)-popcount_u64(composite));
					}
				}
			}

			total_acc.fetch_add(local_total,std::memory_order_relaxed);
		});
	}

	for(auto&worker : workers){
		worker.join();
	}
	return total_acc.load(std::memory_order_relaxed);
}

template<typename MarkT>
std::uint64_t count_with_kernel(const WheelBitmapKernel&kernel,
								std::uint64_t from,std::uint64_t to,
								unsigned threads,const SegmentConfig&config,
								const std::vector<std::uint32_t>&base_primes){
	std::uint64_t total=0;
	for(std::uint32_t prime : kernel.prefix_primes){
		std::uint64_t p64=static_cast<std::uint64_t>(prime);
		if(p64>=from&&p64<to){
			++total;
		}
	}

	std::uint64_t low=std::max<std::uint64_t>(from,kernel.min_value);
	if(to<=low){
		return total;
	}

	std::uint64_t block_begin=low/kernel.modulus;
	std::uint64_t block_end=ceil_div_u64(to,kernel.modulus);
	if(block_end<=block_begin){
		return total;
	}

	std::uint64_t segment_blocks=config.segment_bytes/sizeof(MarkT);
	if(segment_blocks==0){
		segment_blocks=1;
	}

	std::uint64_t total_blocks=block_end-block_begin;
	std::uint64_t total_segments=ceil_div_u64(total_blocks,segment_blocks);
	if(total_segments==0){
		return total;
	}

	unsigned worker_count=threads?threads:1;
	if(static_cast<std::uint64_t>(worker_count)>total_segments){
		worker_count=static_cast<unsigned>(total_segments);
	}
	if(worker_count==0){
		worker_count=1;
	}

	std::atomic<std::uint64_t> total_acc{total};
	std::vector<std::thread> workers;
	workers.reserve(worker_count);

	for(unsigned worker=0;worker<worker_count;++worker){
		workers.emplace_back([&,worker](){
			std::uint64_t seg_first=
				(total_segments*worker)/worker_count;
			std::uint64_t seg_last=
				(total_segments*(worker+1ULL))/worker_count;
			if(seg_first>=seg_last){
				return;
			}

			std::uint64_t chunk_block_begin=
				block_begin+seg_first*segment_blocks;
			std::uint64_t chunk_block_end=std::min<std::uint64_t>(
				block_end,block_begin+seg_last*segment_blocks);
			std::uint64_t chunk_low=
				chunk_block_begin*kernel.modulus;
			std::uint64_t chunk_high=
				(chunk_block_end==block_end)?to:
											 chunk_block_end*kernel.modulus;

			std::vector<PrimeState> states;
			states.reserve(base_primes.size());
			for(std::uint32_t prime : base_primes){
				if(prime<=kernel.wheel_factor_max){
					continue;
				}
				if(chunk_high>0&&
				   static_cast<std::uint64_t>(prime)>
					   (chunk_high-1ULL)/static_cast<std::uint64_t>(prime)){
					break;
				}
				PrimeState state{};
				if(initialize_state(kernel,prime,chunk_low,chunk_high,state)){
					states.push_back(state);
				}
			}

			std::vector<MarkT> marks(static_cast<std::size_t>(segment_blocks),0);
			std::uint64_t local_total=0;
			const bool low_aligned=(low%kernel.modulus)==0ULL;
			const bool high_aligned=(to%kernel.modulus)==0ULL;
			const MarkT full_mask=static_cast<MarkT>(kernel.full_mask_u64);
			for(std::uint64_t seg=seg_first;seg<seg_last;++seg){
				std::uint64_t seg_block_begin=
					block_begin+seg*segment_blocks;
				std::uint64_t seg_block_end=std::min<std::uint64_t>(
					block_end,seg_block_begin+segment_blocks);
				std::uint64_t block_count=seg_block_end-seg_block_begin;
				if(block_count==0){
					continue;
				}
				std::fill_n(marks.begin(),static_cast<std::size_t>(block_count),
							static_cast<MarkT>(0));

				for(auto&state : states){
					while(state.block<seg_block_begin){
						advance_state_once(kernel,state);
					}
					while(state.block<seg_block_end){
						std::uint8_t phase=state.phase;
						std::uint8_t bit=state.bit_for_phase[phase];
						marks[static_cast<std::size_t>(state.block-
													 seg_block_begin)]|=
							static_cast<MarkT>(static_cast<MarkT>(1)<<bit);
						advance_state_once(kernel,state);
					}
				}

				bool lower_full=(seg_block_begin>block_begin)||
								(low_aligned&&seg_block_begin==block_begin);
				bool upper_full=(seg_block_end<block_end)||
								(high_aligned&&seg_block_end==block_end);
				bool segment_full=lower_full&&upper_full;

				if(segment_full){
					for(std::uint64_t i=0;i<block_count;++i){
						MarkT composite=static_cast<MarkT>(
							marks[static_cast<std::size_t>(i)]&full_mask);
						local_total+=static_cast<std::uint64_t>(
							kernel.phase_count-popcount_mark(composite));
					}
				}else{
					for(std::uint64_t i=0;i<block_count;++i){
						std::uint64_t block_index=seg_block_begin+i;
						std::uint64_t block_base=block_index*kernel.modulus;
						MarkT valid=build_valid_mask<MarkT>(
							kernel,block_base,low,to);
						if(valid==static_cast<MarkT>(0)){
							continue;
						}

						MarkT composite=static_cast<MarkT>(
							marks[static_cast<std::size_t>(i)]&valid);
						local_total+=static_cast<std::uint64_t>(
							popcount_mark(valid)-popcount_mark(composite));
					}
				}
			}

			total_acc.fetch_add(local_total,std::memory_order_relaxed);
		});
	}

	for(auto&worker : workers){
		worker.join();
	}
	return total_acc.load(std::memory_order_relaxed);
}

} // namespace

bool supports_wheel_bitmap_count(WheelType wheel) noexcept{
	return wheel==WheelType::Mod30||wheel==WheelType::Mod210;
}

std::uint64_t count_primes_wheel_bitmap(std::uint64_t from,std::uint64_t to,
										unsigned threads,WheelType wheel,
										const SegmentConfig&config,
										const std::vector<std::uint32_t>&
											base_primes,
										bool user_segment_override){
	if(to<=from||to<2){
		return 0;
	}
	switch(wheel){
	case WheelType::Mod30:
		return count_with_kernel30_fast(kernel30(),from,to,threads,config,
										user_segment_override,base_primes);
	case WheelType::Mod210:
		return count_with_kernel210_fast(kernel210(),from,to,threads,config,
										 user_segment_override,base_primes);
	case WheelType::Mod1155:
		break;
	}
	return 0;
}

} // namespace calcprime

