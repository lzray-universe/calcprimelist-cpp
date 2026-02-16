#include "popcnt.h"

#include<cstddef>
#include<cstdint>

#if defined(__AVX512F__)&&defined(__AVX512VPOPCNTDQ__)
#include<immintrin.h>
#elif defined(__AVX2__)
#include<immintrin.h>
#endif

namespace calcprime{
namespace{

#if defined(__AVX512F__)&&defined(__AVX512VPOPCNTDQ__)

std::uint64_t popcount_words_u64_avx512(const std::uint64_t*words,
										std::size_t word_count) noexcept{
	constexpr std::size_t kStride=8;
	std::size_t i=0;
	__m512i acc=_mm512_setzero_si512();
	for(;i+kStride<=word_count;i+=kStride){
		__m512i data=
			_mm512_loadu_si512(reinterpret_cast<const void*>(words+i));
		__m512i pop=_mm512_popcnt_epi64(data);
		acc=_mm512_add_epi64(acc,pop);
	}
	alignas(64) std::uint64_t lanes[kStride];
	_mm512_store_si512(reinterpret_cast<void*>(lanes),acc);
	std::uint64_t total=0;
	for(std::size_t lane=0;lane<kStride;++lane){
		total+=lanes[lane];
	}
	for(;i<word_count;++i){
		total+=popcount_u64(words[i]);
	}
	return total;
}

std::uint64_t popcount_words_u64_masked_avx512(const std::uint64_t*words,
											   std::size_t word_count,
											   std::uint64_t mask) noexcept{
	if(mask==~0ULL){
		return popcount_words_u64_avx512(words,word_count);
	}
	constexpr std::size_t kStride=8;
	const __m512i mask_vec=
		_mm512_set1_epi64(static_cast<long long>(mask));
	std::size_t i=0;
	__m512i acc=_mm512_setzero_si512();
	for(;i+kStride<=word_count;i+=kStride){
		__m512i data=
			_mm512_loadu_si512(reinterpret_cast<const void*>(words+i));
		data=_mm512_and_si512(data,mask_vec);
		__m512i pop=_mm512_popcnt_epi64(data);
		acc=_mm512_add_epi64(acc,pop);
	}
	alignas(64) std::uint64_t lanes[kStride];
	_mm512_store_si512(reinterpret_cast<void*>(lanes),acc);
	std::uint64_t total=0;
	for(std::size_t lane=0;lane<kStride;++lane){
		total+=lanes[lane];
	}
	for(;i<word_count;++i){
		total+=popcount_u64(words[i]&mask);
	}
	return total;
}

#elif defined(__AVX2__)

std::uint64_t popcount_words_u64_avx2(const std::uint64_t*words,
									  std::size_t word_count) noexcept{
	const __m256i nibble_mask=_mm256_set1_epi8(0x0F);
	const __m256i lookup=_mm256_setr_epi8(
		0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4,
		0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4);
	const __m256i zero=_mm256_setzero_si256();
	constexpr std::size_t kStride=4;
	constexpr std::size_t kUnroll=16;

	std::size_t i=0;
	__m256i acc0=_mm256_setzero_si256();
	__m256i acc1=_mm256_setzero_si256();
	__m256i acc2=_mm256_setzero_si256();
	__m256i acc3=_mm256_setzero_si256();

	auto accumulate_lane=[&](const std::uint64_t*ptr,__m256i&accumulator){
		__m256i data=
			_mm256_loadu_si256(reinterpret_cast<const __m256i*>(ptr));
		__m256i lo=_mm256_and_si256(data,nibble_mask);
		__m256i hi=
			_mm256_and_si256(_mm256_srli_epi16(data,4),nibble_mask);
		__m256i pop8=
			_mm256_add_epi8(_mm256_shuffle_epi8(lookup,lo),
							_mm256_shuffle_epi8(lookup,hi));
		__m256i lane_sum=_mm256_sad_epu8(pop8,zero);
		accumulator=_mm256_add_epi64(accumulator,lane_sum);
	};

	for(;i+kUnroll<=word_count;i+=kUnroll){
		accumulate_lane(words+i+0,acc0);
		accumulate_lane(words+i+4,acc1);
		accumulate_lane(words+i+8,acc2);
		accumulate_lane(words+i+12,acc3);
	}
	for(;i+kStride<=word_count;i+=kStride){
		accumulate_lane(words+i,acc0);
	}

	alignas(32) std::uint64_t lanes[kStride];
	__m256i acc=_mm256_add_epi64(_mm256_add_epi64(acc0,acc1),
								 _mm256_add_epi64(acc2,acc3));
	_mm256_store_si256(reinterpret_cast<__m256i*>(lanes),acc);
	std::uint64_t total=lanes[0]+lanes[1]+lanes[2]+lanes[3];
	for(;i<word_count;++i){
		total+=popcount_u64(words[i]);
	}
	return total;
}

std::uint64_t popcount_words_u64_masked_avx2(const std::uint64_t*words,
											 std::size_t word_count,
											 std::uint64_t mask) noexcept{
	if(mask==~0ULL){
		return popcount_words_u64_avx2(words,word_count);
	}
	const __m256i nibble_mask=_mm256_set1_epi8(0x0F);
	const __m256i lookup=_mm256_setr_epi8(
		0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4,
		0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4);
	const __m256i zero=_mm256_setzero_si256();
	const __m256i mask_vec=
		_mm256_set1_epi64x(static_cast<long long>(mask));
	constexpr std::size_t kStride=4;
	constexpr std::size_t kUnroll=16;

	std::size_t i=0;
	__m256i acc0=_mm256_setzero_si256();
	__m256i acc1=_mm256_setzero_si256();
	__m256i acc2=_mm256_setzero_si256();
	__m256i acc3=_mm256_setzero_si256();

	auto accumulate_lane=[&](const std::uint64_t*ptr,__m256i&accumulator){
		__m256i data=
			_mm256_loadu_si256(reinterpret_cast<const __m256i*>(ptr));
		data=_mm256_and_si256(data,mask_vec);
		__m256i lo=_mm256_and_si256(data,nibble_mask);
		__m256i hi=
			_mm256_and_si256(_mm256_srli_epi16(data,4),nibble_mask);
		__m256i pop8=
			_mm256_add_epi8(_mm256_shuffle_epi8(lookup,lo),
							_mm256_shuffle_epi8(lookup,hi));
		__m256i lane_sum=_mm256_sad_epu8(pop8,zero);
		accumulator=_mm256_add_epi64(accumulator,lane_sum);
	};

	for(;i+kUnroll<=word_count;i+=kUnroll){
		accumulate_lane(words+i+0,acc0);
		accumulate_lane(words+i+4,acc1);
		accumulate_lane(words+i+8,acc2);
		accumulate_lane(words+i+12,acc3);
	}
	for(;i+kStride<=word_count;i+=kStride){
		accumulate_lane(words+i,acc0);
	}

	alignas(32) std::uint64_t lanes[kStride];
	__m256i acc=_mm256_add_epi64(_mm256_add_epi64(acc0,acc1),
								 _mm256_add_epi64(acc2,acc3));
	_mm256_store_si256(reinterpret_cast<__m256i*>(lanes),acc);
	std::uint64_t total=lanes[0]+lanes[1]+lanes[2]+lanes[3];
	for(;i<word_count;++i){
		total+=popcount_u64(words[i]&mask);
	}
	return total;
}

#endif

} // namespace

std::uint64_t popcount_u64(std::uint64_t x) noexcept{
#if defined(_MSC_VER)&&defined(_M_X64)
	return static_cast<std::uint64_t>(__popcnt64(x));
#elif defined(__GNUC__)||defined(__clang__)
	return static_cast<std::uint64_t>(__builtin_popcountll(x));
#else
	// fallback to software popcount
	x=x-((x>>1)&0x5555555555555555ULL);
	x=(x&0x3333333333333333ULL)+((x>>2)&0x3333333333333333ULL);
	x=(x+(x>>4))&0x0F0F0F0F0F0F0F0FULL;
	return (x*0x0101010101010101ULL)>>56;
#endif
}

std::uint64_t popcount_words_u64(const std::uint64_t*words,
								 std::size_t word_count) noexcept{
	if(words==nullptr||word_count==0){
		return 0;
	}
#if defined(__AVX512F__)&&defined(__AVX512VPOPCNTDQ__)
	return popcount_words_u64_avx512(words,word_count);
#elif defined(__AVX2__)
	return popcount_words_u64_avx2(words,word_count);
#else
	std::uint64_t total=0;
	for(std::size_t i=0;i<word_count;++i){
		total+=popcount_u64(words[i]);
	}
	return total;
#endif
}

std::uint64_t popcount_words_u64_masked(const std::uint64_t*words,
										std::size_t word_count,
										std::uint64_t mask) noexcept{
	if(words==nullptr||word_count==0){
		return 0;
	}
#if defined(__AVX512F__)&&defined(__AVX512VPOPCNTDQ__)
	return popcount_words_u64_masked_avx512(words,word_count,mask);
#elif defined(__AVX2__)
	return popcount_words_u64_masked_avx2(words,word_count,mask);
#else
	std::uint64_t total=0;
	for(std::size_t i=0;i<word_count;++i){
		total+=popcount_u64(words[i]&mask);
	}
	return total;
#endif
}

std::uint64_t count_zero_bits(const std::uint64_t*bits,
							  std::size_t bit_count) noexcept{
	if(bits==nullptr||bit_count==0){
		return 0;
	}
	std::size_t full_words=bit_count/64;
	std::size_t rem_bits=bit_count%64;

	std::uint64_t ones=popcount_words_u64(bits,full_words);
	std::uint64_t total=full_words*64ULL-ones;

	if(rem_bits){
		std::uint64_t mask=(1ULL<<rem_bits)-1ULL;
		total+=rem_bits-popcount_words_u64_masked(bits+full_words,1,mask);
	}
	return total;
}

} // namespace calcprime
