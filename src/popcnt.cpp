#include "popcnt.h"

#include<array>
#include<cstddef>
#include<cstdint>

#if defined(__AVX512F__)&&defined(__AVX512VPOPCNTDQ__)
#include<immintrin.h>
#elif defined(__AVX2__)
#include<immintrin.h>
#endif

namespace calcprime{

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

std::uint64_t count_zero_bits(const std::uint64_t*bits,
							  std::size_t bit_count) noexcept{
	std::uint64_t total=0;
	std::size_t full_words=bit_count/64;
	std::size_t rem_bits=bit_count%64;

#if defined(__AVX512F__)&&defined(__AVX512VPOPCNTDQ__)
	constexpr std::size_t stride=8;
	std::size_t i=0;
	for(;i+stride<=full_words;i+=stride){
		__m512i data=_mm512_loadu_si512(reinterpret_cast<const void*>(bits+i));
		__m512i pop=_mm512_popcnt_epi64(data);
		alignas(64) std::array<std::uint64_t,stride> buf;
		_mm512_store_si512(reinterpret_cast<void*>(buf.data()),pop);
		for(std::size_t j=0;j<stride;++j){
			total+=64-buf[j];
		}
	}
	for(;i<full_words;++i){
		total+=64-popcount_u64(bits[i]);
	}
#elif defined(__AVX2__)
	std::uint64_t ones=0;
	std::size_t i=0;
	for(;i+4<=full_words;i+=4){
		ones+=popcount_u64(bits[i]);
		ones+=popcount_u64(bits[i+1]);
		ones+=popcount_u64(bits[i+2]);
		ones+=popcount_u64(bits[i+3]);
	}
	for(;i<full_words;++i){
		ones+=popcount_u64(bits[i]);
	}
	total+=full_words*64-ones;
#else
	for(std::size_t i=0;i<full_words;++i){
		total+=64-popcount_u64(bits[i]);
	}
#endif

	if(rem_bits){
		std::uint64_t mask=(rem_bits==64)?~0ULL:((1ULL<<rem_bits)-1);
		total+=rem_bits-popcount_u64(bits[full_words]&mask);
	}
	return total;
}

} // namespace calcprime
