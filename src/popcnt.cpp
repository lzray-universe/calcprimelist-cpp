#include "popcnt.h"

#include <array>
#include <cstddef>
#include <cstdint>

#if defined(__AVX512F__) && defined(__AVX512VPOPCNTDQ__)
#include <immintrin.h>
#elif defined(__AVX2__)
#include <immintrin.h>
#endif

namespace calcprime {

std::uint64_t popcount_u64(std::uint64_t x) noexcept {
#if defined(_MSC_VER) && defined(_M_X64)
    return static_cast<std::uint64_t>(__popcnt64(x));
#elif defined(__GNUC__) || defined(__clang__)
    return static_cast<std::uint64_t>(__builtin_popcountll(x));
#else

    x=x-((x>>1)&0x5555555555555555ULL);
    x=(x&0x3333333333333333ULL)+((x>>2)&0x3333333333333333ULL);
    x=(x+(x>>4))&0x0F0F0F0F0F0F0F0FULL;
    return (x*0x0101010101010101ULL)>>56;
#endif
}

std::uint64_t count_zero_bits(const std::uint64_t*bits,std::size_t bit_count) noexcept {
    std::uint64_t total=0;
    std::size_t full_words=bit_count/64;
    std::size_t rem_bits=bit_count%64;

#if defined(__AVX512F__) && defined(__AVX512VPOPCNTDQ__)
    constexpr std::size_t stride=8;
    std::size_t i=0;
    for(;i+stride<=full_words;i+=stride) {
        __m512i data=_mm512_loadu_si512(reinterpret_cast<const void*>(bits+i));
        __m512i pop=_mm512_popcnt_epi64(data);
        alignas(64) std::array<std::uint64_t,stride>buf;
        _mm512_store_si512(reinterpret_cast<void*>(buf.data()),pop);
        for(std::size_t j=0;j<stride;++j) {
            total+=64-buf[j];
        }
    }
    for(;i<full_words;++i) {
        total+=64-popcount_u64(bits[i]);
    }
#elif defined(__AVX2__)
    constexpr std::size_t stride=4;
    const __m256i low_mask=_mm256_set1_epi8(0x0F);
    const __m256i nibble_popcnt=_mm256_setr_epi8(0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4,
        0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4);
    const __m256i zero=_mm256_setzero_si256();
    std::uint64_t ones=0;
    std::size_t i=0;
    for(;i+stride<=full_words;i+=stride) {
        __m256i data=_mm256_loadu_si256(reinterpret_cast<const __m256i*>(bits+i));
        __m256i lo=_mm256_and_si256(data,low_mask);
        __m256i hi=_mm256_and_si256(_mm256_srli_epi16(data,4),low_mask);
        __m256i popcnt=_mm256_add_epi8(_mm256_shuffle_epi8(nibble_popcnt,lo),
            _mm256_shuffle_epi8(nibble_popcnt,hi));
        __m256i sad=_mm256_sad_epu8(popcnt,zero);
        alignas(32) std::array<std::uint64_t,stride>buf{};
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(buf.data()),sad);
        ones+=buf[0]+buf[1]+buf[2]+buf[3];
    }
    for(;i<full_words;++i) {
        ones+=popcount_u64(bits[i]);
    }
    total+=full_words*64-ones;
#else
    for(std::size_t i=0;i<full_words;++i) {
        total+=64-popcount_u64(bits[i]);
    }
#endif

    if(rem_bits) {
        std::uint64_t mask=(rem_bits==64)?~0ULL:((1ULL<<rem_bits)-1);
        total+=rem_bits-popcount_u64(bits[full_words]&mask);
    }
    return total;
}

}
