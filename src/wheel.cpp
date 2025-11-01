#include "wheel.h"

#include <algorithm>
#include <bit>
#include <numeric>

namespace calcprime {
namespace {

SmallPrimePattern build_small_pattern(std::uint32_t prime) {
    SmallPrimePattern pattern{};
    pattern.prime=prime;
    pattern.phase_count=prime;
    pattern.masks.resize(prime);
    pattern.next_phase.resize(prime);
    std::uint32_t word_stride=static_cast<std::uint32_t>(128%prime);
    pattern.word_stride=word_stride;
    std::uint32_t inv2=(prime+1)/2;
    for(std::size_t bit=0;bit<pattern.start_phase.size();++bit) {
        std::uint32_t twice=static_cast<std::uint32_t>(bit<<1);
        while(twice>=prime) {
            twice-=prime;
        }
        std::uint32_t phase=(prime-twice);
        if(phase==prime) {
            phase=0;
        }
        pattern.start_phase[bit]=static_cast<std::uint8_t>(phase);
    }
    for(std::uint32_t residue=0;residue<prime;++residue) {
        std::uint64_t mask=0;
        std::uint32_t offset=static_cast<std::uint32_t>(((prime-residue)%prime)*static_cast<std::uint64_t>(inv2)%prime);
        while(offset<64) {
            mask|=(1ULL<<offset);
            offset+=prime;
        }
        pattern.masks[residue]=mask;
        pattern.next_phase[residue]=static_cast<std::uint32_t>((residue+word_stride)%prime);
    }
    return pattern;
}

Wheel build_wheel(std::uint32_t modulus,WheelType type) {
    Wheel wheel;
    wheel.type=type;
    wheel.modulus=modulus;
    wheel.allowed.assign(modulus,0);

    for(std::uint32_t r=0;r<modulus;++r) {
        if(std::gcd(r,modulus)==1) {
            wheel.allowed[r]=1;
            wheel.residues.push_back(static_cast<std::uint16_t>(r));
        }
    }
    if(!wheel.residues.empty()) {
        wheel.steps.reserve(wheel.residues.size());
        for(std::size_t i=0;i<wheel.residues.size();++i) {
            std::uint32_t a=wheel.residues[i];
            std::uint32_t b=wheel.residues[(i+1)%wheel.residues.size()];
            std::uint32_t step=(b+modulus-a)%modulus;
            if(step==0) {
                step=modulus;
            }
            wheel.steps.push_back(static_cast<std::uint16_t>(step));
        }
    }

    static const std::uint32_t kSmallPrimes[]={3,5,7,11,13,17,19,
                                                 23,29,31,37,41,43,47};
    std::uint32_t small_limit=29u;
    switch(type) {
    case WheelType::Mod30:
        small_limit=29u;
        break;
    case WheelType::Mod210:
        small_limit=47u;
        break;
    case WheelType::Mod1155:
        small_limit=47u;
        break;
    }
    for(std::uint32_t prime : kSmallPrimes) {
        if(prime>small_limit) {
            break;
        }
        wheel.small_patterns.push_back(build_small_pattern(prime));
    }
    return wheel;
}

}

const Wheel&get_wheel(WheelType type) {
    static const Wheel wheel30=build_wheel(30,WheelType::Mod30);
    static const Wheel wheel210=build_wheel(210,WheelType::Mod210);
    static const Wheel wheel1155=build_wheel(1155,WheelType::Mod1155);
    switch(type) {
    case WheelType::Mod30:
        return wheel30;
    case WheelType::Mod210:
        return wheel210;
    case WheelType::Mod1155:
        return wheel1155;
    }
    return wheel30;
}

void Wheel::apply_presieve(std::uint64_t start_value,std::size_t bit_count,std::uint64_t*bits) const {
    if(allowed.empty()) {
        return;
    }
    std::uint32_t rem=static_cast<std::uint32_t>(start_value%modulus);
    for(std::size_t idx=0;idx<bit_count;++idx) {
        if(!allowed[rem]) {
            std::size_t word=idx/64;
            std::size_t bit=idx%64;
            bits[word]|=(1ULL<<bit);
        }
        rem+=2;
        if(rem>=modulus) {
            rem-=modulus;
        }
    }
}

}
