#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace calcprime {

enum class WheelType {
    Mod30,
    Mod210,
    Mod1155,
};

struct SmallPrimePattern {
    std::uint32_t prime;
    std::uint32_t phase_count;
    std::uint32_t word_stride;
    std::vector<std::uint64_t>masks;
    std::vector<std::uint32_t>next_phase;
    std::array<std::uint8_t,64>start_phase;
};

struct Wheel {
    WheelType type;
    std::uint32_t modulus;
    std::vector<std::uint8_t>allowed;
    std::vector<std::uint16_t>residues;
    std::vector<std::uint16_t>steps;
    std::vector<SmallPrimePattern>small_patterns;

    void apply_presieve(std::uint64_t start_value,std::size_t bit_count,std::uint64_t*bits) const;
};

const Wheel&get_wheel(WheelType type);

}
