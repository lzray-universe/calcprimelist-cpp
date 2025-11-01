#pragma once

#include <cstddef>
#include <cstdint>

namespace calcprime {

std::uint64_t popcount_u64(std::uint64_t x) noexcept;
std::uint64_t count_zero_bits(const std::uint64_t*bits,std::size_t bit_count) noexcept;

}
