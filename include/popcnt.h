#pragma once

#include<cstddef>
#include<cstdint>

namespace calcprime{

std::uint64_t popcount_u64(std::uint64_t x) noexcept;
std::uint64_t popcount_words_u64(const std::uint64_t*words,
								 std::size_t word_count) noexcept;
std::uint64_t popcount_words_u64_masked(const std::uint64_t*words,
										std::size_t word_count,
										std::uint64_t mask) noexcept;
std::uint64_t count_zero_bits(const std::uint64_t*bits,
							  std::size_t bit_count) noexcept;

} // namespace calcprime
