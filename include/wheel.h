#pragma once

#include<array>
#include<cstddef>
#include<cstdint>
#include<vector>

namespace calcprime{

enum class WheelType{
	Mod30,
	Mod210,
	Mod1155,
};

struct SmallPrimePattern{
	std::uint32_t prime;
	std::uint32_t phase_count;
	std::uint32_t word_stride;
	std::vector<std::uint64_t> masks; // masks[phase] covers one 64-bit word
	std::vector<std::uint32_t> next_phase;
	std::array<std::uint8_t,64> start_phase;
};

struct Wheel{
	WheelType type;
	std::uint32_t modulus;
	std::uint32_t presieve_modulus;
	std::vector<std::uint8_t> allowed; // allowed residues modulo modulus
	std::vector<std::uint16_t> residues;
	std::vector<std::uint16_t> steps;
	std::vector<std::uint16_t> presieved_primes;
	std::vector<SmallPrimePattern> small_patterns;
	std::vector<std::uint64_t> presieve_word_masks;
	std::vector<std::uint16_t> presieve_next_phase;
	std::vector<std::uint64_t> presieve_block_masks; // flattened [phase][word]
	std::vector<std::uint16_t> presieve_next_block_phase;

	void fill_presieve(std::uint64_t start_value,std::size_t bit_count,
					   std::uint64_t*bits) const;
	void apply_presieve(std::uint64_t start_value,std::size_t bit_count,
						std::uint64_t*bits) const;
};

const Wheel&get_wheel(WheelType type);

} // namespace calcprime
