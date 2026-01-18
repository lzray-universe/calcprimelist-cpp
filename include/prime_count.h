#pragma once

#include<cstdint>
#include<vector>

namespace calcprime{

// Compute the number of primes in the half-open interval [from, to) using
// the Meissel–Lehmer algorithm. The provided prime table must contain all
// primes up to at least sqrt(to).
std::uint64_t meissel_count(std::uint64_t from,std::uint64_t to,
							const std::vector<std::uint32_t>&primes,
							unsigned threads=0);

// Deterministic Miller–Rabin primality check for 64-bit integers.
bool miller_rabin_is_prime(std::uint64_t n);

} // namespace calcprime
