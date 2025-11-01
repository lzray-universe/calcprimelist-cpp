#pragma once

#include <cstdint>
#include <vector>

namespace calcprime {

std::uint64_t meissel_count(std::uint64_t from,std::uint64_t to,const std::vector<std::uint32_t>&primes,unsigned threads=0);

bool miller_rabin_is_prime(std::uint64_t n);

}
