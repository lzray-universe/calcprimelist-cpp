#include "base_sieve.h"

#include <cmath>
#include <vector>

namespace calcprime {

std::vector<std::uint32_t>simple_sieve(std::uint64_t limit) {
    if(limit<2) {
        return {};
    }
    std::uint64_t root=static_cast<std::uint64_t>(std::sqrt(static_cast<long double>(limit)));
    std::uint64_t max=std::max(limit,root+1);
    std::uint64_t size=(max+1)/2;
    std::vector<bool>is_composite(size,false);
    std::uint64_t bound=(static_cast<std::uint64_t>(std::sqrt(static_cast<long double>(max)))+1)/2;
    for(std::uint64_t i=1;i<=bound;++i) {
        if(!is_composite[i]) {
            std::uint64_t p=2*i+1;
            std::uint64_t start=(p*p)/2;
            for(std::uint64_t j=start;j<size;j+=p) {
                is_composite[j]=true;
            }
        }
    }
    std::vector<std::uint32_t>primes;
    primes.push_back(2);
    for(std::uint64_t i=1;i<size&&(2*i+1)<=limit;++i) {
        if(!is_composite[i]) {
            primes.push_back(static_cast<std::uint32_t>(2*i+1));
        }
    }
    return primes;
}

}
