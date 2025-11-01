#include "prime_count.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <future>
#include <map>
#include <mutex>
#include <thread>
#include <utility>
#include <vector>

namespace calcprime {
namespace {

std::uint64_t integer_sqrt(std::uint64_t n) {
    long double root_ld=std::sqrt(static_cast<long double>(n));
    std::uint64_t root=static_cast<std::uint64_t>(root_ld);
    auto square_le=[n](std::uint64_t v) {
        if(v==0) {
            return true;
        }
        return v<=n/v;
    };
    while(square_le(root+1)) {
        ++root;
    }
    while(!square_le(root)) {
        --root;
    }
    return root;
}

std::uint64_t integer_cuberoot(std::uint64_t n) {
    long double root_ld=std::cbrt(static_cast<long double>(n));
    std::uint64_t root=static_cast<std::uint64_t>(root_ld);
    auto cube_le=[n](std::uint64_t v) {
        if(v==0) {
            return true;
        }
        if(v>n/v) {
            return false;
        }
        std::uint64_t q=n/v;
        return v<=q/v;
    };
    while(cube_le(root+1)) {
        ++root;
    }
    while(!cube_le(root)) {
        --root;
    }
    return root;
}

std::uint64_t integer_fourth_root(std::uint64_t n) {
    if(n==0) {
        return 0;
    }
    std::uint64_t root=integer_sqrt(integer_sqrt(n));
    auto fourth_le=[n](std::uint64_t v) {
        if(v==0) {
            return true;
        }
        if(v>n/v) {
            return false;
        }
        std::uint64_t q=n/v;
        if(v>q/v) {
            return false;
        }
        q=q/v;
        return v<=q/v;
    };
    while(fourth_le(root+1)) {
        ++root;
    }
    while(!fourth_le(root)) {
        --root;
    }
    return root;
}

class MeisselCalculator {
public:
    explicit MeisselCalculator(const std::vector<std::uint32_t>&primes)
        : primes_(primes),max_prime_(primes.empty() ? 0 : primes.back()) {}

    std::uint64_t pi(std::uint64_t n,unsigned threads=1) {
        if(n<2) {
            return 0;
        }
        if(primes_.empty()) {
            return small_pi(n);
        }
        if(n<=max_prime_) {
            auto it=std::upper_bound(primes_.begin(),primes_.end(),static_cast<std::uint32_t>(n));
            return static_cast<std::uint64_t>(std::distance(primes_.begin(),it));
        }
        {
            std::lock_guard<std::mutex>lock(pi_mutex_);
            auto cached=pi_cache_.find(n);
            if(cached!=pi_cache_.end()) {
                return cached->second;
            }
        }

        std::uint64_t a=pi(integer_fourth_root(n),1);
        std::uint64_t b=pi(integer_sqrt(n),1);
        std::uint64_t c=pi(integer_cuberoot(n),1);

        std::uint64_t result=phi(n,static_cast<std::size_t>(a));
        if(b+a>=2) {
            std::uint64_t left=b+a-2;
            std::uint64_t right=b-a+1;
            result+=(left*right)/2;
        }

        std::uint64_t effective_b=std::min<std::uint64_t>(b,static_cast<std::uint64_t>(primes_.size()));
        std::uint64_t iteration_count=0;
        if(effective_b>a) {
            iteration_count=effective_b-a;
        }

        auto compute_range=[this,n,c](std::uint64_t start,std::uint64_t end) {
            std::uint64_t subtotal=0;
            for(std::uint64_t i=start;i<end;++i) {
                std::uint64_t index=i-1;
                if(index>=primes_.size()) {
                    break;
                }
                std::uint64_t p=primes_[static_cast<std::size_t>(index)];
                std::uint64_t w=n/p;
                subtotal+=this->pi(w,1);
                if(i<=c) {
                    std::uint64_t limit=this->pi(integer_sqrt(w),1);
                    for(std::uint64_t j=i;j<=limit;++j) {
                        std::uint64_t j_index=j-1;
                        if(j_index>=primes_.size()) {
                            break;
                        }
                        std::uint64_t pj=primes_[static_cast<std::size_t>(j_index)];
                        subtotal+=this->pi(w/pj,1)-(j-1);
                    }
                }
            }
            return subtotal;
        };

        if(iteration_count>0) {
            if(threads<=1||iteration_count==1) {
                result-=compute_range(a+1,effective_b+1);
            } else {
                unsigned worker_count=std::min<std::uint64_t>(threads,iteration_count);
                if(worker_count==0) {
                    worker_count=1;
                }
                std::vector<std::future<std::uint64_t>>futures;
                futures.reserve(worker_count);
                std::uint64_t chunk=iteration_count/worker_count;
                std::uint64_t remainder=iteration_count%worker_count;
                std::uint64_t current=a+1;
                for(unsigned w=0;w<worker_count;++w) {
                    std::uint64_t size=chunk+(w<remainder ? 1 : 0);
                    if(size==0) {
                        continue;
                    }
                    std::uint64_t chunk_start=current;
                    std::uint64_t chunk_end=chunk_start+size;
                    current=chunk_end;
                    futures.emplace_back(std::async(std::launch::async,
                                                    [compute_range,chunk_start,chunk_end]() {
                                                        return compute_range(chunk_start,chunk_end);
                                                    }));
                }
                std::uint64_t subtract_total=0;
                for(auto&fut : futures) {
                    subtract_total+=fut.get();
                }
                result-=subtract_total;
            }
        }

        {
            std::lock_guard<std::mutex>lock(pi_mutex_);
            auto inserted=pi_cache_.emplace(n,result);
            if(!inserted.second) {
                result=inserted.first->second;
            }
        }
        return result;
    }

private:
    std::uint64_t phi(std::uint64_t x,std::size_t s) {
        if(s==0) {
            return x;
        }
        if(s==1) {
            return (x+1)>>1;
        }
        if(s>primes_.size()) {
            return phi(x,primes_.size());
        }
        auto key=std::make_pair(x,s);
        {
            std::lock_guard<std::mutex>lock(phi_mutex_);
            auto cached=phi_cache_.find(key);
            if(cached!=phi_cache_.end()) {
                return cached->second;
            }
        }
        std::uint64_t result=phi(x,s-1);
        std::uint32_t p=primes_[s-1];
        result-=phi(x/p,s-1);
        {
            std::lock_guard<std::mutex>lock(phi_mutex_);
            phi_cache_.emplace(key,result);
        }
        return result;
    }

    static std::uint64_t small_pi(std::uint64_t n) {
        static constexpr std::array<std::uint32_t,12>small_primes{
            2,3,5,7,11,13,17,19,23,29,31,37};
        auto it=std::upper_bound(small_primes.begin(),small_primes.end(),static_cast<std::uint32_t>(n));
        return static_cast<std::uint64_t>(std::distance(small_primes.begin(),it));
    }

    const std::vector<std::uint32_t>&primes_;
    std::uint64_t max_prime_;
    std::map<std::pair<std::uint64_t,std::size_t>,std::uint64_t>phi_cache_;
    std::map<std::uint64_t,std::uint64_t>pi_cache_;
    mutable std::mutex phi_mutex_;
    mutable std::mutex pi_mutex_;
};

std::uint64_t count_small_range(std::uint64_t from,std::uint64_t to) {
    auto count_up_to=[](std::uint64_t bound) {
        if(bound<2) {
            return std::uint64_t{0};
        }
        static constexpr std::array<std::uint32_t,12>sp{
            2,3,5,7,11,13,17,19,23,29,31,37};
        auto it=std::upper_bound(sp.begin(),sp.end(),static_cast<std::uint32_t>(bound));
        return static_cast<std::uint64_t>(std::distance(sp.begin(),it));
    };
    std::uint64_t upper=to==0 ? 0 : count_up_to(to-1);
    std::uint64_t lower=from==0 ? 0 : count_up_to(from-1);
    return upper>=lower ? upper-lower : 0;
}

}

std::uint64_t meissel_count(std::uint64_t from,std::uint64_t to,const std::vector<std::uint32_t>&primes,unsigned threads) {
    if(to<=from) {
        return 0;
    }
    if(primes.empty()) {
        return count_small_range(from,to);
    }
    unsigned effective_threads=threads;
    if(effective_threads==0) {
        effective_threads=std::thread::hardware_concurrency();
    }
    if(effective_threads==0) {
        effective_threads=1;
    }
    MeisselCalculator calc(primes);
    std::uint64_t upper=calc.pi(to-1,effective_threads);
    std::uint64_t lower=from==0?0:calc.pi(from-1,effective_threads);
    return upper>=lower?upper-lower:0;
}

namespace {
std::uint64_t mul_mod(std::uint64_t a,std::uint64_t b,std::uint64_t mod) {
    std::uint64_t result=0;
    while(b>0) {
        if(b&1ULL) {
            result=(result+a)%mod;
        }
        a=(a<<1U)%mod;
        b>>=1U;
    }
    return result;
}

std::uint64_t mod_pow(std::uint64_t base,std::uint64_t exp,std::uint64_t mod) {
    std::uint64_t result=1%mod;
    std::uint64_t b=base%mod;
    while(exp>0) {
        if(exp&1ULL) {
            result=mul_mod(result,b,mod);
        }
        b=mul_mod(b,b,mod);
        exp>>=1ULL;
    }
    return static_cast<std::uint64_t>(result);
}

bool check_composite(std::uint64_t n,std::uint64_t a,std::uint64_t d,unsigned r) {
    std::uint64_t x=mod_pow(a,d,n);
    if(x==1||x==n-1) {
        return false;
    }
    for(unsigned i=1;i<r;++i) {
        x=static_cast<std::uint64_t>(mul_mod(x,x,n));
        if(x==n-1) {
            return false;
        }
    }
    return true;
}
}

bool miller_rabin_is_prime(std::uint64_t n) {
    if(n<2) {
        return false;
    }
    static constexpr std::array<std::uint64_t,12>bases{
        2,3,5,7,11,13,17,19,23,29,31,37};
    for(std::uint64_t p : {2ULL,3ULL,5ULL,7ULL,11ULL,13ULL,17ULL,19ULL,23ULL,29ULL,31ULL,37ULL}) {
        if(n==p) {
            return true;
        }
        if(n%p==0) {
            return false;
        }
    }
    std::uint64_t d=n-1;
    unsigned r=0;
    while((d&1ULL)==0) {
        d>>=1ULL;
        ++r;
    }
    for(std::uint64_t a : bases) {
        if(a%n==0) {
            continue;
        }
        if(check_composite(n,a,d,r)) {
            return false;
        }
    }
    return true;
}

}
