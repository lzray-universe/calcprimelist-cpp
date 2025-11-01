#pragma once

#include <cstdint>
#include <vector>

namespace calcprime {

struct LargePrimeState;

struct BucketEntry {
    std::uint32_t prime;
    std::uint64_t next_index;
    std::uint64_t offset;
    std::uint64_t value;
    LargePrimeState*owner;
};

class BucketRing {
public:
    BucketRing();
    void reset(std::uint64_t start_segment);
    void push(std::uint64_t segment,BucketEntry entry);
    std::vector<BucketEntry>take(std::uint64_t segment);

private:
    void ensure_capacity(std::uint64_t segment);
    void rehash(std::size_t new_size);

    std::uint64_t base_segment_;
    std::size_t mask_;
    std::vector<std::vector<BucketEntry>>buckets_;
};

}
