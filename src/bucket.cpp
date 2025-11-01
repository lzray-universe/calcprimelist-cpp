#include "bucket.h"

#include <algorithm>
#include <iterator>

namespace calcprime {

BucketRing::BucketRing() : base_segment_(0),mask_(0) {}

void BucketRing::reset(std::uint64_t start_segment) {
    base_segment_=start_segment;
    mask_=0;
    buckets_.clear();
}

void BucketRing::ensure_capacity(std::uint64_t segment) {
    if(segment<base_segment_) {
        return;
    }
    if(buckets_.empty()) {
        buckets_.resize(1024);
        mask_=buckets_.size()-1;
        return;
    }
    while((segment-base_segment_)>mask_) {
        rehash(buckets_.size()*2);
    }
}

void BucketRing::rehash(std::size_t new_size) {
    if((new_size&(new_size-1))!=0) {

        std::size_t pow2=1;
        while(pow2<new_size) {
            pow2<<=1;
        }
        new_size=pow2;
    }
    std::vector<std::vector<BucketEntry>>new_buckets(new_size);
    std::size_t new_mask=new_size-1;
    for(auto&bucket : buckets_) {
        for(auto&entry : bucket) {
            new_buckets[entry.next_index&new_mask].push_back(entry);
        }
    }
    buckets_=std::move(new_buckets);
    mask_=new_mask;
}

void BucketRing::push(std::uint64_t segment,BucketEntry entry) {
    ensure_capacity(segment);
    if(buckets_.empty()) {
        buckets_.resize(1024);
        mask_=buckets_.size()-1;
    }
    buckets_[segment&mask_].push_back(entry);
}

std::vector<BucketEntry>BucketRing::take(std::uint64_t segment) {
    ensure_capacity(segment);
    std::vector<BucketEntry>hits;
    if(buckets_.empty()) {
        return hits;
    }
    auto&bucket=buckets_[segment&mask_];
    auto first_hit=std::partition(bucket.begin(),bucket.end(),
        [segment](const BucketEntry&entry) { return entry.next_index!=segment;});

    if(first_hit!=bucket.end()) {
        hits.reserve(static_cast<std::size_t>(bucket.end()-first_hit));
        hits.insert(hits.end(),std::make_move_iterator(first_hit),
                    std::make_move_iterator(bucket.end()));
        bucket.erase(first_hit,bucket.end());
    }
    if(segment>=base_segment_) {
        base_segment_=segment+1;
    }
    return hits;
}

}
