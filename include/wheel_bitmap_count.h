#pragma once

#include "segmenter.h"
#include "wheel.h"

#include<cstdint>
#include<vector>

namespace calcprime{

bool supports_wheel_bitmap_count(WheelType wheel) noexcept;

std::uint64_t count_primes_wheel_bitmap(std::uint64_t from,std::uint64_t to,
										unsigned threads,WheelType wheel,
										const SegmentConfig&config,
										const std::vector<std::uint32_t>&
											base_primes,
										bool user_segment_override);

} // namespace calcprime
