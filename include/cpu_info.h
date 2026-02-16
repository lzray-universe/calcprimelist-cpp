#pragma once

#include<cstddef>
#include<cstdint>
#include<string>

namespace calcprime{

struct CpuInfo{
	unsigned logical_cpus=1;
	unsigned physical_cpus=1;
	unsigned performance_logical_cpus=1;
	unsigned efficiency_logical_cpus=0;
	std::size_t l1_data_bytes=32*1024;
	std::size_t l2_bytes=1024*1024;
	std::size_t performance_l1_data_bytes=32*1024;
	std::size_t performance_l2_bytes=1024*1024;
	std::size_t efficiency_l1_data_bytes=32*1024;
	std::size_t efficiency_l2_bytes=1024*1024;
	std::size_t l2_total_bytes=0;
	bool has_smt=false;
	bool has_hybrid=false;
};

enum class CoreSchedulingMode{
	Auto,
	BigOnly,
	AllCores,
	Legacy,
};

CpuInfo detect_cpu_info();

unsigned effective_thread_count(const CpuInfo&info);

unsigned choose_thread_count(const CpuInfo&info,unsigned requested_threads,
							 std::uint64_t range_span,
							 CoreSchedulingMode mode);

bool is_performance_worker(const CpuInfo&info,unsigned worker_index,
						   unsigned thread_count,CoreSchedulingMode mode);

std::uint32_t choose_worker_segment_batch(const CpuInfo&info,
										  unsigned worker_index,
										  unsigned thread_count,
										  std::uint64_t range_span,
										  CoreSchedulingMode mode);

} // namespace calcprime
