#pragma once

#include <cstddef>
#include <cstdint>
#include <string>

namespace calcprime {

struct CpuInfo {
    unsigned logical_cpus=1;
    unsigned physical_cpus=1;
    std::size_t l1_data_bytes=32*1024;
    std::size_t l2_bytes=1024*1024;
    std::size_t l2_total_bytes=0;
    bool has_smt=false;
};

CpuInfo detect_cpu_info();

unsigned effective_thread_count(const CpuInfo&info);

}
