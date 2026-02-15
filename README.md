
# calcprimelist
[![License: MIT](https://img.shields.io/badge/License-MIT-2ea44f)](./LICENSE)© 2025 lzray

[![Lines of Code](https://img.shields.io/endpoint?url=https%3A%2F%2Ftokei.kojix2.net%2Fbadge%2Fgithub%2Flzray-universe%2Fcalcprimelist-cpp%2Flines)](https://tokei.kojix2.net/github/lzray-universe/calcprimelist-cpp)

[![Top Language](https://img.shields.io/endpoint?url=https%3A%2F%2Ftokei.kojix2.net%2Fbadge%2Fgithub%2Flzray-universe%2Fcalcprimelist-cpp%2Flanguage)](https://tokei.kojix2.net/github/lzray-universe/calcprimelist-cpp)

English README: [![English](https://img.shields.io/badge/lang-English-informational)](README_ENG.md)

[![Build (Linux + Windows) for Release](https://github.com/lzray-universe/calcprimelist-cpp/actions/workflows/release.yml/badge.svg)](https://github.com/lzray-universe/calcprimelist-cpp/actions/workflows/release_v3.0.yml)

一个可并行的**分段筛素数**工具与库（C++20）。支持：

* 区间计数 `π(B)−π(A)`、区间打印、区间内第 *K* 个素数
* 轮因子（wheel）预筛：`mod 30 / 210 / 1155`
* 自动依据 CPU 缓存与线程数选取分段/分块尺寸
* 三种基础输出（`text` / `binary` / `delta16`）与可选 Zstd 压缩
* Meissel–Lehmer 质数计数
* Miller–Rabin 素性测试
* C++ 静态库、跨语言 C ABI（DLL/.so）调用



---

## 目录

* [快速开始](#快速开始)
* [命令行用法与运行示例](#命令行用法与运行示例)
* [输出格式](#输出格式)
* [库集成（CMake，C++）](#库集成cmakec)
* [DLL / C ABI 调用说明](#dll--c-abi-调用说明)
* [算法与数据结构](#算法与数据结构)
* [说明与提示](#说明与提示)

---

## 快速开始

### 编译

依赖：CMake ≥ 3.16、C++20 编译器（GCC/Clang/MSVC）；**Zstd 开发包是可选依赖**（仅在启用 `--zstd` 时需要）

安装 Zstd（任选其一）：
- Ubuntu / Debian：`sudo apt-get install -y libzstd-dev`
- Fedora / CentOS / RHEL：`sudo dnf install -y libzstd-devel`（或 `yum install -y libzstd-devel`）
- Arch：`sudo pacman -S zstd`
- macOS（Homebrew）：`brew install zstd`
- Windows：
  - vcpkg：`vcpkg install zstd:x64-windows`
    - CMake：`-DCMAKE_TOOLCHAIN_FILE=C:/vcpkg/scripts/buildsystems/vcpkg.cmake`
  - 或手工下载预编译包，并在 CMake 配置时显式指定路径（见下）

若要启用 Zstd 且 CMake 未能自动找到头文件/库，可显式传入：
```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DZSTD_INCLUDE_DIR=/path/to/include \
  -DZSTD_LIBRARY=/path/to/libzstd.{a,so,dylib,lib}
````

> 若未安装 Zstd，项目仍可正常构建（仅禁用 `--zstd` 功能）。

```bash
# Linux / macOS
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j

# Windows（x64，VS 生成器）
cmake -S . -B build -A x64 -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release
```

生成目标：

* 可执行程序：`calcprimelist`（Windows: `calcprimelist.exe`）
* 可执行程序（静态链接项目静态库）：`calcprimelist-static`
* 静态库：`calcprime`（以及同名的 `calcprime_static`）
* 共享库（导出 C ABI，便于 DLL/FFI）：`calcprime-cli`

> 说明：`--help` 输出里的命令前缀仍显示为 `prime-sieve`，在本仓库直接构建后请使用 `calcprimelist`（Windows 下为 `calcprimelist.exe`）。

### 快速试用

```bash
# 计数：2 ≤ n < 100000 的素数个数（应为 9592）
./build/calcprimelist --to 100000 --count --time

# 打印：1 ≤ n < 100 的素数（逐行）
./build/calcprimelist --from 1 --to 100 --print

# 第 K 个素数：区间 [1, 100000) 内的第 100 个素数
./build/calcprimelist --from 1 --to 100000 --nth 100

# 保存到文件（文本）：每行一个素数
./build/calcprimelist --to 1000000 --print --out primes.txt

# 保存到文件（二进制 little-endian uint64）
./build/calcprimelist --to 1000000 --print --out primes.bin --out-format binary

# 保存到文件（delta16 + Zstd；更省空间/带宽）
./build/calcprimelist --to 10000000 --print --out primes.zst --out-format delta16 --zstd
```

> 支持科学计数法、十六进制与大小后缀：`--to 1e8`、`--to 0x1e5`、`--segment 1M`、`--tile 256K` 等。

---

## 命令行用法与运行示例

```
calcprimelist --from A --to B [options]

  主功能（三选一，默认 --count）：
  --count             统计区间 [A, B) 的素数个数
  --print             打印区间素数（输出到 stdout 或文件）
  --nth K             返回区间 [A, B) 内第 K 个素数

  性能/正确性相关：
  --threads N         指定线程数（缺省 0=自动，取决于 CPU）
  --wheel 30|210|1155 轮因子选择（默认 30）
  --segment BYTES     覆盖分段大小（默认依据缓存自适应）
  --tile BYTES        覆盖分块大小（默认依据缓存自适应）

  输出与统计：
  --out PATH          将输出写入文件（默认 stdout）
  --out-format FMT    text（默认）| binary | delta16
  --zstd              对输出字节流做 zstd 压缩（若构建支持）
  --time              打印耗时（微秒）
  --stats             打印配置统计（线程、缓存、分段等）

  其他：
  --ml                用 Meissel-Lehmer 做计数（仅 --count）
  --wheel-bitmap      强制使用 wheel-bitmap 计数路径
  --stest             自动基准测试（1e6..1e11，每点 10 次）
  --test N            对 N 做 Miller-Rabin 素性测试
  --help/-h           打印帮助
```

**示例**：

```bash
# 1) 计数
./calcprimelist --to 1e6 --count --time
# 期望输出（示例）：78498
# Elapsed: 123456 us

# 2) 打印并压缩（delta16 + Zstd）
./calcprimelist --from 1 --to 1e7 --print --out primes.zst --out-format delta16 --zstd

# 3) 在大区间内找第 1e5 个素数（建议单线程以降低内存占用峰值）
./calcprimelist --from 1 --to 1e8 --nth 100000 --threads 1 --time

# 4) 指定更强的轮因子与更大分段（高吞吐场景）
./calcprimelist --to 1e9 --count --wheel 210 --segment 8M --tile 256K --time

# 5) 一键性能基准（自动识别主机信息并输出 Markdown 表）
./calcprimelist --stest

# 6) 单点素性测试（无需 --to）
./calcprimelist --test 1000000007
```

---

## 输出格式

### `text`（默认）

* 一行一个十进制素数，以换行 `\n` 分隔。
* 适合人读或简单的管道处理（例如 `wc -l` 计数）。

### `binary`

* **每个素数一个 `uint64_t`（little-endian）连续写入**，无头部。
* 读取示例（Python）：

  ```python
  import struct
  primes = []
  with open('primes.bin','rb') as f:
      data = f.read()
      for i in range(0, len(data), 8):
          (x,) = struct.unpack('<Q', data[i:i+8])
          primes.append(x)
  ```

### `delta16`

* 开头写第一个素数（`uint64_t` little-endian）。
* 从第二个素数开始，每个值写为与前一个素数的差（`int16_t` little-endian，正数）。
* 若差值超过 `INT16_MAX`，会直接报错（不使用逃逸码/变长编码）。

### `--zstd`（可选压缩开关）

* `--zstd` 不再是输出格式，而是对上述 `text` / `binary` / `delta16` 任一格式的字节流进行 zstd frame 流式压缩。
* 若当前构建不支持 zstd，传入 `--zstd` 会报错：`zstd not supported in this build`。
* 兼容别名：`--out-format zstd` 或 `--out-format zstd+delta` 等价于 `--out-format delta16 --zstd`（已弃用）。

---

## 库集成（CMake，C++）

### 作为子工程引入

在你的 `CMakeLists.txt`：

```cmake
add_subdirectory(calcprimelist)
target_link_libraries(your_target PRIVATE calcprime)
target_include_directories(your_target PRIVATE calcprimelist/include)
```

### 直接使用 C++ 接口（头文件）

常用头文件与函数（命名空间 `calcprime`）：

* `#include <base_sieve.h>`
  `std::vector<uint32_t> simple_sieve(uint64_t limit);`
* `#include <prime_count.h>`
  `uint64_t meissel_count(uint64_t from, uint64_t to, const std::vector<uint32_t>& primes, unsigned threads=0);`
  `bool miller_rabin_is_prime(uint64_t n);`
* `#include <popcnt.h>`
  `uint64_t popcount_u64(uint64_t);`
  `uint64_t count_zero_bits(const uint64_t* bits, size_t bit_count);`

> `meissel_count` 需要提供覆盖 `sqrt(to)` 的素数表，可先用 `simple_sieve` 生成。

示例：

```cpp
#include "base_sieve.h"
#include "prime_count.h"
#include <cmath>
#include <iostream>

int main() {
    auto primes = calcprime::simple_sieve(100);
    std::cout << "π(100) = " << primes.size() << "\n";

    // 计数 [1, 1e8)
    std::uint64_t to = 100000000;
    auto base = calcprime::simple_sieve(
        static_cast<std::uint64_t>(std::sqrt(static_cast<long double>(to))) + 1);
    std::uint64_t cnt = calcprime::meissel_count(1, to, base);
    std::cout << cnt << "\n";
}
```

---

## DLL / C ABI 调用说明

共享库目标名：**`calcprime-cli`**（Windows: `.dll`+`.lib`，Linux: `.so`，macOS: `.dylib`）。所有 C 接口声明位于：`include/calcprime/api.h`。

### 关键类型（摘录）

```c
typedef struct calcprime_cpu_info {
    unsigned    logical_cpus;
    unsigned    physical_cpus;
    size_t      l1_data_bytes;
    size_t      l2_bytes;
    size_t      l2_total_bytes;
    int         has_smt;
} calcprime_cpu_info;

typedef struct calcprime_segment_config {
    size_t      segment_bytes;
    size_t      tile_bytes;
    size_t      segment_bits;
    size_t      tile_bits;
    uint64_t    segment_span;  // 本段覆盖的数值跨度（仅奇数）
    uint64_t    tile_span;
} calcprime_segment_config;

typedef enum calcprime_status {
    CALCPRIME_STATUS_SUCCESS           = 0,
    CALCPRIME_STATUS_INVALID_ARGUMENT  = 1,
    CALCPRIME_STATUS_CANCELLED         = 2,
    CALCPRIME_STATUS_IO_ERROR          = 3,
    CALCPRIME_STATUS_INTERNAL_ERROR    = 4
} calcprime_status;

typedef enum calcprime_wheel_type {
    CALCPRIME_WHEEL_MOD30   = 0,
    CALCPRIME_WHEEL_MOD210  = 1,
    CALCPRIME_WHEEL_MOD1155 = 2
} calcprime_wheel_type;

typedef enum calcprime_output_format {
    CALCPRIME_OUTPUT_TEXT        = 0,
    CALCPRIME_OUTPUT_BINARY      = 1,
    CALCPRIME_OUTPUT_DELTA16     = 2,
    CALCPRIME_OUTPUT_ZSTD_DELTA  = CALCPRIME_OUTPUT_DELTA16 // deprecated alias
} calcprime_output_format;

struct calcprime_cancel_token;
typedef struct calcprime_cancel_token calcprime_cancel_token;

typedef int (*calcprime_prime_chunk_callback)(const uint64_t* primes,
                                              size_t count,
                                              void* user_data);

typedef int (*calcprime_progress_callback)(double progress, void* user_data);
```

**运行选项（核心）**：

```c
typedef struct calcprime_range_options {
    uint64_t    from;
    uint64_t    to;
    unsigned    threads;         // 0 = 自动
    calcprime_wheel_type wheel;  // 30/210/1155
    size_t      segment_bytes;   // 0 = 自动
    size_t      tile_bytes;      // 0 = 自动
    uint64_t    nth_index;       // 0 表示不启用 --nth
    int         collect_primes;  // 1=收集各段素数(或通过回调流式获取)
    int         use_meissel;     // 1=用 Meissel-Lehmer 做计数
    int         write_to_file;   // 1=写文件，配合 output_path/output_format
    calcprime_output_format output_format;
    const char* output_path;     // NULL/空=stdout（当 write_to_file=1 时）
    calcprime_prime_chunk_callback  prime_callback;       // 可选：流式回调
    void*       prime_callback_user_data;
    calcprime_progress_callback     progress_callback;    // 可选：进度回调
    void*       progress_user_data;
    calcprime_cancel_token*        cancel_token;         // 可选：可取消
    int         compress_zstd;      // 1=启用 zstd 压缩（若构建支持）
} calcprime_range_options;
```

**结果与统计（摘录）**：

```c
typedef struct calcprime_range_stats {
    uint64_t    from, to;
    uint64_t    elapsed_us;
    unsigned    threads;
    calcprime_cpu_info     cpu;
    calcprime_segment_config segment;
    calcprime_wheel_type   wheel;
    calcprime_output_format output_format;
    size_t      segments_total, segments_processed;
    uint64_t    prime_count;
    uint64_t    nth_index;
    int         nth_found;
    int         use_meissel;
    int         completed;
    int         cancelled;
} calcprime_range_stats;
```

### 常用函数

```c
// CPU/配置
calcprime_cpu_info        calcprime_detect_cpu_info(void);
unsigned                  calcprime_effective_thread_count(const calcprime_cpu_info*);
calcprime_segment_config  calcprime_choose_segment_config(const calcprime_cpu_info*,
                                                          unsigned threads,
                                                          size_t requested_segment_bytes,
                                                          size_t requested_tile_bytes,
                                                          uint64_t range_length);

// 运行
int  calcprime_run_cli(int argc, char** argv);
int  calcprime_range_options_init(calcprime_range_options* options);
calcprime_cancel_token*  calcprime_cancel_token_create(void);
void calcprime_cancel_token_destroy(calcprime_cancel_token*);
void calcprime_cancel_token_request(calcprime_cancel_token*);

typedef struct calcprime_range_run_result calcprime_range_run_result;
calcprime_status calcprime_run_range(const calcprime_range_options*,
                                     calcprime_range_run_result** out);

// 查询结果
calcprime_status  calcprime_range_result_status(const calcprime_range_run_result*);
const char*       calcprime_range_result_error_message(const calcprime_range_run_result*);
uint64_t          calcprime_range_result_count(const calcprime_range_run_result*);
int               calcprime_range_result_nth_prime(const calcprime_range_run_result*, uint64_t* out_value);
int               calcprime_range_result_stats(const calcprime_range_run_result*, calcprime_range_stats* out);

// 按段访问与汇总拷贝
size_t            calcprime_range_result_segment_count(const calcprime_range_run_result*);
// A) 直接返回内部指针（只读，生命周期受 result 管理）
int               calcprime_range_result_segment(const calcprime_range_run_result*, size_t index,
                                                 const uint64_t** out_primes, size_t* out_count);
// B) 拷贝所有已收集素数到用户缓冲区（buffer=NULL 时可仅探测 out_written）
int               calcprime_range_result_copy_primes(const calcprime_range_run_result*,
                                                     uint64_t* buffer, size_t capacity, size_t* out_written);

void              calcprime_range_result_release(calcprime_range_run_result*);
```

### 典型用法（C）

```c
#include "calcprime/api.h"
#include <stdio.h>

static int on_chunk(const uint64_t* primes, size_t n, void* ud) {
    FILE* fp = (FILE*)ud;
    for (size_t i = 0; i < n; ++i) fprintf(fp, "%llu\n", (unsigned long long)primes[i]);
    return 0; // 返回非零可中断
}

int main() {
    calcprime_range_options opt;
    calcprime_range_options_init(&opt);
    opt.from = 1; opt.to = 1000000;
    opt.threads = 0; // 自动
    opt.wheel = CALCPRIME_WHEEL_MOD210;
    opt.output_format = CALCPRIME_OUTPUT_TEXT;
    opt.collect_primes = 0;      // 不保存在 result 中
    opt.prime_callback = on_chunk;  // 改为流式回调输出
    opt.prime_callback_user_data = stdout;

    calcprime_range_run_result* res = NULL;
    calcprime_status st = calcprime_run_range(&opt, &res);
    if (st != CALCPRIME_STATUS_SUCCESS) {
        fprintf(stderr, "error: %s\n", calcprime_range_result_error_message(res));
    } else {
        calcprime_range_stats stats;
        calcprime_range_result_stats(res, &stats);
        fprintf(stderr, "count=%llu, elapsed=%llu us\n",
                (unsigned long long)stats.prime_count,
                (unsigned long long)stats.elapsed_us);
    }
    if (res) calcprime_range_result_release(res);
}
```

> 说明：也提供 `calcprime_simple_sieve/…_release_u32_buffer` 与 `calcprime_meissel_count` / `calcprime_miller_rabin_is_prime` 等函数，可独立调用。

---

## 算法与数据结构

本仓库实现的是**位图分段筛**（仅表示奇数），配合**轮因子预筛**、**分块（tile）**、**多级素因子分类**与**桶环（Bucket Ring）**调度大素因子命中。整体结构与数据流如下：

### 1. 数轴映射与位图表示

* 仅存储奇数：对区间 `[L, R)`（内部会调整到奇数边界），构造位图 `bitset` ，**1 表示合数，0 表示素数**。
* 映射关系：`index i  <->  value = seg_low + 2*i`。
* 每段（segment）长度为 `segment_span`，以字节对齐（源码中常见对齐到 128B）；段内再划分为小块（tile），更好地命中 **L1D**。

相关代码：`segmenter.*` / `marker.*` / `popcnt.*`

### 2. 轮因子（Wheel）预筛

* 轮类型：`mod 30 / 210 / 1155`（默认 30）。
* 预先构建 `Wheel`：

  * `allowed[ residue ]`：该余数是否可能为素（与 modulus 互素）。
  * `residues`：所有允许的奇数余数集合（加速遍历）。
  * `steps`：从一个允许余数到下一个允许余数需要前进的 2 的倍数（用于快速跳过不可能的位置）。
  * `small_prime_patterns`：对**小素因子**的周期性掩码（见下一节）。

`Wheel::apply_presieve(start, bit_count, bits)` 会**一次性将所有不允许的余数标记为合数**（位图置 1），极大减少后续标记工作。

相关代码：`wheel.*`

### 3. 小/中/大素因子三段式标记

在分段筛中，按素因子大小划分处理方式：

* **Small primes（小素因子）**
  对于较小的 `p`，构建**相位模式（phase）**与**掩码表**，实现**批量按位标记**。这避免了逐个加 `p` 的分支与访存。结构体（简化）：

  ```cpp
  struct SmallPrimePattern {
      uint32_t prime;
      uint32_t phase_count;           // = p
      std::vector<uint64_t> masks;    // 每个相位对应的一组 64-bit 掩码
      std::vector<uint32_t> next_phase;
  };
  ```

  `apply_small_primes` 会根据当前 tile 的起始相位选择合适的掩码，直接对位图做 `OR`。

* **Medium primes（中素因子）**
  仍然在一个段/块内能“较密集”命中，但不再适合完全模板化。做法是：

  1. 计算在本段的**首命中** `first_hit(p, start)`；
  2. 采用**允许余数步进**（来自 wheel 的 `steps`）进行快速跳跃标记；
  3. 在 tile 边界对齐时复用偏移以减少整除与取模。

* **Large primes（大素因子）**
  大素因子在当前段内命中很稀疏，且可能跨多个段。使用**桶环（BucketRing）**将“下一次命中”投递到未来的某个段，并在该段到来时批量处理：

  ```cpp
  struct LargePrimeState {
      uint32_t    prime;
      uint64_t    next_value;   // 下一次命中值
      uint64_t    stride;       // 前进步长（结合 wheel）
  };

  struct BucketEntry {
      uint32_t    prime;
      uint64_t    next_index;   // 下一命中在该段内的 index
      uint64_t    offset;       // 相对段起点的偏移
      uint64_t    value;        // 下一命中的绝对值（调试/断言用）
      LargePrimeState* owner;   // 回写 next_value / next_index
  };

  class BucketRing {
      uint64_t base_segment;            // 环的起始段号
      size_t   mask;                    // 环大小 - 1，按 2^k 管理
      std::vector<std::vector<BucketEntry>> buckets;
      // push(segment, entry), take(segment) ...
  };
  ```

  这样每个段只处理**正好命中到该段**的那些大素因子，大幅减少跨段扫描开销。

相关代码：`marker.*` / `bucket.*` / `wheel.*`

### 4. 分段/分块与任务调度

* **SegmentWorkQueue**：全局原子段号 `next_segment_`，工作线程调用 `next(...)` 领取下一个待处理段。
* **尺寸**：`choose_segment_config(cpu, requested_segment, requested_tile, range_length)` 综合 L1D/L2/线程数等信息给出 `segment_bytes/tile_bytes/…`；也可用命令行覆盖。
* **多线程**：每个线程独立持有临时位图与本地桶结构，避免共享写冲突，仅在**结果**与**进度**上用条件变量/原子做同步。

相关代码：`segmenter.*` / `cpu_info.*`

### 5. 计数与输出

* **计数**：位图就绪后调用 `count_zero_bits(bits, bit_count)`，配合 AVX2/AVX-512（如可用）的 `popcnt` 变体优化。
* **输出**：`PrimeWriter` 维护一个 I/O 线程与**块队列**（`Chunk`），前端将 `text`/`binary`/`delta16` 编码后的块入队；后端顺序写文件/stdout，且仅在 writer 线程中按需执行 zstd 流式压缩。

相关代码：`popcnt.*` / `writer.*`

### 6. Meissel–Lehmer 计数（`--ml`）

* 目标：计算 `π(N)`（或区间计数拆解成端点差）。
* 经典公式：
  `π(n) = φ(n, a) + a - 1 - P2(n, a)`，其中常取 `a = π(n^{1/4})`；
  `P2` 项可由 `b = π(√n)`、`c = π(n^{1/3})` 等辅助分界高效计算。
* 实现要点：

  * `phi(n, a)` 采用**递归/记忆化**与**分块**；
  * `π(x)` 小范围使用普通筛，大范围调用上述分段筛；
  * 线程并行拆分子任务（未来可能受限于递归依赖）。

相关代码：`prime_count.*`

### 7. Miller–Rabin 素性测试

* 针对 64 位整数的确定版底集合实现（典型多个固定 bases），对单点 `N` 的快速判定。
* 命令行：`--test N`，输出 `prime` 或 `composite`。

相关代码：`prime_count.*`

---
## 性能基准（π(x) 计数，多线程）

> 对每个区间上界 `--to ∈ {1e6, 1e7, 1e8, 1e9, 1e10, 1e11}`，运行
> `calcprimelist --to {N} --count --time` 重复 10 次，统计 **best / median / mean**。
> Windows 等价命令：`.\calcprimelist.exe --to {N} --count --time`。
> 注：`--time` 为工具内计时（微秒），若不可用则以外部 wall-clock 计时。具体环境（CPU/内存/OS/编译器）会影响绝对值。
> 也可直接运行：`calcprimelist --stest` 自动输出同格式表格。

**测试数据（10 次/每点）：**（本表来自 `Intel® Core™ Ultra 7 265K @ P 3.90 / E 3.30 GHz · 20C/20T (8P+12E)` `Microsoft Windows 11 Pro for Workstations, Version 25H2 (OS Build 26220.7052), x64`）

|     N (`--to`) | runs | ok | best (s) | median (s) | mean (s) | throughput median (range/s) |
| -------------: | ---: | -: | -------: | ---------: | -------: | --------------------------: |
|    1,000,000 |   10 | 10 |  0.001010 |    0.001169 |  0.001196 |                855,066,267 |
|   10,000,000 |   10 | 10 |  0.001332 |    0.001587 |  0.001622 |              6,303,183,107 |
|  100,000,000 |   10 | 10 |  0.003483 |    0.003680 |  0.003706 |             27,173,913,043 |
| 1,000,000,000 |   10 | 10 |  0.019459 |    0.019804 |  0.020066 |             50,494,849,525 |
| 10,000,000,000 |   10 | 10 |  0.190158 |    0.192444 |  0.193137 |             51,963,168,506 |
| 100,000,000,000 |   10 | 10 |  1.901580 |    1.924440 |  1.931370 |             51,963,168,506 |

---

## 说明与提示

* **线程数**：`--threads 0`（默认）将依据物理/逻辑核与 SMT 自动选择；也可用 `calcprime_detect_cpu_info` / `calcprime_effective_thread_count` 在应用层拿到推荐值。
* **轮因子**：`--wheel 210` 在大区间往往更快；`1155` 预筛最强，但掩码/步进表更大，小区间未必划算。
* **分段/分块**：若清楚目标平台缓存，可手动设定 `--segment / --tile`；一般保证 **tile ≤ L1D，segment 近似 L2** 会有较好效果。
* **寻找第 K 个素数**：若内存紧/更稳定，可用 `--threads 1`；并行情况下内部会以段计数推进，也能找到，但需要额外同步与（可能）二次扫描某些段。
* **输出吞吐**：批量写文件时，优先 `--out-format binary`，或 `--out-format delta16 --zstd`。文本输出人类友好但对磁盘/带宽不友好。
* **边界**：所有计算在 `uint64_t` 范围内进行；请确保 `--from/--to` 满足 `0 ≤ from < to` 且上界不溢出。内部仅标记奇数，`2` 会在前缀处理中单独考虑。
* **测试**：`ctest` 中含有示例（如 `--to 100000 --count --time`）。
