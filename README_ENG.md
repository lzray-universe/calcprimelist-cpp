# calcprimelist

[![License: MIT](https://img.shields.io/badge/License-MIT-2ea44f)](./LICENSE)© 2025 lzray

A parallelizable **segmented prime sieve** tool & library (C++20). It supports:

* Interval counting `π(B) − π(A)`, printing primes in a range, and the *K*-th prime within a range
* Wheel pre-sieving: `mod 30 / 210 / 1155`
* Auto-tuned segment/tile sizes based on CPU cache & thread count
* Binary output and optional Zstd + Δ encoding
* Meissel–Lehmer prime counting
* Miller–Rabin primality testing
* C++ static library and a cross-language C ABI (DLL/.so)

---

## Table of Contents

* [Quick Start](#quick-start)
* [CLI Usage & Examples](#cli-usage--examples)
* [Output Formats](#output-formats)
* [Library Integration (CMake, C++)](#library-integration-cmake-c)
* [DLL / C ABI](#dll--c-abi)
* [Algorithms & Data Structures](#algorithms--data-structures)
* [Notes & Tips](#notes--tips)

---

## Quick Start

### Build

Dependencies: CMake ≥ 3.16, a C++20 compiler (GCC/Clang/MSVC), and **Zstd development files (≥ 1.5)**

Install Zstd:
- Ubuntu / Debian: `sudo apt-get install -y libzstd-dev`
- Fedora / CentOS / RHEL: `sudo dnf install -y libzstd-devel`
- Arch: `sudo pacman -S zstd`
- macOS (Homebrew): `brew install zstd`
- Windows:
  - vcpkg: `vcpkg install zstd:x64-windows`
    - CMake: `-DCMAKE_TOOLCHAIN_FILE=C:/vcpkg/scripts/buildsystems/vcpkg.cmake`
  - or use a prebuilt package and pass paths explicitly:

If auto-detection fails, pass the paths explicitly:
```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DZSTD_INCLUDE_DIR=/path/to/include \
  -DZSTD_LIBRARY=/path/to/libzstd.{a,so,dylib,lib}
````

> Common error: `Failed to locate zstd headers or library` — install Zstd dev files or provide the paths above.

```bash
# Linux / macOS
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j

# Windows (x64, VS generator)
cmake -S . -B build -A x64 -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release
```

Build targets:

* Executable: `prime-sieve`
* Executable (linked against the static project library): `prime-sieve-static`
* Static library: `calcprime` (and a same-named `calcprime_static`)
* Shared library (exports a C ABI for DLL/FFI): `calcprime-cli`

### Quick Try

```bash
# Count: number of primes with 2 ≤ n < 100000 (should be 9592)
./build/prime-sieve --to 100000 --count --time

# Print: primes with 1 ≤ n < 100 (one per line)
./build/prime-sieve --from 1 --to 100 --print

# K-th prime: the 100th prime within [1, 100000)
./build/prime-sieve --from 1 --to 100000 --nth 100

# Save to file (text): one prime per line
./build/prime-sieve --to 1000000 --print --out primes.txt

# Save to file (binary little-endian uint64)
./build/prime-sieve --to 1000000 --print --out primes.bin --out-format binary

# Save to file (Zstd + Δ; space/bandwidth friendly)
./build/prime-sieve --to 10000000 --print --out primes.zst --out-format zstd
```

> Scientific notation and size suffixes are supported: `--to 1e8`, `--segment 1M`, `--tile 256K`, etc.

---

## CLI Usage & Examples

```
prime-sieve --from A --to B [options]

  Main modes (choose one; default is --count):
  --count             Count primes in [A, B)
  --print             Print primes in the interval (to stdout or file)
  --nth K             Return the K-th prime within [A, B)

  Performance / correctness:
  --threads N         Number of threads (default 0 = auto, based on CPU)
  --wheel 30|210|1155 Wheel selection (default 30)
  --segment BYTES     Override segment size (default: cache-aware)
  --tile BYTES        Override tile size (default: cache-aware)

  Output & stats:
  --out PATH          Write output to file (default stdout)
  --out-format FMT    text (default) | binary | zstd
  --time              Print elapsed time (microseconds)
  --stats             Print configuration stats (threads, cache, segments, etc.)

  Misc:
  --ml                Use Meissel–Lehmer for counting (only with --count)
  --test N            Miller–Rabin primality test for N
  --help/-h           Show help
```

**Examples:**

```bash
# 1) Counting
./prime-sieve --to 1e6 --count --time
# Expected (example): 78498
# Elapsed: 123456 us

# 2) Print and compress (Zstd + Δ)
./prime-sieve --from 1 --to 1e7 --print --out primes.zst --out-format zstd

# 3) Find the 1e5-th prime in a large interval (single-thread to reduce peak memory)
./prime-sieve --from 1 --to 1e8 --nth 100000 --threads 1 --time

# 4) Stronger wheel and larger segment sizes (throughput-oriented)
./prime-sieve --to 1e9 --count --wheel 210 --segment 8M --tile 256K --time
```

---

## Output Formats

### `text` (default)

* One decimal prime per line, separated by `\n`.
* Human-friendly, easy to pipe (e.g., `wc -l` for counting).

### `binary`

* **Each prime is written as a `uint64_t` (little-endian) consecutively**, no header.
* Reading example (Python):

  ```python
  import struct
  primes = []
  with open('primes.bin','rb') as f:
      data = f.read()
      for i in range(0, len(data), 8):
          (x,) = struct.unpack('<Q', data[i:i+8])
          primes.append(x)
  ```

### `zstd` (Zstd + Δ)

* Primes in ascending order are **delta-encoded (Δ)**, then compressed with Zstd (implemented as `ZstdDelta`).
* Decoding: Zstd-decompress to obtain a little-endian `uint64_t` Δ sequence; then prefix-sum to recover values.
* Suited for **large-scale** export (bandwidth/storage friendly).

> Note: the `writer` module encapsulates Δ encoding; whether Zstd is enabled depends on build settings/environment. For long-term archiving, prefer `--out-format zstd`.

---

## Library Integration (CMake, C++)

### Add as a subproject

In your `CMakeLists.txt`:

```cmake
add_subdirectory(calcprimelist)
target_link_libraries(your_target PRIVATE calcprime)
target_include_directories(your_target PRIVATE calcprimelist/include)
```

### Use the C++ API directly (headers)

Common headers & functions (namespace `calcprime`):

* `#include <base_sieve.h>`
  `std::vector<uint32_t> simple_sieve(uint64_t limit);`
* `#include <prime_count.h>`
  `uint64_t meissel_count(uint64_t from, uint64_t to, unsigned threads=0);`
  `bool miller_rabin_is_prime(uint64_t n);`
* `#include <popcnt.h>`
  `uint64_t popcount_u64(uint64_t);`
  `uint64_t count_zero_bits(const uint64_t* bits, size_t bit_count);`

Example:

```cpp
#include "base_sieve.h"
#include "prime_count.h"
#include <iostream>

int main() {
    auto primes = calcprime::simple_sieve(100);
    std::cout << "π(100) = " << primes.size() << "\n";

    // Count [1, 1e8)
    std::uint64_t cnt = calcprime::meissel_count(1, 100000000);
    std::cout << cnt << "\n";
}
```

---

## DLL / C ABI

Shared library target name: **`calcprime-cli`** (Windows: `.dll` + `.lib`, Linux: `.so`, macOS: `.dylib`). All C interfaces are declared in: `include/calcprime/api.h`.

### Key types (excerpt)

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
    uint64_t    segment_span;  // numeric span covered by this segment (odd-only)
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
    CALCPRIME_OUTPUT_ZSTD_DELTA  = 2
} calcprime_output_format;

struct calcprime_cancel_token;
typedef struct calcprime_cancel_token calcprime_cancel_token;

typedef int (*calcprime_prime_chunk_callback)(const uint64_t* primes,
                                              size_t count,
                                              void* user_data);

typedef int (*calcprime_progress_callback)(double progress, void* user_data);
```

**Runtime options (core):**

```c
typedef struct calcprime_range_options {
    uint64_t    from;
    uint64_t    to;
    unsigned    threads;         // 0 = auto
    calcprime_wheel_type wheel;  // 30/210/1155
    size_t      segment_bytes;   // 0 = auto
    size_t      tile_bytes;      // 0 = auto
    uint64_t    nth_index;       // 0 = disable --nth
    int         collect_primes;  // 1 = collect per-segment primes (or stream via callback)
    int         use_meissel;     // 1 = use Meissel–Lehmer for counting
    int         write_to_file;   // 1 = write to file with output_path/output_format
    calcprime_output_format output_format;
    const char* output_path;     // NULL/empty = stdout (when write_to_file=1)
    calcprime_prime_chunk_callback  prime_callback;       // optional: streaming callback
    void*       prime_callback_user_data;
    calcprime_progress_callback     progress_callback;    // optional: progress callback
    void*       progress_user_data;
    calcprime_cancel_token*        cancel_token;         // optional: cancellable
} calcprime_range_options;
```

**Results & stats (excerpt):**

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

### Common functions

```c
// CPU / configuration
calcprime_cpu_info        calcprime_detect_cpu_info(void);
unsigned                  calcprime_effective_thread_count(const calcprime_cpu_info*);
calcprime_segment_config  calcprime_choose_segment_config(const calcprime_cpu_info*,
                                                          size_t requested_segment_bytes,
                                                          size_t requested_tile_bytes,
                                                          uint64_t range_length);

// Run
int  calcprime_range_options_init(calcprime_range_options* options);
calcprime_cancel_token*  calcprime_cancel_token_create(void);
void calcprime_cancel_token_destroy(calcprime_cancel_token*);
void calcprime_cancel_token_request(calcprime_cancel_token*);

typedef struct calcprime_range_run_result calcprime_range_run_result;
calcprime_status calcprime_run_range(const calcprime_range_options*,
                                     calcprime_range_run_result** out);

// Query results
calcprime_status  calcprime_range_result_status(const calcprime_range_run_result*);
const char*       calcprime_range_result_error_message(const calcprime_range_run_result*);
uint64_t          calcprime_range_result_count(const calcprime_range_run_result*);
int               calcprime_range_result_nth_prime(const calcprime_range_run_result*, uint64_t* out_value);
int               calcprime_range_result_stats(const calcprime_range_run_result*, calcprime_range_stats* out);

// Per-segment access (two ways)
size_t            calcprime_range_result_segment_count(const calcprime_range_run_result*);
// A) Return an internal pointer (read-only; lifetime owned by result)
int               calcprime_range_result_segment(const calcprime_range_run_result*, size_t index,
                                                 const uint64_t** out_primes, size_t* out_count);
// B) Copy into a user buffer
int               calcprime_range_result_copy_primes(const calcprime_range_run_result*, size_t index,
                                                     uint64_t* buffer, size_t capacity, size_t* out_written);

void              calcprime_range_result_release(calcprime_range_run_result*);
```

### Typical usage (C)

```c
#include "calcprime/api.h"
#include <stdio.h>

static int on_chunk(const uint64_t* primes, size_t n, void* ud) {
    FILE* fp = (FILE*)ud;
    for (size_t i = 0; i < n; ++i) fprintf(fp, "%llu\n", (unsigned long long)primes[i]);
    return 0; // return non-zero to interrupt
}

int main() {
    calcprime_range_options opt;
    calcprime_range_options_init(&opt);
    opt.from = 1; opt.to = 1000000;
    opt.threads = 0; // auto
    opt.wheel = CALCPRIME_WHEEL_MOD210;
    opt.output_format = CALCPRIME_OUTPUT_TEXT;
    opt.collect_primes = 0;      // do not store in the result
    opt.prime_callback = on_chunk;  // stream via callback
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

> Note: standalone helpers like `calcprime_simple_sieve/…_release_u32_buffer`, `calcprime_meissel_count`, and `calcprime_miller_rabin_is_prime` are also provided.

---

## Algorithms & Data Structures

This repo implements an **odd-only bitset segmented sieve**, with **wheel pre-sieving**, **tiling**, **multi-tier sieving primes**, and a **Bucket Ring** scheduler for large primes.

### 1. Number-line mapping & bitset

* Store only odd numbers: for range `[L, R)` (internally snapped to odd bounds), build a `bitset` where **1 means composite, 0 means prime**.
* Mapping: `index i  <->  value = seg_low + 2*i`.
* Each segment has length `segment_span`, byte-aligned (often aligned to 128B); the segment is split into small **tiles** for better **L1D** locality.

Relevant code: `segmenter.*` / `marker.*` / `popcnt.*`

### 2. Wheel pre-sieving

* Wheels: `mod 30 / 210 / 1155` (default 30).
* Precompute `Wheel`:

  * `allowed[residue]`: whether a residue can be prime (coprime to the modulus).
  * `residues`: all allowed odd residues (for fast iteration).
  * `steps`: delta (in multiples of 2) to jump from one allowed residue to the next.
  * `small_prime_patterns`: periodic masks for **small sieving primes** (see next section).

`Wheel::apply_presieve(start, bit_count, bits)` **marks all disallowed residues composite** in one go, drastically reducing later work.

Relevant code: `wheel.*`

### 3. Three-tier marking: small / medium / large sieving primes

* **Small primes**
  For small `p`, build **phase patterns** and **mask tables** to perform **batch bit-OR marking**, avoiding per-hit branching and memory traffic:

  ```cpp
  struct SmallPrimePattern {
      uint32_t prime;
      uint32_t phase_count;           // = p
      std::vector<uint64_t> masks;    // masks per phase (64-bit chunks)
      std::vector<uint32_t> next_phase;
  };
  ```

  `apply_small_primes` picks masks based on the tile’s starting phase and ORs them into the bitset.

* **Medium primes**
  Still dense enough within a segment/tile, but no longer worth full templating:

  1. Compute the **first hit** in the current segment: `first_hit(p, start)`;
  2. Use **allowed-residue stepping** (`steps` from the wheel) to jump quickly;
  3. Reuse offsets at tile boundaries to reduce divisions/mods.

* **Large primes**
  Hits are sparse and may cross segments. Use a **BucketRing** to enqueue the **next hit** into a future segment; when that segment arrives, process all scheduled entries:

  ```cpp
  struct LargePrimeState {
      uint32_t    prime;
      uint64_t    next_value;   // next absolute hit
      uint64_t    stride;       // step (wheel-aware)
  };

  struct BucketEntry {
      uint32_t    prime;
      uint64_t    next_index;   // index of next hit within the segment
      uint64_t    offset;       // offset relative to segment base
      uint64_t    value;        // absolute value of the next hit (debug/assert)
      LargePrimeState* owner;   // to update next_value / next_index
  };

  class BucketRing {
      uint64_t base_segment;
      size_t   mask;                    // ring size - 1 (power of two)
      std::vector<std::vector<BucketEntry>> buckets;
      // push(segment, entry), take(segment) ...
  };
  ```

  Each segment handles only the large primes **that actually hit this segment**, cutting cross-segment scanning.

Relevant code: `marker.*` / `bucket.*` / `wheel.*`

### 4. Segmentation/tiling & task scheduling

* **SegmentWorkQueue**: a global atomic segment counter `next_segment_`; worker threads call `next(...)` to fetch work.
* **Sizing**: `choose_segment_config(cpu, requested_segment, requested_tile, range_length)` uses L1D/L2/thread info to choose `segment_bytes/tile_bytes/...`; CLI can override.
* **Multithreading**: each thread owns its local bitset and bucket structures to avoid shared writes; only **results** and **progress** use condition vars/atomics.

Relevant code: `segmenter.*` / `cpu_info.*`

### 5. Counting & output

* **Counting**: after the bitset is ready, call `count_zero_bits(bits, bit_count)`, with AVX2/AVX-512 `popcnt` variants when available.
* **Output**: `PrimeWriter` maintains an I/O thread and a **chunk queue**; front-end enqueues encoded text/binary chunks; back-end writes file/stdout sequentially to reduce I/O impact. In `ZstdDelta` mode, Δ-encode first, then compress/assemble.

Relevant code: `popcnt.*` / `writer.*`

### 6. Meissel–Lehmer counting (`--ml`)

* Goal: compute `π(N)` (or use endpoint differences for interval counts).
* Classic formula:
  `π(n) = φ(n, a) + a − 1 − P2(n, a)`, commonly with `a = π(n^{1/4})`;
  `P2` uses helpers `b = π(√n)`, `c = π(n^{1/3})`, etc.
* Implementation notes:

  * `phi(n, a)` via **recursion/memoization** and **blocking**;
  * small ranges use a plain sieve; large ranges call the segmented sieve above;
  * parallelize subproblems where possible (recursion dependencies may limit this).

Relevant code: `prime_count.*`

### 7. Miller–Rabin primality test

* Deterministic base set for 64-bit integers (fixed bases), fast single-point testing.
* CLI: `--test N`, outputs `prime` or `composite`.

Relevant code: `prime_count.*`

---

## Benchmarks (π(x) counting, multithreaded)

> For each upper bound `--to ∈ {1e6, 1e7, 1e8, 1e9, 1e10}`, run
> `prime-sieve --to {N} --count --time` 10 times and report **best / median / mean**.
> Windows equivalent: `.\prime-sieve.exe --to {N} --count --time`.
> Note: `--time` uses the tool’s internal timing (µs); if unavailable, wall-clock is used. Absolute values depend on CPU/RAM/OS/compiler.

**Test data (10 runs per point):** (table from `Intel® Core™ Ultra 7 265K @ P 3.90 / E 3.30 GHz · 20C/20T (8P+12E)` `Microsoft Windows 11 Pro for Workstations, Version 25H2 (OS Build 26220.7052), x64`)

|     N (`--to`) | runs | ok | best (s) | median (s) | mean (s) | throughput median (range/s) |
| -------------: | ---: | -: | -------: | ---------: | -------: | --------------------------: |
|    1,000,000 |   10 | 10 |  0.001010 |    0.001169 |  0.001196 |                855,066,267 |
|   10,000,000 |   10 | 10 |  0.001332 |    0.001587 |  0.001622 |              6,303,183,107 |
|  100,000,000 |   10 | 10 |  0.003483 |    0.003680 |  0.003706 |             27,173,913,043 |
| 1,000,000,000 |   10 | 10 |  0.019459 |    0.019804 |  0.020066 |             50,494,849,525 |
| 10,000,000,000 |   10 | 10 |  0.190158 |    0.192444 |  0.193137 |             51,963,168,506 |

---

## Notes & Tips

* **Threads**: `--threads 0` (default) auto-selects based on physical/logical cores and SMT; you can also use `calcprime_detect_cpu_info` / `calcprime_effective_thread_count` to obtain recommended values in your app.
* **Wheel**: `--wheel 210` is often faster for large ranges; `1155` pre-sieves the most but has bigger masks/step tables and may not pay off for small ranges.
* **Segments/tiles**: if you know the target cache hierarchy, set `--segment / --tile` manually; as a rule of thumb, **tile ≤ L1D, segment ≈ L2** performs well.
* **Finding the K-th prime**: if memory is tight or you want predictable peaks, consider `--threads 1`. In parallel mode, the tool advances by segment counts and can still find it, with extra synchronization and potential re-scans for some segments.
* **Output throughput**: for bulk export, prefer `--out-format binary` or `--out-format zstd`. Text is human-friendly but not storage/bandwidth-friendly.
* **Bounds**: all computations use `uint64_t`. Ensure `0 ≤ from < to` and the upper bound doesn’t overflow. Only odd numbers are marked; `2` is handled separately in a prefix step.
* **Tests**: `ctest` includes examples (e.g., `--to 100000 --count --time`).
