#pragma once

#include <atomic>
#include <condition_variable>
#include <cstdint>
#include <cstdio>
#include <deque>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

namespace calcprime {

enum class PrimeOutputFormat {
    Text,
    Binary,
    ZstdDelta,
};

class PrimeWriter {
public:
PrimeWriter(bool enabled,
            const std::string &path="",
            PrimeOutputFormat format=PrimeOutputFormat::Text);
~PrimeWriter();

bool enabled() const {
    return enabled_;
}
void write_segment(const std::vector<std::uint64_t> &primes);
void write_value(std::uint64_t value);
void flush();
void finish();

private:
struct Chunk {
    std::string data;
    bool flush=false;
};

void enqueue_chunk(Chunk &&chunk);
void writer_loop();
void flush_buffer();
void check_io_error() const;
void set_error(const std::string &message);
std::string encode_deltas(const std::vector<std::uint64_t> &primes);
std::string encode_delta_value(std::uint64_t value);

bool enabled_;
std::FILE *file_;
bool owns_file_;
std::thread writer_thread_;

std::mutex queue_mutex_;
std::condition_variable queue_not_empty_;
std::condition_variable queue_not_full_;
std::deque<Chunk> queue_;
std::size_t queue_capacity_;
bool stop_requested_;

std::string buffer_;
std::size_t buffer_threshold_;

PrimeOutputFormat format_;
std::uint64_t previous_prime_;

mutable std::mutex error_mutex_;
std::atomic<bool> io_error_;
std::string error_message_;
};

} // namespace calcprime

