#pragma once

#include<atomic>
#include<condition_variable>
#include<cstdint>
#include<cstdio>
#include<deque>
#include<mutex>
#include<string>
#include<thread>
#include<vector>

namespace calcprime{

enum class PrimeOutputFormat{
	Text,
	Binary,
	Delta16,
};

class PrimeWriter{
  public:
	PrimeWriter(bool enabled,const std::string&path="",
				PrimeOutputFormat format=PrimeOutputFormat::Text,
				bool use_zstd=false);
	~PrimeWriter();

	bool enabled() const{ return enabled_; }
	void write_segment(const std::vector<std::uint64_t>&primes);
	void write_value(std::uint64_t value);
	void flush();
	void finish();

  private:
	struct Chunk{
		std::string data;
		bool flush=false;
	};

	void enqueue_chunk(Chunk&&chunk);
	void writer_loop();
	void flush_buffer();
	void check_io_error() const;
	void set_error(const std::string&message);
	std::string encode_delta16(const std::vector<std::uint64_t>&primes);
	std::string encode_delta16_value(std::uint64_t value);
	void write_file_bytes(const char*data,std::size_t size);
#if defined(CALCPRIME_HAS_ZSTD)
	void flush_zstd_stream(bool final_frame);
#endif

	bool enabled_;
	std::FILE*file_;
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
	bool use_zstd_;
	bool has_first_prime_;
	std::uint64_t previous_prime_;
	void*zstd_cctx_;
	std::string zstd_out_buffer_;

	mutable std::mutex error_mutex_;
	std::atomic<bool> io_error_;
	std::string error_message_;
};

} // namespace calcprime
