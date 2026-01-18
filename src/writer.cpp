#include "writer.h"

#include<cerrno>
#include<charconv>
#include<cstring>
#include<exception>
#include<stdexcept>

#if defined(_MSC_VER)
#include<intrin.h>
#endif

namespace calcprime{

namespace{

constexpr std::size_t kDefaultFileBuffer=8u<<20; // 8 MiB
constexpr std::size_t kDefaultQueueCapacity=8;
constexpr std::size_t kDefaultBufferThreshold=8u<<20; // 8 MiB

inline std::uint64_t to_little_endian(std::uint64_t value){
#if defined(__BYTE_ORDER__)&&(__BYTE_ORDER__==__ORDER_BIG_ENDIAN__)
	return __builtin_bswap64(value);
#elif defined(_MSC_VER)&&defined(_WIN32)
	return _byteswap_uint64(value);
#else
	return value;
#endif
}

} // namespace

PrimeWriter::PrimeWriter(bool enabled,const std::string&path,
						 PrimeOutputFormat format)
	: enabled_(enabled),file_(nullptr),owns_file_(false),
	  queue_capacity_(kDefaultQueueCapacity),stop_requested_(false),
	  buffer_threshold_(kDefaultBufferThreshold),format_(format),
	  previous_prime_(0),io_error_(false){
	if(!enabled_){
		return;
	}

	if(path.empty()){
		file_=stdout;
		owns_file_=false;
		std::fprintf(stderr,"[calcprime] warning: writing primes to stdout may "
							"stall large outputs."
							" Consider using --out <path>.\n");
	}else{
		file_=std::fopen(path.c_str(),"wb");
		if(!file_){
			throw std::runtime_error("Failed to open output file");
		}
		owns_file_=true;
	}

	if(!file_){
		throw std::runtime_error("Invalid output handle");
	}

	if(std::setvbuf(file_,nullptr,_IOFBF,kDefaultFileBuffer)!=0){
		throw std::runtime_error("Failed to set file buffer");
	}

	buffer_.reserve(buffer_threshold_);
	queue_.clear();

	writer_thread_=std::thread(&PrimeWriter::writer_loop,this);
}

PrimeWriter::~PrimeWriter(){
	try{
		finish();
	}catch(...){
		std::terminate();
	}
}

void PrimeWriter::write_segment(const std::vector<std::uint64_t>&primes){
	if(!enabled_){
		return;
	}
	if(primes.empty()){
		return;
	}

	switch(format_){
	case PrimeOutputFormat::Text:{
		std::string chunk;
		chunk.reserve(primes.size()*24);
		char local[32];
		for(std::uint64_t value : primes){
			auto result=std::to_chars(local,local+sizeof(local),value);
			if(result.ec!=std::errc()){
				throw std::runtime_error("Failed to convert prime to string");
			}
			chunk.append(local,result.ptr);
			chunk.push_back('\n');
		}
		enqueue_chunk(Chunk{std::move(chunk),false});
		break;
	}
	case PrimeOutputFormat::Binary:{
		std::string chunk;
		chunk.resize(primes.size()*sizeof(std::uint64_t));
		char*dest=chunk.data();
		for(std::uint64_t value : primes){
			std::uint64_t encoded=to_little_endian(value);
			std::memcpy(dest,&encoded,sizeof(encoded));
			dest+=sizeof(encoded);
		}
		enqueue_chunk(Chunk{std::move(chunk),false});
		break;
	}
	case PrimeOutputFormat::ZstdDelta:{
		std::string data=encode_deltas(primes);
		if(!data.empty()){
			enqueue_chunk(Chunk{std::move(data),false});
		}
		break;
	}
	}
}

void PrimeWriter::write_value(std::uint64_t value){
	if(!enabled_){
		return;
	}
	switch(format_){
	case PrimeOutputFormat::Text:{
		char local[32];
		auto result=std::to_chars(local,local+sizeof(local),value);
		if(result.ec!=std::errc()){
			throw std::runtime_error("Failed to convert prime to string");
		}
		std::string chunk(local,result.ptr);
		chunk.push_back('\n');
		enqueue_chunk(Chunk{std::move(chunk),false});
		break;
	}
	case PrimeOutputFormat::Binary:{
		std::uint64_t encoded=to_little_endian(value);
		std::string chunk(reinterpret_cast<const char*>(&encoded),
						  sizeof(encoded));
		enqueue_chunk(Chunk{std::move(chunk),false});
		break;
	}
	case PrimeOutputFormat::ZstdDelta:{
		std::string data=encode_delta_value(value);
		if(!data.empty()){
			enqueue_chunk(Chunk{std::move(data),false});
		}
		break;
	}
	}
}

void PrimeWriter::flush(){
	if(!enabled_){
		return;
	}
	enqueue_chunk(Chunk{{},true});
}

void PrimeWriter::finish(){
	if(!enabled_){
		return;
	}

	bool already_stopped=false;
	{
		std::lock_guard<std::mutex> lock(queue_mutex_);
		already_stopped=stop_requested_;
	}

	std::exception_ptr flush_error;
	if(!already_stopped){
		try{
			flush();
		}catch(...){
			flush_error=std::current_exception();
		}
		{
			std::lock_guard<std::mutex> lock(queue_mutex_);
			stop_requested_=true;
		}
		queue_not_empty_.notify_one();
	}

	if(writer_thread_.joinable()){
		writer_thread_.join();
	}

	if(file_){
		if(owns_file_){
			if(std::fclose(file_)!=0){
				if(!flush_error){
					flush_error=std::make_exception_ptr(
						std::runtime_error("Failed to close output file"));
				}
			}
		}else{
			if(std::fflush(file_)!=0){
				if(!flush_error){
					flush_error=std::make_exception_ptr(
						std::runtime_error("Failed to flush output stream"));
				}
			}
		}
		file_=nullptr;
	}

	if(flush_error){
		std::rethrow_exception(flush_error);
	}

	check_io_error();
}

void PrimeWriter::enqueue_chunk(Chunk&&chunk){
	if(!enabled_){
		return;
	}

	check_io_error();

	std::unique_lock<std::mutex> lock(queue_mutex_);
	queue_not_full_.wait(
		lock,[&]{ return queue_.size()<queue_capacity_||stop_requested_; });
	if(stop_requested_){
		throw std::runtime_error("Writer has been stopped");
	}
	queue_.push_back(std::move(chunk));
	lock.unlock();
	queue_not_empty_.notify_one();
}

void PrimeWriter::writer_loop(){
	for(;;){
		Chunk chunk;
		{
			std::unique_lock<std::mutex> lock(queue_mutex_);
			queue_not_empty_.wait(
				lock,[&]{ return stop_requested_||!queue_.empty(); });
			if(queue_.empty()){
				if(stop_requested_){
					break;
				}
				continue;
			}
			chunk=std::move(queue_.front());
			queue_.pop_front();
			queue_not_full_.notify_one();
		}

		if(!chunk.data.empty()){
			buffer_.append(chunk.data);
			if(buffer_.size()>=buffer_threshold_){
				flush_buffer();
			}
		}
		if(chunk.flush){
			flush_buffer();
			if(file_&&std::fflush(file_)!=0){
				set_error(std::strerror(errno));
			}
		}
	}

	flush_buffer();
	if(file_&&std::fflush(file_)!=0){
		set_error(std::strerror(errno));
	}
}

void PrimeWriter::flush_buffer(){
	if(!file_||buffer_.empty()){
		return;
	}
	const char*data=buffer_.data();
	std::size_t remaining=buffer_.size();
	while(remaining>0){
		std::size_t written=std::fwrite(data,1,remaining,file_);
		if(written==0){
			if(std::ferror(file_)){
				set_error(std::strerror(errno));
			}
			break;
		}
		data+=written;
		remaining-=written;
	}
	if(remaining==0){
		buffer_.clear();
	}
}

void PrimeWriter::check_io_error() const{
	if(!io_error_.load(std::memory_order_acquire)){
		return;
	}
	std::lock_guard<std::mutex> lock(error_mutex_);
	throw std::runtime_error(error_message_.empty()?"I/O error":error_message_);
}

void PrimeWriter::set_error(const std::string&message){
	bool expected=false;
	if(io_error_.compare_exchange_strong(expected,true,
										 std::memory_order_acq_rel)){
		std::lock_guard<std::mutex> lock(error_mutex_);
		error_message_=message;
	}
}

std::string PrimeWriter::encode_deltas(const std::vector<std::uint64_t>&primes){
	if(format_!=PrimeOutputFormat::ZstdDelta){
		return {};
	}
	std::string raw;
	raw.resize(primes.size()*sizeof(std::uint64_t));
	char*dest=raw.data();
	for(std::uint64_t value : primes){
		if(value<previous_prime_){
			throw std::runtime_error(
				"Primes must be non-decreasing for delta encoding");
		}
		std::uint64_t delta=value-previous_prime_;
		previous_prime_=value;
		std::uint64_t encoded=to_little_endian(delta);
		std::memcpy(dest,&encoded,sizeof(encoded));
		dest+=sizeof(encoded);
	}
	return raw;
}

std::string PrimeWriter::encode_delta_value(std::uint64_t value){
	if(format_!=PrimeOutputFormat::ZstdDelta){
		return {};
	}
	if(value<previous_prime_){
		throw std::runtime_error(
			"Primes must be non-decreasing for delta encoding");
	}
	std::uint64_t delta=value-previous_prime_;
	previous_prime_=value;
	std::uint64_t encoded=to_little_endian(delta);
	return std::string(reinterpret_cast<const char*>(&encoded),sizeof(encoded));
}

} // namespace calcprime
