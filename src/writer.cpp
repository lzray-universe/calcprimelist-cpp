#include "writer.h"

#include<cerrno>
#include<charconv>
#include<climits>
#include<cstring>
#include<exception>
#include<stdexcept>

#if defined(_MSC_VER)
#include<intrin.h>
#endif

#if defined(CALCPRIME_HAS_ZSTD)
#include<zstd.h>
#endif

namespace calcprime{

namespace{

constexpr std::size_t kDefaultFileBuffer=8u<<20; // 8 MiB
constexpr std::size_t kDefaultQueueCapacity=8;
constexpr std::size_t kDefaultBufferThreshold=8u<<20; // 8 MiB

inline std::uint64_t to_little_endian_u64(std::uint64_t value){
#if defined(__BYTE_ORDER__)&&(__BYTE_ORDER__==__ORDER_BIG_ENDIAN__)
	return __builtin_bswap64(value);
#else
	return value;
#endif
}

inline std::uint16_t to_little_endian_u16(std::uint16_t value){
#if defined(__BYTE_ORDER__)&&(__BYTE_ORDER__==__ORDER_BIG_ENDIAN__)
	return static_cast<std::uint16_t>((value>>8)|(value<<8));
#else
	return value;
#endif
}

} // namespace

PrimeWriter::PrimeWriter(bool enabled,const std::string&path,
						 PrimeOutputFormat format,bool use_zstd)
	: enabled_(enabled),file_(nullptr),owns_file_(false),
	  queue_capacity_(kDefaultQueueCapacity),stop_requested_(false),
	  buffer_threshold_(kDefaultBufferThreshold),format_(format),
	  use_zstd_(use_zstd),has_first_prime_(false),previous_prime_(0),
	  zstd_cctx_(nullptr),io_error_(false){
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

	if(use_zstd_){
#if defined(CALCPRIME_HAS_ZSTD)
		ZSTD_CCtx*cctx=ZSTD_createCCtx();
		if(!cctx){
			throw std::runtime_error("Failed to create zstd context");
		}
		std::size_t configured=ZSTD_CCtx_setParameter(
			cctx,ZSTD_c_compressionLevel,1);
		if(ZSTD_isError(configured)){
			std::string message="Failed to configure zstd context: ";
			message.append(ZSTD_getErrorName(configured));
			ZSTD_freeCCtx(cctx);
			throw std::runtime_error(message);
		}
		zstd_cctx_=cctx;
		zstd_out_buffer_.resize(ZSTD_CStreamOutSize());
#else
		throw std::runtime_error("zstd not supported in this build");
#endif
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
			std::uint64_t encoded=to_little_endian_u64(value);
			std::memcpy(dest,&encoded,sizeof(encoded));
			dest+=sizeof(encoded);
		}
		enqueue_chunk(Chunk{std::move(chunk),false});
		break;
	}
	case PrimeOutputFormat::Delta16:{
		std::string data=encode_delta16(primes);
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
		std::uint64_t encoded=to_little_endian_u64(value);
		std::string chunk(reinterpret_cast<const char*>(&encoded),
						  sizeof(encoded));
		enqueue_chunk(Chunk{std::move(chunk),false});
		break;
	}
	case PrimeOutputFormat::Delta16:{
		std::string data=encode_delta16_value(value);
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

#if defined(CALCPRIME_HAS_ZSTD)
	if(zstd_cctx_){
		ZSTD_freeCCtx(static_cast<ZSTD_CCtx*>(zstd_cctx_));
		zstd_cctx_=nullptr;
	}
#endif

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
#if defined(CALCPRIME_HAS_ZSTD)
			if(use_zstd_){
				flush_zstd_stream(false);
			}
#endif
			if(file_&&std::fflush(file_)!=0){
				set_error(std::strerror(errno));
			}
		}
	}

	flush_buffer();
#if defined(CALCPRIME_HAS_ZSTD)
	if(use_zstd_){
		flush_zstd_stream(true);
	}
#endif
	if(file_&&std::fflush(file_)!=0){
		set_error(std::strerror(errno));
	}
}

void PrimeWriter::flush_buffer(){
	if(!file_||buffer_.empty()){
		return;
	}

	if(!use_zstd_){
		write_file_bytes(buffer_.data(),buffer_.size());
		if(!io_error_.load(std::memory_order_acquire)){
			buffer_.clear();
		}
		return;
	}

#if defined(CALCPRIME_HAS_ZSTD)
	if(!zstd_cctx_){
		set_error("zstd context is not initialized");
		return;
	}
	ZSTD_inBuffer input{buffer_.data(),buffer_.size(),0};
	while(input.pos<input.size){
		ZSTD_outBuffer output{zstd_out_buffer_.data(),zstd_out_buffer_.size(),0};
		std::size_t code=ZSTD_compressStream2(
			static_cast<ZSTD_CCtx*>(zstd_cctx_),&output,&input,ZSTD_e_continue);
		if(ZSTD_isError(code)){
			std::string message="zstd compress error: ";
			message.append(ZSTD_getErrorName(code));
			set_error(message);
			break;
		}
		if(output.pos>0){
			write_file_bytes(zstd_out_buffer_.data(),output.pos);
			if(io_error_.load(std::memory_order_acquire)){
				break;
			}
		}
	}
	if(!io_error_.load(std::memory_order_acquire)){
		buffer_.clear();
	}
#else
	set_error("zstd not supported in this build");
#endif
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

std::string PrimeWriter::encode_delta16(
	const std::vector<std::uint64_t>&primes){
	if(format_!=PrimeOutputFormat::Delta16||primes.empty()){
		return {};
	}

	bool local_has_first=has_first_prime_;
	std::uint64_t local_previous=previous_prime_;
	std::size_t bytes=primes.size()*sizeof(std::int16_t);
	if(!local_has_first){
		bytes=sizeof(std::uint64_t);
		if(primes.size()>1){
			bytes+=(primes.size()-1)*sizeof(std::int16_t);
		}
	}

	std::string encoded;
	encoded.resize(bytes);
	char*dest=encoded.data();

	std::size_t index=0;
	if(!local_has_first){
		std::uint64_t first_prime=primes[0];
		std::uint64_t first_encoded=to_little_endian_u64(first_prime);
		std::memcpy(dest,&first_encoded,sizeof(first_encoded));
		dest+=sizeof(first_encoded);
		local_previous=first_prime;
		local_has_first=true;
		index=1;
	}

	for(;index<primes.size();++index){
		std::uint64_t value=primes[index];
		if(value<local_previous){
			throw std::runtime_error(
				"Primes must be non-decreasing for delta16 encoding");
		}
		std::uint64_t delta=value-local_previous;
		if(delta==0){
			throw std::runtime_error(
				"Prime delta must be positive for delta16 encoding");
		}
		if(delta>static_cast<std::uint64_t>(INT16_MAX)){
			throw std::runtime_error(
				"Prime delta exceeds int16 range in delta16 output");
		}
		std::int16_t signed_delta=static_cast<std::int16_t>(delta);
		std::uint16_t delta_encoded=
			to_little_endian_u16(static_cast<std::uint16_t>(signed_delta));
		std::memcpy(dest,&delta_encoded,sizeof(delta_encoded));
		dest+=sizeof(delta_encoded);
		local_previous=value;
	}

	has_first_prime_=local_has_first;
	previous_prime_=local_previous;
	return encoded;
}

std::string PrimeWriter::encode_delta16_value(std::uint64_t value){
	if(format_!=PrimeOutputFormat::Delta16){
		return {};
	}

	if(!has_first_prime_){
		has_first_prime_=true;
		previous_prime_=value;
		std::uint64_t encoded=to_little_endian_u64(value);
		return std::string(reinterpret_cast<const char*>(&encoded),
						   sizeof(encoded));
	}

	if(value<previous_prime_){
		throw std::runtime_error(
			"Primes must be non-decreasing for delta16 encoding");
	}

	std::uint64_t delta=value-previous_prime_;
	if(delta==0){
		throw std::runtime_error(
			"Prime delta must be positive for delta16 encoding");
	}
	if(delta>static_cast<std::uint64_t>(INT16_MAX)){
		throw std::runtime_error(
			"Prime delta exceeds int16 range in delta16 output");
	}

	previous_prime_=value;
	std::int16_t signed_delta=static_cast<std::int16_t>(delta);
	std::uint16_t encoded=
		to_little_endian_u16(static_cast<std::uint16_t>(signed_delta));
	return std::string(reinterpret_cast<const char*>(&encoded),sizeof(encoded));
}

void PrimeWriter::write_file_bytes(const char*data,std::size_t size){
	if(!file_||size==0){
		return;
	}
	const char*cursor=data;
	std::size_t remaining=size;
	while(remaining>0){
		std::size_t written=std::fwrite(cursor,1,remaining,file_);
		if(written==0){
			if(std::ferror(file_)){
				set_error(std::strerror(errno));
			}
			break;
		}
		cursor+=written;
		remaining-=written;
	}
}

#if defined(CALCPRIME_HAS_ZSTD)
void PrimeWriter::flush_zstd_stream(bool final_frame){
	if(!use_zstd_){
		return;
	}
	if(!zstd_cctx_){
		set_error("zstd context is not initialized");
		return;
	}
	ZSTD_EndDirective mode=final_frame?ZSTD_e_end:ZSTD_e_flush;
	ZSTD_inBuffer input{nullptr,0,0};
	for(;;){
		ZSTD_outBuffer output{zstd_out_buffer_.data(),zstd_out_buffer_.size(),0};
		std::size_t code=ZSTD_compressStream2(
			static_cast<ZSTD_CCtx*>(zstd_cctx_),&output,&input,mode);
		if(ZSTD_isError(code)){
			std::string message="zstd compress error: ";
			message.append(ZSTD_getErrorName(code));
			set_error(message);
			return;
		}
		if(output.pos>0){
			write_file_bytes(zstd_out_buffer_.data(),output.pos);
			if(io_error_.load(std::memory_order_acquire)){
				return;
			}
		}
		if(code==0){
			return;
		}
	}
}
#endif

} // namespace calcprime
