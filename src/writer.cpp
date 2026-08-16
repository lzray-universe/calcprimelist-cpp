#include "writer.h"
#include "parquet_format.h"

#include<algorithm>
#include<cerrno>
#include<charconv>
#include<climits>
#include<cstring>
#include<exception>
#include<limits>
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
constexpr std::size_t kParquetRowsPerGroup=1u<<20;
constexpr char kParquetMagic[]={'P','A','R','1'};

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
						 PrimeOutputFormat format,bool use_zstd,
						 ParquetEncoding parquet_encoding,
						 std::size_t parquet_delta_block_values)
	: enabled_(enabled),file_(nullptr),owns_file_(false),
	  queue_capacity_(kDefaultQueueCapacity),stop_requested_(false),
	  buffer_threshold_(kDefaultBufferThreshold),format_(format),
	  use_zstd_(use_zstd),parquet_encoding_(parquet_encoding),
	  parquet_delta_block_values_(parquet_delta_block_values),
	  has_first_prime_(false),previous_prime_(0),
	  zstd_cctx_(nullptr),file_offset_(0),parquet_num_rows_(0),
	  parquet_footer_written_(false),io_error_(false){
	if(!enabled_){
		return;
	}
	if(parquet_encoding_==ParquetEncoding::DeltaBinaryPacked&&
	   (parquet_delta_block_values_==0||
		(parquet_delta_block_values_%128U)!=0U||
		parquet_delta_block_values_>kParquetRowsPerGroup)){
		throw std::invalid_argument(
			"Parquet delta block values must be a multiple of 128 between 128 and 1048576");
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
	if(format_==PrimeOutputFormat::Parquet){
		write_file_bytes(kParquetMagic,sizeof(kParquetMagic));
		check_io_error();
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
	case PrimeOutputFormat::Parquet:{
		for(std::size_t offset=0;offset<primes.size();){
			std::size_t count=std::min(kParquetRowsPerGroup,
								 primes.size()-offset);
			std::string data;
			if(parquet_encoding_==ParquetEncoding::DeltaBinaryPacked){
				data=parquet::encode_delta_binary_packed(
					primes.data()+offset,count,parquet_delta_block_values_);
			}else{
				data.resize(count*sizeof(std::uint64_t));
				char*dest=data.data();
				for(std::size_t i=0;i<count;++i){
					std::uint64_t encoded=
						to_little_endian_u64(primes[offset+i]);
					std::memcpy(dest,&encoded,sizeof(encoded));
					dest+=sizeof(encoded);
				}
			}
			enqueue_chunk(Chunk{std::move(data),false,
								static_cast<std::uint64_t>(count)});
			offset+=count;
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
	case PrimeOutputFormat::Parquet:{
		std::string chunk;
		if(parquet_encoding_==ParquetEncoding::DeltaBinaryPacked){
			chunk=parquet::encode_delta_binary_packed(
				&value,1,parquet_delta_block_values_);
		}else{
			std::uint64_t encoded=to_little_endian_u64(value);
			chunk.assign(reinterpret_cast<const char*>(&encoded),sizeof(encoded));
		}
		enqueue_chunk(Chunk{std::move(chunk),false,1});
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

		if(format_==PrimeOutputFormat::Parquet&&!chunk.data.empty()){
			write_parquet_chunk(chunk);
		}else if(!chunk.data.empty()){
			buffer_.append(chunk.data);
			if(buffer_.size()>=buffer_threshold_){
				flush_buffer();
			}
		}
		if(chunk.flush){
			if(format_!=PrimeOutputFormat::Parquet){
				flush_buffer();
#if defined(CALCPRIME_HAS_ZSTD)
				if(use_zstd_){
					flush_zstd_stream(false);
				}
#endif
			}
			if(file_&&std::fflush(file_)!=0){
				set_error(std::strerror(errno));
			}
		}
	}

	if(format_==PrimeOutputFormat::Parquet){
		write_parquet_footer();
	}else{
		flush_buffer();
	}
#if defined(CALCPRIME_HAS_ZSTD)
	if(use_zstd_&&format_!=PrimeOutputFormat::Parquet){
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
		file_offset_+=static_cast<std::uint64_t>(written);
	}
}

void PrimeWriter::write_parquet_chunk(const Chunk&chunk){
	if(chunk.value_count==0||chunk.data.empty()){
		return;
	}
	if(chunk.value_count>
	   static_cast<std::uint64_t>(std::numeric_limits<std::int32_t>::max())||
	   chunk.data.size()>
		   static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max())){
		set_error("Parquet data page exceeds the supported 2 GiB limit");
		return;
	}

	const std::string*payload=&chunk.data;
	std::string compressed;
#if defined(CALCPRIME_HAS_ZSTD)
	if(use_zstd_){
		std::size_t bound=ZSTD_compressBound(chunk.data.size());
		compressed.resize(bound);
		std::size_t result=ZSTD_compress2(
			static_cast<ZSTD_CCtx*>(zstd_cctx_),compressed.data(),
			compressed.size(),chunk.data.data(),chunk.data.size());
		if(ZSTD_isError(result)){
			std::string message="zstd compress error: ";
			message.append(ZSTD_getErrorName(result));
			set_error(message);
			return;
		}
		compressed.resize(result);
		payload=&compressed;
	}
#else
	if(use_zstd_){
		set_error("zstd not supported in this build");
		return;
	}
#endif

	if(payload->size()>
	   static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max())){
		set_error("Compressed Parquet data page exceeds the supported 2 GiB limit");
		return;
	}

	std::string header=parquet::make_data_page_header(
		static_cast<std::int32_t>(chunk.value_count),
		static_cast<std::int32_t>(chunk.data.size()),
		static_cast<std::int32_t>(payload->size()),
		parquet_encoding_==ParquetEncoding::DeltaBinaryPacked
			?parquet::ValueEncoding::DeltaBinaryPacked
			:parquet::ValueEncoding::Plain);
	std::uint64_t page_offset=file_offset_;
	write_file_bytes(header.data(),header.size());
	write_file_bytes(payload->data(),payload->size());
	if(io_error_.load(std::memory_order_acquire)){
		return;
	}

	ParquetRowGroup row_group;
	row_group.data_page_offset=static_cast<std::int64_t>(page_offset);
	row_group.num_values=static_cast<std::int64_t>(chunk.value_count);
	row_group.total_uncompressed_size=static_cast<std::int64_t>(
		header.size()+chunk.data.size());
	row_group.total_compressed_size=static_cast<std::int64_t>(
		header.size()+payload->size());
	parquet_row_groups_.push_back(row_group);
	parquet_num_rows_+=chunk.value_count;
}

void PrimeWriter::write_parquet_footer(){
	if(parquet_footer_written_||io_error_.load(std::memory_order_acquire)){
		return;
	}
	std::vector<parquet::RowGroupMetadata> metadata;
	metadata.reserve(parquet_row_groups_.size());
	for(const ParquetRowGroup&row_group : parquet_row_groups_){
		metadata.push_back(parquet::RowGroupMetadata{
			row_group.data_page_offset,row_group.num_values,
			row_group.total_uncompressed_size,row_group.total_compressed_size});
	}
	std::string footer=parquet::make_file_metadata(
		metadata,static_cast<std::int64_t>(parquet_num_rows_),use_zstd_,
		parquet_encoding_==ParquetEncoding::DeltaBinaryPacked
			?parquet::ValueEncoding::DeltaBinaryPacked
			:parquet::ValueEncoding::Plain);
	if(footer.size()>std::numeric_limits<std::uint32_t>::max()){
		set_error("Parquet footer exceeds the 4 GiB format limit");
		return;
	}
	write_file_bytes(footer.data(),footer.size());
	std::uint32_t footer_size=static_cast<std::uint32_t>(footer.size());
	char length[4]={
		static_cast<char>(footer_size&0xffU),
		static_cast<char>((footer_size>>8U)&0xffU),
		static_cast<char>((footer_size>>16U)&0xffU),
		static_cast<char>((footer_size>>24U)&0xffU),
	};
	write_file_bytes(length,sizeof(length));
	write_file_bytes(kParquetMagic,sizeof(kParquetMagic));
	parquet_footer_written_=!io_error_.load(std::memory_order_acquire);
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
