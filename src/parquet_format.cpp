#include "parquet_format.h"

#include<algorithm>
#include<bit>
#include<cstddef>
#include<cstdint>
#include<limits>
#include<stdexcept>
#include<string>
#include<vector>

namespace calcprime::parquet{

namespace{

// Parquet metadata and page headers use the Thrift Compact Protocol.  Keeping
// this tiny encoder local avoids adding Arrow/Thrift as build dependencies for
// the single required uint64 column emitted by calcprimelist.
enum CompactType : std::uint8_t{
	CompactStop=0x00,
	CompactBooleanTrue=0x01,
	CompactBooleanFalse=0x02,
	CompactByte=0x03,
	CompactI16=0x04,
	CompactI32=0x05,
	CompactI64=0x06,
	CompactBinary=0x08,
	CompactList=0x09,
	CompactStruct=0x0c,
};

void append_uvarint(std::string&out,std::uint64_t value){
	while(value>=0x80U){
		out.push_back(static_cast<char>((value&0x7fU)|0x80U));
		value>>=7U;
	}
	out.push_back(static_cast<char>(value));
}

std::uint64_t zigzag_i64(std::int64_t value){
	return (static_cast<std::uint64_t>(value)<<1U)^
		   static_cast<std::uint64_t>(-(value<0));
}

void append_field_header(std::string&out,std::int16_t&last_field,
						 std::int16_t field,CompactType type){
	std::int16_t delta=static_cast<std::int16_t>(field-last_field);
	if(delta>0&&delta<=15){
		out.push_back(static_cast<char>((delta<<4U)|type));
	}else{
		out.push_back(static_cast<char>(type));
		append_uvarint(out,zigzag_i64(field));
	}
	last_field=field;
}

void append_i32_field(std::string&out,std::int16_t&last_field,
					  std::int16_t field,std::int32_t value){
	append_field_header(out,last_field,field,CompactI32);
	append_uvarint(out,zigzag_i64(value));
}

void append_i64_field(std::string&out,std::int16_t&last_field,
					  std::int16_t field,std::int64_t value){
	append_field_header(out,last_field,field,CompactI64);
	append_uvarint(out,zigzag_i64(value));
}

void append_byte_field(std::string&out,std::int16_t&last_field,
					   std::int16_t field,std::uint8_t value){
	append_field_header(out,last_field,field,CompactByte);
	out.push_back(static_cast<char>(value));
}

void append_bool_field(std::string&out,std::int16_t&last_field,
					   std::int16_t field,bool value){
	append_field_header(out,last_field,field,
						value?CompactBooleanTrue:CompactBooleanFalse);
}

void append_binary_field(std::string&out,std::int16_t&last_field,
						 std::int16_t field,const std::string&value){
	append_field_header(out,last_field,field,CompactBinary);
	append_uvarint(out,static_cast<std::uint64_t>(value.size()));
	out.append(value);
}

void append_list_header(std::string&out,std::size_t size,
						CompactType element_type){
	if(size<=14){
		out.push_back(static_cast<char>((size<<4U)|element_type));
	}else{
		out.push_back(static_cast<char>(0xf0U|element_type));
		append_uvarint(out,static_cast<std::uint64_t>(size));
	}
}

void append_stop(std::string&out){
	out.push_back(static_cast<char>(CompactStop));
}

void append_int_type(std::string&out){
	std::int16_t last=0;
	append_byte_field(out,last,1,64);
	append_bool_field(out,last,2,false);
	append_stop(out);
}

void append_logical_type(std::string&out){
	std::int16_t last=0;
	append_field_header(out,last,10,CompactStruct); // INTEGER
	append_int_type(out);
	append_stop(out);
}

void append_root_schema(std::string&out){
	std::int16_t last=0;
	append_binary_field(out,last,4,"schema");
	append_i32_field(out,last,5,1);
	append_stop(out);
}

void append_prime_schema(std::string&out){
	std::int16_t last=0;
	append_i32_field(out,last,1,2); // Type::INT64
	append_i32_field(out,last,3,0); // FieldRepetitionType::REQUIRED
	append_binary_field(out,last,4,"prime");
	append_i32_field(out,last,6,14); // ConvertedType::UINT_64
	append_field_header(out,last,10,CompactStruct);
	append_logical_type(out);
	append_stop(out);
}

void append_column_metadata(std::string&out,const RowGroupMetadata&row_group,
							bool use_zstd,ValueEncoding encoding){
	std::int16_t last=0;
	append_i32_field(out,last,1,2); // Type::INT64

	append_field_header(out,last,2,CompactList);
	append_list_header(out,2,CompactI32);
	append_uvarint(out,zigzag_i64(static_cast<std::int32_t>(encoding)));
	append_uvarint(out,zigzag_i64(3)); // Encoding::RLE (level encoding)

	append_field_header(out,last,3,CompactList);
	append_list_header(out,1,CompactBinary);
	append_uvarint(out,5);
	out.append("prime");

	append_i32_field(out,last,4,use_zstd?6:0); // CompressionCodec
	append_i64_field(out,last,5,row_group.num_values);
	append_i64_field(out,last,6,row_group.total_uncompressed_size);
	append_i64_field(out,last,7,row_group.total_compressed_size);
	append_i64_field(out,last,9,row_group.data_page_offset);
	append_stop(out);
}

void append_column_chunk(std::string&out,const RowGroupMetadata&row_group,
						 bool use_zstd,ValueEncoding encoding){
	std::int16_t last=0;
	// parquet.thrift requires this deprecated field and recommends zero when
	// column metadata lives only in the footer.
	append_i64_field(out,last,2,0);
	append_field_header(out,last,3,CompactStruct);
	append_column_metadata(out,row_group,use_zstd,encoding);
	append_stop(out);
}

void append_row_group(std::string&out,const RowGroupMetadata&row_group,
					  bool use_zstd,ValueEncoding encoding){
	std::int16_t last=0;
	append_field_header(out,last,1,CompactList);
	append_list_header(out,1,CompactStruct);
	append_column_chunk(out,row_group,use_zstd,encoding);
	append_i64_field(out,last,2,row_group.total_uncompressed_size);
	append_i64_field(out,last,3,row_group.num_values);
	append_i64_field(out,last,5,row_group.data_page_offset);
	append_i64_field(out,last,6,row_group.total_compressed_size);
	append_stop(out);
}

void append_data_page_header(std::string&out,std::int32_t num_values,
							 ValueEncoding encoding){
	std::int16_t last=0;
	append_i32_field(out,last,1,num_values);
	append_i32_field(out,last,2,static_cast<std::int32_t>(encoding));
	append_i32_field(out,last,3,3); // Encoding::RLE
	append_i32_field(out,last,4,3); // Encoding::RLE
	append_stop(out);
}

unsigned bit_width(std::uint64_t value){
	return value==0?0U:64U-std::countl_zero(value);
}

void append_bit_packed(std::string&out,const std::uint64_t*values,
					   std::size_t count,unsigned width){
	if(width==0){
		return;
	}
	std::uint64_t buffer=0;
	unsigned buffered_bits=0;
	for(std::size_t i=0;i<count;++i){
		std::uint64_t value=values[i];
		unsigned remaining=width;
		while(remaining>0){
			unsigned room=64U-buffered_bits;
			unsigned take=std::min(remaining,room);
			std::uint64_t mask=
				(take==64U)?std::numeric_limits<std::uint64_t>::max()
						   :((std::uint64_t{1}<<take)-1U);
			buffer|=(value&mask)<<buffered_bits;
			if(take<64U){
				value>>=take;
			}
			buffered_bits+=take;
			remaining-=take;
			while(buffered_bits>=8U){
				out.push_back(static_cast<char>(buffer&0xffU));
				buffer>>=8U;
				buffered_bits-=8U;
			}
		}
	}
	if(buffered_bits!=0U){
		out.push_back(static_cast<char>(buffer&0xffU));
	}
}

template<class Delta>
void append_delta_block(std::string&out,const std::vector<Delta>&deltas,
						std::size_t mini_block_count,
						std::size_t values_per_mini_block){
	Delta min_delta=*std::min_element(deltas.begin(),deltas.end());
	append_uvarint(out,zigzag_i64(static_cast<std::int64_t>(min_delta)));

	std::vector<std::uint8_t> widths(mini_block_count,0);
	std::size_t used_mini_blocks=
		(deltas.size()+values_per_mini_block-1U)/values_per_mini_block;
	for(std::size_t mini=0;mini<used_mini_blocks;++mini){
		std::size_t begin=mini*values_per_mini_block;
		std::size_t end=std::min(begin+values_per_mini_block,deltas.size());
		std::uint64_t max_adjusted=0;
		for(std::size_t i=begin;i<end;++i){
			max_adjusted=std::max(
				max_adjusted,static_cast<std::uint64_t>(deltas[i]-min_delta));
		}
		widths[mini]=static_cast<std::uint8_t>(bit_width(max_adjusted));
	}
	out.append(reinterpret_cast<const char*>(widths.data()),widths.size());

	std::vector<std::uint64_t> adjusted(values_per_mini_block,0);
	for(std::size_t mini=0;mini<used_mini_blocks;++mini){
		std::fill(adjusted.begin(),adjusted.end(),0);
		std::size_t begin=mini*values_per_mini_block;
		std::size_t end=std::min(begin+values_per_mini_block,deltas.size());
		for(std::size_t i=begin;i<end;++i){
			adjusted[i-begin]=
				static_cast<std::uint64_t>(deltas[i]-min_delta);
		}
		append_bit_packed(out,adjusted.data(),adjusted.size(),widths[mini]);
	}
}

} // namespace

std::string make_data_page_header(std::int32_t num_values,
								  std::int32_t uncompressed_size,
								  std::int32_t compressed_size,
								  ValueEncoding encoding){
	std::string out;
	out.reserve(32);
	std::int16_t last=0;
	append_i32_field(out,last,1,0); // PageType::DATA_PAGE
	append_i32_field(out,last,2,uncompressed_size);
	append_i32_field(out,last,3,compressed_size);
	append_field_header(out,last,5,CompactStruct);
	append_data_page_header(out,num_values,encoding);
	append_stop(out);
	return out;
}

std::string make_file_metadata(const std::vector<RowGroupMetadata>&row_groups,
							   std::int64_t num_rows,bool use_zstd,
							   ValueEncoding encoding){
	std::string out;
	out.reserve(128+row_groups.size()*64);
	std::int16_t last=0;
	append_i32_field(out,last,1,1); // File format version: always 1

	append_field_header(out,last,2,CompactList);
	append_list_header(out,2,CompactStruct);
	append_root_schema(out);
	append_prime_schema(out);

	append_i64_field(out,last,3,num_rows);
	append_field_header(out,last,4,CompactList);
	append_list_header(out,row_groups.size(),CompactStruct);
	for(const RowGroupMetadata&row_group : row_groups){
		append_row_group(out,row_group,use_zstd,encoding);
	}

	append_binary_field(out,last,6,"calcprimelist-cpp");
	append_stop(out);
	return out;
}

std::string encode_delta_binary_packed(const std::uint64_t*values,
									   std::size_t count,
									   std::size_t block_value_count){
	if(block_value_count==0||(block_value_count%128U)!=0U){
		throw std::invalid_argument(
			"Parquet delta block value count must be a positive multiple of 128");
	}
	if(count==0){
		return {};
	}
	if(!values){
		throw std::invalid_argument("Parquet delta values pointer is null");
	}
	if(values[0]>
	   static_cast<std::uint64_t>(std::numeric_limits<std::int64_t>::max())){
		throw std::overflow_error(
			"DELTA_BINARY_PACKED requires values within signed INT64 range");
	}

	constexpr std::size_t values_per_mini_block=32;
	std::size_t mini_block_count=block_value_count/values_per_mini_block;
	std::string out;
	out.reserve(count*2U+32U);
	append_uvarint(out,static_cast<std::uint64_t>(block_value_count));
	append_uvarint(out,static_cast<std::uint64_t>(mini_block_count));
	append_uvarint(out,static_cast<std::uint64_t>(count));
	append_uvarint(out,zigzag_i64(static_cast<std::int64_t>(values[0])));

	std::size_t offset=1;
	while(offset<count){
		std::size_t delta_count=std::min(block_value_count,count-offset);
		bool fits_u16=true;
		std::vector<std::uint16_t> deltas16;
		deltas16.reserve(delta_count);
		std::uint64_t previous=values[offset-1U];
		for(std::size_t i=0;i<delta_count;++i){
			std::uint64_t current=values[offset+i];
			if(current>static_cast<std::uint64_t>(
						 std::numeric_limits<std::int64_t>::max())){
				throw std::overflow_error(
					"DELTA_BINARY_PACKED requires values within signed INT64 range");
			}
			if(current<previous){
				throw std::invalid_argument(
					"Parquet DELTA_BINARY_PACKED values must be non-decreasing");
			}
			std::uint64_t delta=current-previous;
			if(delta>static_cast<std::uint64_t>(
						 std::numeric_limits<std::int64_t>::max())){
				throw std::overflow_error(
					"Parquet DELTA_BINARY_PACKED delta exceeds INT64 range");
			}
			if(delta>std::numeric_limits<std::uint16_t>::max()){
				fits_u16=false;
			}
			if(fits_u16){
				deltas16.push_back(static_cast<std::uint16_t>(delta));
			}
			previous=current;
		}

		if(fits_u16){
			append_delta_block(out,deltas16,mini_block_count,
							   values_per_mini_block);
		}else{
			std::vector<std::uint64_t> deltas64;
			deltas64.reserve(delta_count);
			previous=values[offset-1U];
			for(std::size_t i=0;i<delta_count;++i){
				std::uint64_t current=values[offset+i];
				deltas64.push_back(current-previous);
				previous=current;
			}
			append_delta_block(out,deltas64,mini_block_count,
							   values_per_mini_block);
		}
		offset+=delta_count;
	}
	return out;
}

} // namespace calcprime::parquet
