#include "parquet_format.h"

#include<cstddef>
#include<cstdint>
#include<string>

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
							bool use_zstd){
	std::int16_t last=0;
	append_i32_field(out,last,1,2); // Type::INT64

	append_field_header(out,last,2,CompactList);
	append_list_header(out,2,CompactI32);
	append_uvarint(out,zigzag_i64(0)); // Encoding::PLAIN
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
						 bool use_zstd){
	std::int16_t last=0;
	// parquet.thrift requires this deprecated field and recommends zero when
	// column metadata lives only in the footer.
	append_i64_field(out,last,2,0);
	append_field_header(out,last,3,CompactStruct);
	append_column_metadata(out,row_group,use_zstd);
	append_stop(out);
}

void append_row_group(std::string&out,const RowGroupMetadata&row_group,
					  bool use_zstd){
	std::int16_t last=0;
	append_field_header(out,last,1,CompactList);
	append_list_header(out,1,CompactStruct);
	append_column_chunk(out,row_group,use_zstd);
	append_i64_field(out,last,2,row_group.total_uncompressed_size);
	append_i64_field(out,last,3,row_group.num_values);
	append_i64_field(out,last,5,row_group.data_page_offset);
	append_i64_field(out,last,6,row_group.total_compressed_size);
	append_stop(out);
}

void append_data_page_header(std::string&out,std::int32_t num_values){
	std::int16_t last=0;
	append_i32_field(out,last,1,num_values);
	append_i32_field(out,last,2,0); // Encoding::PLAIN
	append_i32_field(out,last,3,3); // Encoding::RLE
	append_i32_field(out,last,4,3); // Encoding::RLE
	append_stop(out);
}

} // namespace

std::string make_data_page_header(std::int32_t num_values,
								  std::int32_t uncompressed_size,
								  std::int32_t compressed_size){
	std::string out;
	out.reserve(32);
	std::int16_t last=0;
	append_i32_field(out,last,1,0); // PageType::DATA_PAGE
	append_i32_field(out,last,2,uncompressed_size);
	append_i32_field(out,last,3,compressed_size);
	append_field_header(out,last,5,CompactStruct);
	append_data_page_header(out,num_values);
	append_stop(out);
	return out;
}

std::string make_file_metadata(const std::vector<RowGroupMetadata>&row_groups,
							   std::int64_t num_rows,bool use_zstd){
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
		append_row_group(out,row_group,use_zstd);
	}

	append_binary_field(out,last,6,"calcprimelist-cpp");
	append_stop(out);
	return out;
}

} // namespace calcprime::parquet
