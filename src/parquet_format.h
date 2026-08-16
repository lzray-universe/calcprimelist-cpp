#pragma once

#include<cstdint>
#include<cstddef>
#include<string>
#include<vector>

namespace calcprime::parquet{

enum class ValueEncoding : std::int32_t{
	Plain=0,
	DeltaBinaryPacked=5,
};

struct RowGroupMetadata{
	std::int64_t data_page_offset=0;
	std::int64_t num_values=0;
	std::int64_t total_uncompressed_size=0;
	std::int64_t total_compressed_size=0;
};

std::string make_data_page_header(std::int32_t num_values,
								  std::int32_t uncompressed_size,
								  std::int32_t compressed_size,
								  ValueEncoding encoding);

std::string make_file_metadata(const std::vector<RowGroupMetadata>&row_groups,
							   std::int64_t num_rows,bool use_zstd,
							   ValueEncoding encoding);

std::string encode_delta_binary_packed(
	const std::uint64_t*values,std::size_t count,
	std::size_t block_value_count);

} // namespace calcprime::parquet
