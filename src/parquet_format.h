#pragma once

#include<cstdint>
#include<string>
#include<vector>

namespace calcprime::parquet{

struct RowGroupMetadata{
	std::int64_t data_page_offset=0;
	std::int64_t num_values=0;
	std::int64_t total_uncompressed_size=0;
	std::int64_t total_compressed_size=0;
};

std::string make_data_page_header(std::int32_t num_values,
								  std::int32_t uncompressed_size,
								  std::int32_t compressed_size);

std::string make_file_metadata(const std::vector<RowGroupMetadata>&row_groups,
							   std::int64_t num_rows,bool use_zstd);

} // namespace calcprime::parquet
