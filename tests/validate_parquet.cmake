if(NOT DEFINED CALCPRIME_EXE OR NOT DEFINED OUTPUT_FILE)
    message(FATAL_ERROR "CALCPRIME_EXE and OUTPUT_FILE are required")
endif()

set(extra_args)
if(USE_ZSTD)
    list(APPEND extra_args --zstd)
endif()
if(USE_DELTA)
    list(APPEND extra_args --parquet-encoding delta)
    if(DEFINED DELTA_BLOCK_VALUES)
        list(APPEND extra_args --parquet-delta-block-values
             "${DELTA_BLOCK_VALUES}")
    endif()
endif()

execute_process(
    COMMAND "${CALCPRIME_EXE}" --from 1 --to 100 --print
            --out "${OUTPUT_FILE}" --out-format parquet ${extra_args}
    RESULT_VARIABLE result
    ERROR_VARIABLE error_output)
if(NOT result EQUAL 0)
    message(FATAL_ERROR "Parquet export failed: ${error_output}")
endif()

file(SIZE "${OUTPUT_FILE}" file_size)
if(file_size LESS 12)
    message(FATAL_ERROR "Parquet output is too small: ${file_size} bytes")
endif()

file(READ "${OUTPUT_FILE}" leading_magic OFFSET 0 LIMIT 4 HEX)
math(EXPR trailing_offset "${file_size} - 4")
file(READ "${OUTPUT_FILE}" trailing_magic OFFSET ${trailing_offset} LIMIT 4 HEX)
if(NOT leading_magic STREQUAL "50415231" OR
   NOT trailing_magic STREQUAL "50415231")
    message(FATAL_ERROR "Parquet output is missing PAR1 magic bytes")
endif()
