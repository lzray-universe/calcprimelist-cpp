cmake_minimum_required(VERSION 3.16)
project(calcprimelist LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CALCPRIME_SOURCES
    src/cpu_info.cpp
    src/wheel.cpp
    src/base_sieve.cpp
    src/bucket.cpp
    src/marker.cpp
    src/popcnt.cpp
    src/prime_count.cpp
    src/segmenter.cpp
    src/writer.cpp
)

option(CALCPRIME_WITH_ZSTD "Enable zstd compression if available" ON)

set(CALCPRIME_HAS_ZSTD FALSE)
set(CALCPRIME_ZSTD_TARGET "")
set(CALCPRIME_ZSTD_INCLUDE_DIR "")
set(CALCPRIME_ZSTD_LIBRARY "")

if(CALCPRIME_WITH_ZSTD)
    find_package(zstd CONFIG QUIET)
    if(TARGET zstd::libzstd_static)
        set(CALCPRIME_ZSTD_TARGET zstd::libzstd_static)
    elseif(TARGET zstd::libzstd_shared)
        set(CALCPRIME_ZSTD_TARGET zstd::libzstd_shared)
    elseif(TARGET zstd::libzstd)
        set(CALCPRIME_ZSTD_TARGET zstd::libzstd)
    endif()

    if(NOT CALCPRIME_ZSTD_TARGET)
        find_path(ZSTD_INCLUDE_DIR zstd.h)
        find_library(ZSTD_LIBRARY NAMES zstd libzstd_static zstd_static libzstd)
        if(ZSTD_INCLUDE_DIR AND ZSTD_LIBRARY)
            set(CALCPRIME_ZSTD_INCLUDE_DIR "${ZSTD_INCLUDE_DIR}")
            set(CALCPRIME_ZSTD_LIBRARY "${ZSTD_LIBRARY}")
        endif()
    endif()

    if(CALCPRIME_ZSTD_TARGET OR
       (CALCPRIME_ZSTD_INCLUDE_DIR AND CALCPRIME_ZSTD_LIBRARY))
        set(CALCPRIME_HAS_ZSTD TRUE)
        message(STATUS "zstd support enabled")
    else()
        message(STATUS "zstd not found, building without compression")
    endif()
else()
    message(STATUS "CALCPRIME_WITH_ZSTD=OFF, building without compression")
endif()

function(calcprime_configure_library target)
    target_include_directories(${target} PUBLIC include)
    target_compile_definitions(${target} PRIVATE _USE_MATH_DEFINES)
    set_target_properties(${target} PROPERTIES POSITION_INDEPENDENT_CODE ON)

    if(MSVC)
        target_compile_options(${target} PRIVATE /arch:AVX2 /W4)
    else()
        target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic)
        if(CMAKE_SYSTEM_PROCESSOR MATCHES "(x86_64|AMD64|i.86)")
            target_compile_options(${target} PRIVATE -mavx2)
        endif()
    endif()
endfunction()

function(calcprime_enable_zstd target)
    if(CALCPRIME_HAS_ZSTD)
        target_compile_definitions(${target} PUBLIC CALCPRIME_HAS_ZSTD)
        if(CALCPRIME_ZSTD_TARGET)
            target_link_libraries(${target} PUBLIC ${CALCPRIME_ZSTD_TARGET})
        else()
            target_include_directories(${target} PUBLIC ${CALCPRIME_ZSTD_INCLUDE_DIR})
            target_link_libraries(${target} PUBLIC ${CALCPRIME_ZSTD_LIBRARY})
        endif()
    endif()
endfunction()

add_library(calcprime STATIC ${CALCPRIME_SOURCES})
calcprime_configure_library(calcprime)
calcprime_enable_zstd(calcprime)

add_library(calcprime_static STATIC ${CALCPRIME_SOURCES})
calcprime_configure_library(calcprime_static)
calcprime_enable_zstd(calcprime_static)
if(MSVC)
    set_property(TARGET calcprime_static PROPERTY MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
endif()

add_executable(prime-sieve src/main.cpp)

target_link_libraries(prime-sieve PRIVATE calcprime)

add_executable(prime-sieve-static src/main.cpp)
target_link_libraries(prime-sieve-static PRIVATE calcprime_static)
if(MSVC)
    set_property(TARGET prime-sieve-static PROPERTY MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
    target_compile_definitions(prime-sieve-static PRIVATE _CRT_SECURE_NO_WARNINGS)
endif()

add_library(calcprime-cli SHARED src/main.cpp src/api.cpp)
target_link_libraries(calcprime-cli PRIVATE calcprime)
target_compile_definitions(calcprime-cli PRIVATE CALCPRIME_DLL_BUILD CALCPRIME_DLL_EXPORT)
target_include_directories(calcprime-cli PUBLIC include)
calcprime_enable_zstd(calcprime-cli)

enable_testing()

add_test(NAME prime_sieve_time_100k
    COMMAND $<TARGET_FILE:prime-sieve> --to 100000 --count --time)
set_tests_properties(prime_sieve_time_100k
    PROPERTIES PASS_REGULAR_EXPRESSION "Elapsed: [0-9]+ us")

add_test(NAME prime_sieve_scientific_to_100k
    COMMAND $<TARGET_FILE:prime-sieve> --to 1e5 --count --time)
set_tests_properties(prime_sieve_scientific_to_100k
    PROPERTIES PASS_REGULAR_EXPRESSION "9592;Elapsed: [0-9]+ us")
