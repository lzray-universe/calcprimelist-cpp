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

add_library(calcprime STATIC ${CALCPRIME_SOURCES})
calcprime_configure_library(calcprime)

add_library(calcprime_static STATIC ${CALCPRIME_SOURCES})
calcprime_configure_library(calcprime_static)
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

enable_testing()

add_test(NAME prime_sieve_time_100k
    COMMAND $<TARGET_FILE:prime-sieve> --to 100000 --count --time)
set_tests_properties(prime_sieve_time_100k
    PROPERTIES PASS_REGULAR_EXPRESSION "Elapsed: [0-9]+ us")

add_test(NAME prime_sieve_scientific_to_100k
    COMMAND $<TARGET_FILE:prime-sieve> --to 1e5 --count --time)
set_tests_properties(prime_sieve_scientific_to_100k
    PROPERTIES PASS_REGULAR_EXPRESSION "9592;Elapsed: [0-9]+ us")

