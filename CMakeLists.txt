cmake_minimum_required(VERSION 3.16.3)
project(Lib90B LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 23)

option(Lib90B_BUILD_TESTS "Build unit tests" ON)

include(FetchContent)

set(BUILD_SHARED_LIBS OFF CACHE BOOL "Build static libraries" FORCE)
set(BUILD_DIVSUFSORT64 ON CACHE BOOL "Enable 64-bit suffix array support" FORCE)

FetchContent_Declare(
    divsufsort
    URL https://github.com/y-256/libdivsufsort/archive/refs/tags/2.0.1.zip
)

FetchContent_MakeAvailable(divsufsort)

# Fetch bzip2
FetchContent_Declare(
    bzip2
    URL https://github.com/libarchive/bzip2/archive/refs/tags/bzip2-1.0.8.zip
)

FetchContent_MakeAvailable(bzip2)

include_directories(
    ${PROJECT_SOURCE_DIR}/include
)

set (NON_IID_SOURCES
    src/non_iid/collision.cpp
    src/non_iid/compression.cpp
    src/non_iid/lag.cpp
    src/non_iid/lz78y.cpp
    src/non_iid/markov.cpp
    src/non_iid/multi_mcw.cpp
    src/non_iid/multi_mmc.cpp
    src/non_iid/t_tuple.cpp
    src/non_iid/non_iid.cpp
)

set (IID_SOURCES
    src/iid/calc_stats.cpp
    src/iid/chi_square.cpp
    src/iid/permutation.cpp
    src/iid/iid.cpp
)

set (SHARED_SOURCES
    src/shared/most_common.cpp
    src/shared/lrs/lrs.cpp
    src/shared/lrs/lrs_algs.cpp
)

set (UTIL_SOURCES
    util/TestRunUtils.h
    util/utils.cpp
)

set (GESTALT_SOURCES
    third_party/gestalt/sha2_256.cpp
)

add_library (${PROJECT_NAME} STATIC 
    ${NON_IID_SOURCES}
    ${IID_SOURCES}
    ${SHARED_SOURCES}
    ${UTIL_SOURCES}
    ${GESTALT_SOURCES}
)

add_library(bzip2_lib STATIC
    ${bzip2_SOURCE_DIR}/blocksort.c
    ${bzip2_SOURCE_DIR}/bzlib.c
    ${bzip2_SOURCE_DIR}/compress.c
    ${bzip2_SOURCE_DIR}/crctable.c
    ${bzip2_SOURCE_DIR}/decompress.c
    ${bzip2_SOURCE_DIR}/huffman.c
    ${bzip2_SOURCE_DIR}/randtable.c
)

target_include_directories(bzip2_lib PUBLIC ${bzip2_SOURCE_DIR})
target_compile_definitions(bzip2_lib PUBLIC BZ_HAVE_BZBUFFTOBUFF)

target_include_directories(${PROJECT_NAME}
    PRIVATE
        ${PROJECT_SOURCE_DIR}/src
        ${PROJECT_SOURCE_DIR}/util
        ${PROJECT_SOURCE_DIR}/third_party
    PUBLIC
        ${PROJECT_SOURCE_DIR}/include
        ${divsufsort_BINARY_DIR}/include
        ${bzip2_SOURCE_DIR}
)

target_link_libraries(${PROJECT_NAME}
    PUBLIC
        divsufsort
        divsufsort64
        bzip2_lib
)

# Enable testing and include GoogleTest if Lib90B_BUILD_TESTS is ON
if(Lib90B_BUILD_TESTS)
    enable_testing()

    add_subdirectory(third_party/googletest)
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(googletest)

    target_link_libraries(${PROJECT_NAME} PRIVATE gtest_main)
    add_subdirectory(tests)
endif()
