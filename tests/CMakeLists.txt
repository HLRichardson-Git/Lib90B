cmake_minimum_required(VERSION 3.16.3)

include(GoogleTest)

set (NON_IID_TEST_SOURCES
    non_iid/test_collision.cpp    
    non_iid/test_compression.cpp
    non_iid/test_lag.cpp
    non_iid/test_lz78y.cpp
    non_iid/test_markov.cpp
    non_iid/test_multi_mcw.cpp
    non_iid/test_multi_mmc.cpp
    non_iid/test_t_tuple.cpp
    non_iid/test_non_iid.cpp
)

set (SHARED_TEST_SOURCES
    shared/test_most_common.cpp
    shared/test_lrs.cpp
)

add_executable (LIB_90B_TESTS
    ${NON_IID_TEST_SOURCES}
    ${SHARED_TEST_SOURCES}
)

target_include_directories(LIB_90B_TESTS
    PRIVATE
        ${PROJECT_SOURCE_DIR}/src
        ${PROJECT_SOURCE_DIR}/util
        ${PROJECT_SOURCE_DIR}/third_party
    PUBLIC
        ${PROJECT_SOURCE_DIR}/include
)

target_link_libraries (LIB_90B_TESTS 
    PRIVATE
        GTest::gtest_main
        Lib90B
)

add_test (
    NAME LIB_90B_TESTS
    COMMAND LIB_90B_TESTS
)

include(GoogleTest)
gtest_discover_tests(LIB_90B_TESTS)