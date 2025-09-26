
#include <gtest/gtest.h>
#include <filesystem>

#include <lib90b/multi_mmc.h>

using namespace lib90b;

inline std::filesystem::path normal_filepath = "../../../tests/vectors/normal.bin";

TEST(MultiMmcTest, NormalBinLiteral) {

    MultiMmcResult literal_result = multiMmcTest(normal_filepath);

    // Reference values taken from normal.res
    const uint64_t expected_literal_C = 19209;
    const int expected_literal_r = 4;
    const uint64_t expected_literal_N = 999998;
    const double expected_literal_P_global = 0.019209038418076835;
    const double expected_literal_P_global_prime = 0.019562594873373633;
    const double expected_literal_min_entropy = 5.6757584410258906;

    // Compare
    EXPECT_EQ(literal_result.C, expected_literal_C);
    EXPECT_EQ(literal_result.r, expected_literal_r);
    EXPECT_EQ(literal_result.N, expected_literal_N);
    EXPECT_DOUBLE_EQ(literal_result.p_global, expected_literal_P_global);
    EXPECT_DOUBLE_EQ(literal_result.p_global_prime, expected_literal_P_global_prime);
    ASSERT_TRUE(literal_result.h_original.has_value());  // make sure it exists
    EXPECT_DOUBLE_EQ(literal_result.h_original.value(), expected_literal_min_entropy);
}

TEST(MultiMmcTest, NormalBinBitstring) {

    MultiMmcResult bitstring_result = multiMmcTest(normal_filepath, SymbolMode::Bitstring);

    // Reference values taken from normal.res
    const uint64_t expected_bitstring_C = 5001029;
    const int expected_bitstring_r = 41;
    const uint64_t expected_bitstring_N = 7999998;
    const double expected_bitstring_P_global = 0.62512878128219529;
    const double expected_bitstring_P_global_prime = 0.62556963850841341;
    const double expected_bitstring_min_entropy = 0.6767576005226027;

    // Compare
    EXPECT_EQ(bitstring_result.C, expected_bitstring_C);
    EXPECT_EQ(bitstring_result.r, expected_bitstring_r);
    EXPECT_EQ(bitstring_result.N, expected_bitstring_N);
    EXPECT_DOUBLE_EQ(bitstring_result.p_global, expected_bitstring_P_global);
    EXPECT_DOUBLE_EQ(bitstring_result.p_global_prime, expected_bitstring_P_global_prime);
    ASSERT_TRUE(bitstring_result.h_bitstring.has_value());  // make sure it exists
    EXPECT_DOUBLE_EQ(bitstring_result.h_bitstring.value(), expected_bitstring_min_entropy);
}
