
#include <gtest/gtest.h>
#include <filesystem>

#include <lib90b/multi_mcw.h>

using namespace lib90b;

inline std::filesystem::path normal_filepath = "../../../tests/vectors/normal.bin";

TEST(MultiMcwTest, NormalBinLiteral) {

    MultiMcwResult literal_result = multiMcwTest(normal_filepath);

    // Reference values taken from normal.res
    const uint64_t expected_literal_C = 19310;
    const int expected_literal_r = 4;
    const uint64_t expected_literal_N = 999937;
    const double expected_literal_P_global = 0.019311216606646218;
    const double expected_literal_P_global_prime = 0.019665704493110798;
    const double expected_literal_min_entropy = 5.6681743202743657;

    // Compare
    EXPECT_EQ(literal_result.C, expected_literal_C);
    EXPECT_EQ(literal_result.r, expected_literal_r);
    EXPECT_EQ(literal_result.N, expected_literal_N);
    EXPECT_DOUBLE_EQ(literal_result.p_global, expected_literal_P_global);
    EXPECT_DOUBLE_EQ(literal_result.p_global_prime, expected_literal_P_global_prime);
    ASSERT_TRUE(literal_result.h_original.has_value());  // make sure it exists
    EXPECT_DOUBLE_EQ(literal_result.h_original.value(), expected_literal_min_entropy);
}

TEST(MultiMcwTest, NormalBinBitstring) {

    MultiMcwResult bitstring_result = multiMcwTest(normal_filepath, SymbolMode::Bitstring);

    // Reference values taken from normal.res
    const uint64_t expected_bitstring_C = 3996310;
    const int expected_bitstring_r = 10;
    const uint64_t expected_bitstring_N = 7999937;
    const double expected_bitstring_P_global = 0.49954268389863571;
    const double expected_bitstring_P_global_prime = 0.49999803212150123;
    const double expected_bitstring_min_entropy = 1.0;

    // Compare
    EXPECT_EQ(bitstring_result.C, expected_bitstring_C);
    EXPECT_EQ(bitstring_result.r, expected_bitstring_r);
    EXPECT_EQ(bitstring_result.N, expected_bitstring_N);
    EXPECT_DOUBLE_EQ(bitstring_result.p_global, expected_bitstring_P_global);
    EXPECT_DOUBLE_EQ(bitstring_result.p_global_prime, expected_bitstring_P_global_prime);
    ASSERT_TRUE(bitstring_result.h_bitstring.has_value());  // make sure it exists
    EXPECT_DOUBLE_EQ(bitstring_result.h_bitstring.value(), expected_bitstring_min_entropy);
}
