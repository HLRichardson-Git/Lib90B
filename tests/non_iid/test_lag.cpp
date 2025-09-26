
#include <gtest/gtest.h>
#include <filesystem>

#include <lib90b/lag.h>

using namespace lib90b;

inline std::filesystem::path normal_filepath = "../../../tests/vectors/normal.bin";

TEST(LagTest, NormalBinLiteral) {
    LagResult literal_result = lagTest(normal_filepath);

    // Reference values taken from normal.res
    const uint64_t expected_C = 14211;
    const uint64_t expected_r = 4;
    const uint64_t expected_N = 999999;
    const double expected_P_global = 0.014211014211014211;
    const double expected_P_global_prime = 0.014515889364164508;
    const double expected_entropy = 6.1062232235998328;

    // Compare
    EXPECT_EQ(literal_result.C, expected_C);
    EXPECT_EQ(literal_result.r, expected_r);
    EXPECT_EQ(literal_result.N, expected_N);
    EXPECT_DOUBLE_EQ(literal_result.p_global, expected_P_global);
    EXPECT_DOUBLE_EQ(literal_result.p_global_prime, expected_P_global_prime);
    ASSERT_TRUE(literal_result.h_original.has_value());  // make sure it exists
    EXPECT_DOUBLE_EQ(literal_result.h_original.value(), expected_entropy);
}

TEST(LagTest, NormalBinBitstring) {
    LagResult bitstring_result = lagTest(normal_filepath, SymbolMode::Bitstring);

    // Reference values taken from normal.res
    const uint64_t expected_C = 4002718;
    const uint64_t expected_r = 27;
    const uint64_t expected_N = 7999999;
    const double expected_P_global = 0.50033981254247661;
    const double expected_P_global_prime = 0.50079515908616445;
    const double expected_entropy = 0.99770747829602213;

    // Compare
    EXPECT_EQ(bitstring_result.C, expected_C);
    EXPECT_EQ(bitstring_result.r, expected_r);
    EXPECT_EQ(bitstring_result.N, expected_N);
    EXPECT_DOUBLE_EQ(bitstring_result.p_global, expected_P_global);
    EXPECT_DOUBLE_EQ(bitstring_result.p_global_prime, expected_P_global_prime);
    ASSERT_TRUE(bitstring_result.h_bitstring.has_value());  // make sure it exists
    EXPECT_DOUBLE_EQ(bitstring_result.h_bitstring.value(), expected_entropy);
}
