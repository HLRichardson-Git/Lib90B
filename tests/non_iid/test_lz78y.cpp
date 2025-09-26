
#include <gtest/gtest.h>
#include <filesystem>

#include <lib90b/lz78y.h>

using namespace lib90b;

inline std::filesystem::path normal_filepath = "../../../tests/vectors/normal.bin";

TEST(LZ78YTest, NormalBinLiteral) {
    Lz78yResult literal_result = lz78yTest(normal_filepath);

    // Reference values taken from normal.res
    const uint64_t expected_C = 19163;
    const uint64_t expected_r = 4;
    const uint64_t expected_N = 999983;
    const double expected_P_global = 0.0191633257765382;
    const double expected_P_global_prime = 0.019516472171868238;
    const double expected_entropy = 5.6791638971261458;

    // Compare
    EXPECT_EQ(literal_result.C, expected_C);
    EXPECT_EQ(literal_result.r, expected_r);
    EXPECT_EQ(literal_result.N, expected_N);
    EXPECT_DOUBLE_EQ(literal_result.p_global, expected_P_global);
    EXPECT_DOUBLE_EQ(literal_result.p_global_prime, expected_P_global_prime);
    ASSERT_TRUE(literal_result.h_original.has_value());  // make sure it exists
    EXPECT_DOUBLE_EQ(literal_result.h_original.value(), expected_entropy);
}

TEST(LZ78YTest, NormalBinBitstring) {
    Lz78yResult bitstring_result = lz78yTest(normal_filepath, SymbolMode::Bitstring);

    // Reference values taken from normal.res
    const uint64_t expected_C = 4017307;
    const uint64_t expected_r = 20;
    const uint64_t expected_N = 7999983;
    const double expected_P_global = 0.50216444209943945;
    const double expected_P_global_prime = 0.50261978493718584;
    const double expected_entropy = 0.99246063284315811;

    // Compare
    EXPECT_EQ(bitstring_result.C, expected_C);
    EXPECT_EQ(bitstring_result.r, expected_r);
    EXPECT_EQ(bitstring_result.N, expected_N);
    EXPECT_DOUBLE_EQ(bitstring_result.p_global, expected_P_global);
    EXPECT_DOUBLE_EQ(bitstring_result.p_global_prime, expected_P_global_prime);
    ASSERT_TRUE(bitstring_result.h_bitstring.has_value());  // make sure it exists
    EXPECT_DOUBLE_EQ(bitstring_result.h_bitstring.value(), expected_entropy);
}
