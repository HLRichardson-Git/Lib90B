
#include <gtest/gtest.h>
#include <filesystem>

#include <lib90b/most_common.h>

using namespace lib90b;

inline std::filesystem::path normal_filepath = "../../../tests/vectors/normal.bin";

TEST(MostCommonTest, NormalBinLiteral) {

    MostCommonResult literal_result = mostCommonTest(normal_filepath);

    // Reference values taken from normal.res
    const uint64_t expected_literal_mode_count = 19943;
    const double expected_literal_p_hat = 0.019942999999999999;
    const double expected_literal_p_u = 0.02030311251014405;
    const double expected_literal_entropy = 5.6221552772047749;

    // Compare
    EXPECT_EQ(literal_result.mode_count, expected_literal_mode_count);
    EXPECT_DOUBLE_EQ(literal_result.p_hat, expected_literal_p_hat);
    EXPECT_DOUBLE_EQ(literal_result.p_u, expected_literal_p_u);
    ASSERT_TRUE(literal_result.h_original.has_value());  // make sure it exists
    EXPECT_DOUBLE_EQ(literal_result.h_original.value(), expected_literal_entropy);
}

TEST(MostCommonTest, NormalBinBitstring) {

    MostCommonResult bitstring_result = mostCommonTest(normal_filepath, SymbolMode::Bitstring);

    // Reference values taken from normal.res
    const uint64_t expected_bitstring_mode_count = 4006586;
    const double expected_bitstring_p_hat = 0.50082325000000005;
    const double expected_bitstring_p_u = 0.5012785960031747;
    const double expected_bitstring_entropy = 0.99631546080565159;

    // Compare
    EXPECT_EQ(bitstring_result.mode_count, expected_bitstring_mode_count);
    EXPECT_DOUBLE_EQ(bitstring_result.p_hat, expected_bitstring_p_hat);
    EXPECT_DOUBLE_EQ(bitstring_result.p_u, expected_bitstring_p_u);
    ASSERT_TRUE(bitstring_result.h_bitstring.has_value());  // make sure it exists
    EXPECT_DOUBLE_EQ(bitstring_result.h_bitstring.value(), expected_bitstring_entropy);
}