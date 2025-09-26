
#include <gtest/gtest.h>
#include <filesystem>

#include <lib90b/markov.h>

using namespace lib90b;

inline std::filesystem::path normal_filepath = "../../../tests/vectors/normal.bin";

TEST(MarkovTest, NormalBinBitstring) {

    MarkovResult bitstring_result = markovTest(normal_filepath);

    // Reference values taken from normal.res
    const double expected_P0 = 0.49917675;
    const double expected_P1 = 0.50082325000000005;
    const double expected_P00 = 0.4969995597751698;
    const double expected_P01 = 0.50300044022483026;
    const double expected_P10 = 0.50134690765327583;
    const double expected_P11 = 0.49865309234672417;
    const double expected_p_hat_max = 5.097359688619317e-39;
    const double expected_entropy = 0.9937925432957424;

    // Compare
    EXPECT_DOUBLE_EQ(bitstring_result.p_0, expected_P0);
    EXPECT_DOUBLE_EQ(bitstring_result.p_1, expected_P1);
    EXPECT_DOUBLE_EQ(bitstring_result.p_00, expected_P00);
    EXPECT_DOUBLE_EQ(bitstring_result.p_01, expected_P01);
    EXPECT_DOUBLE_EQ(bitstring_result.p_10, expected_P10);
    EXPECT_DOUBLE_EQ(bitstring_result.p_11, expected_P11);
    EXPECT_DOUBLE_EQ(bitstring_result.p_hat_max, expected_p_hat_max);
    ASSERT_TRUE(bitstring_result.h_bitstring.has_value());  // make sure it exists
    EXPECT_DOUBLE_EQ(bitstring_result.h_bitstring.value(), expected_entropy);
}
