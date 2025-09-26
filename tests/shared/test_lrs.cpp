
#include <gtest/gtest.h>
#include <filesystem>

#include <lib90b/lrs.h>

using namespace lib90b;

inline std::filesystem::path normal_filepath = "../../../tests/vectors/normal.bin";

TEST(LrsTest, NormalBinLiteral) {

    LrsResult literal_result = lrsTest(normal_filepath);

    // Reference values taken from normal.res
    const int expected_literal_u = 3;
    const int expected_literal_v = 6;
    const double expected_literal_p_hat = 0.014222809025629462;
    const double expected_literal_p_u = 0.01452780869506292;
    const double expected_literal_min_entropy = 6.1050390795897833;

    // Compare
    EXPECT_EQ(literal_result.u, expected_literal_u);
    EXPECT_EQ(literal_result.v, expected_literal_v);
    EXPECT_DOUBLE_EQ(literal_result.p_hat, expected_literal_p_hat);
    EXPECT_DOUBLE_EQ(literal_result.p_u, expected_literal_p_u);
    ASSERT_TRUE(literal_result.h_original.has_value());  // make sure it exists
    EXPECT_DOUBLE_EQ(literal_result.h_original.value(), expected_literal_min_entropy);
}

TEST(LrsTest, NormalBinBitstring) {

    LrsResult bitstring_result = lrsTest(normal_filepath, SymbolMode::Bitstring);

    // Reference values taken from normal.res
    const int expected_bitstring_u = 24;
    const int expected_bitstring_v = 52;
    const double expected_bitstring_p_hat = 0.56270208363310659;
    const double expected_bitstring_p_u = 0.56315383562784926;
    const double expected_bitstring_min_entropy = 0.82839902057212444;

    // Compare
    EXPECT_EQ(bitstring_result.u, expected_bitstring_u);
    EXPECT_EQ(bitstring_result.v, expected_bitstring_v);
    EXPECT_DOUBLE_EQ(bitstring_result.p_hat, expected_bitstring_p_hat);
    EXPECT_DOUBLE_EQ(bitstring_result.p_u, expected_bitstring_p_u);
    ASSERT_TRUE(bitstring_result.h_bitstring.has_value());  // make sure it exists
    EXPECT_DOUBLE_EQ(bitstring_result.h_bitstring.value(), expected_bitstring_min_entropy);
}