
#include <gtest/gtest.h>
#include <filesystem>

#include <lib90b/compression.h>

using namespace lib90b;

inline std::filesystem::path normal_filepath = "../../../tests/vectors/normal.bin";

TEST(CompressionTest, NormalBinBitstring) {

    CompressionResult bitstring_result = compressionTest(normal_filepath);

    // Reference values taken from normal.res
    const double expected_bitstring_x_bar = 5.0306901312440173;
    const double expected_bitstring_sigma_hat = 1.0488632759418641;
    const double expected_bitstring_x_bar_prime = 5.0283495184921581;
    const bool expected_bitstring_found_p = true;
    const double expected_bitstring_p = 0.11866185634166618;
    const double expected_bitstring_entropy = 0.51251197288878858;

    // Compare
    EXPECT_DOUBLE_EQ(bitstring_result.x_bar, expected_bitstring_x_bar);
    EXPECT_DOUBLE_EQ(bitstring_result.sigma_hat, expected_bitstring_sigma_hat);
    EXPECT_DOUBLE_EQ(bitstring_result.x_bar_prime, expected_bitstring_x_bar_prime);
    EXPECT_EQ(bitstring_result.found_p, expected_bitstring_found_p);
    EXPECT_NEAR(bitstring_result.p, expected_bitstring_p, 1e-14);
    EXPECT_TRUE(bitstring_result.h_bitstring.has_value());  // make sure it exists
    EXPECT_NEAR(bitstring_result.h_bitstring.value(), expected_bitstring_entropy, 1e-14);
}