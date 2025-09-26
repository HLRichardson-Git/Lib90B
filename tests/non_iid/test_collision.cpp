
#include <gtest/gtest.h>
#include <filesystem>

#include <lib90b/collision.h>

using namespace lib90b;

inline std::filesystem::path normal_filepath = "../../../tests/vectors/normal.bin";

TEST(CollisionTest, NormalBinBitstring) {

    CollisionResult bitstring_result = collisionTest(normal_filepath);

    // Reference values taken from normal.res
    const uint64_t expected_bitstring_v = 3170586;
    const uint64_t expected_bitstring_sum_ti = 7999999;
    const double expected_bitstring_x_bar = 2.523192558094939;
    const double expected_bitstring_sigma_hat = 0.49946189437149541;
    const double expected_bitstring_x_bar_prime = 2.5224700384317478;
    const bool expected_bitstring_used_lower_bound = false;
    const double expected_bitstring_p = 0.5;
    const double expected_bitstring_entropy = 1.0;

    // Compare
    EXPECT_EQ(bitstring_result.v, expected_bitstring_v);
    EXPECT_EQ(bitstring_result.sum_ti, expected_bitstring_sum_ti);
    EXPECT_DOUBLE_EQ(bitstring_result.x_bar, expected_bitstring_x_bar);
    EXPECT_DOUBLE_EQ(bitstring_result.sigma_hat, expected_bitstring_sigma_hat);
    EXPECT_DOUBLE_EQ(bitstring_result.x_bar_prime, expected_bitstring_x_bar_prime);
    EXPECT_EQ(bitstring_result.used_lower_bound, expected_bitstring_used_lower_bound);
    EXPECT_DOUBLE_EQ(bitstring_result.p, expected_bitstring_p);
    ASSERT_TRUE(bitstring_result.h_bitstring.has_value());  // make sure it exists
    EXPECT_DOUBLE_EQ(bitstring_result.h_bitstring.value(), expected_bitstring_entropy);
}