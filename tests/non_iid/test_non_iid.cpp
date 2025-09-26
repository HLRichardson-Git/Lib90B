
#include <gtest/gtest.h>
#include <filesystem>

#include <lib90b/non_iid.h>

using namespace lib90b;

inline std::filesystem::path normal_filepath = "../../../tests/vectors/normal.bin";

TEST(NonIidTestSuite, NormalBinEntropy) {
    // Run the non-iid test suite on the file
    NonIidResult results = nonIidTestSuite(normal_filepath);

    // Reference values from your test vector
    const double expected_H_bitstring = 0.51251197288878858;
    const double expected_H_original = 5.529117785448844;
    const double expected_H_min_entropy = 4.1000957831103086;

    // Compare
    ASSERT_TRUE(results.H_bitstring.has_value());
    EXPECT_NEAR(results.H_bitstring.value(), expected_H_bitstring, 1e-12);

    ASSERT_TRUE(results.H_original.has_value());
    EXPECT_NEAR(results.H_original.value(), expected_H_original, 1e-12);

    ASSERT_TRUE(results.min_entropy.has_value());
    EXPECT_NEAR(results.min_entropy.value(), expected_H_min_entropy, 1e-12);
}
