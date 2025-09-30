
#include <gtest/gtest.h>
#include <filesystem>

#include <lib90b/iid.h>

using namespace lib90b;

inline std::filesystem::path normal_filepath = "../../../tests/vectors/normal.bin";

TEST(IidTestSuite, NormalBinEntropy) {
    // Run the IID test suite on the file
    IidResult results = iidTestSuite(normal_filepath);

    // Reference values from the original CLI output
    const double expected_H_bitstring = 0.99631546080565159;
    const double expected_H_original  = 5.6221552772047749;
    const double expected_min_entropy = 5.6221552772047749;

    // Compare results
    EXPECT_NEAR(results.H_bitstring, expected_H_bitstring, 1e-12);
    EXPECT_NEAR(results.H_original, expected_H_original, 1e-12);
    EXPECT_NEAR(results.min_entropy, expected_min_entropy, 1e-12);

    // Check that the key statistical tests passed
    EXPECT_TRUE(results.chi_square.passed) << "Chi-square test failed.";
    EXPECT_TRUE(results.lrs.passed) << "Longest repeated substring test failed.";
    EXPECT_TRUE(results.permutation.passed) << "Permutation test failed.";
}
