
#include <gtest/gtest.h>
#include <filesystem>
#include <lib90b/chi_square.h>

using namespace lib90b;

inline std::filesystem::path normal_chi_filepath = "../../../tests/vectors/normal.bin";

TEST(ChiSquareTest, NormalBinLiteral) {
    ChiSquareResult result = chiSquareTest(normal_chi_filepath);

    // Reference values taken from your test vector
    const double expected_independence_score = 11031.596498190238;
    const int expected_independence_df = 11036;
    const double expected_independence_pvalue = 0.51003488778598838;
    const double expected_goodness_score = 1304.8386764863249;
    const int expected_goodness_df = 1314;
    const double expected_goodness_pvalue = 0.56596600311521694;
    const bool expected_passed = true;

    // Compare values
    EXPECT_DOUBLE_EQ(result.independence_score, expected_independence_score);
    EXPECT_EQ(result.independence_df, expected_independence_df);
    EXPECT_DOUBLE_EQ(result.independence_pvalue, expected_independence_pvalue);
    EXPECT_DOUBLE_EQ(result.goodness_of_fit_score, expected_goodness_score);
    EXPECT_EQ(result.goodness_of_fit_df, expected_goodness_df);
    EXPECT_NEAR(result.goodness_of_fit_pvalue, expected_goodness_pvalue, 1e-9);
    EXPECT_EQ(result.passed, expected_passed);
}
