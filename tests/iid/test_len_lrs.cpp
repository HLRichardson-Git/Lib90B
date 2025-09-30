
#include <gtest/gtest.h>
#include <filesystem>

#include <lib90b/lrs.h>

using namespace lib90b;

inline std::filesystem::path normal_filepath = "../../../tests/vectors/normal.bin";

TEST(LenLrsTest, NormalBinFile) {
    // Run the test
    lenLrsResult result = lenLrsTest(normal_filepath);
    
    // Expected values from the reference output
    const double expected_p_col = 0.01408999932799633347377;
    const int expected_W = 6;
    const double expected_prob_x_ge_1 = 0.9800053431248045225752;
    const bool expected_passed = true;
    
    // Compare
    EXPECT_DOUBLE_EQ(result.p_col, expected_p_col);
    EXPECT_EQ(result.W, expected_W);
    EXPECT_DOUBLE_EQ(result.prob_x_ge_1, expected_prob_x_ge_1);
    EXPECT_EQ(result.passed, expected_passed);
}