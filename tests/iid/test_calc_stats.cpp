
#include <gtest/gtest.h>
#include <filesystem>

#include <lib90b/calc_stats.h>

using namespace lib90b;

inline std::filesystem::path normal_filepath = "../../../tests/vectors/normal.bin";

TEST(CalcStatsTest, NormalBinLiteral) {
    CalcStatsResult result = calcStats(normal_filepath);

    // Reference values taken from normal.res
    const double expected_mean   = 126.493376;
    const double expected_median = 90.0;

    // Compare
    EXPECT_EQ(result.mean, expected_mean);
    EXPECT_EQ(result.median, expected_median);
}
