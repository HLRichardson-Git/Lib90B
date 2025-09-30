
#include <gtest/gtest.h>
#include <filesystem>

#include <lib90b/permutation.h>

using namespace lib90b;

inline std::filesystem::path normal_perm_filepath = "../../../tests/vectors/normal.bin";

TEST(PermutationTest, NormalBinLiteral) {
    // Choose a fixed seed for reproducibility
    uint64_t fixed_seed = 123456789ULL;

    PermutationTestResult result = permutationTest(normal_perm_filepath, fixed_seed);

    // Reference baseline stats from the test vector
    const std::vector<long double> expected_stats = {
        27576.77190399914979935L,  // excursion
        666224.0L,                 // numDirectionalRuns
        9.0L,                      // lenDirectionalRuns
        507540.0L,                 // numIncreasesDecreases
        499507.0L,                 // numRunsMedian
        20.0L,                     // lenRunsMedian
        11.32167198786314266101L,  // avgCollision
        38.0L,                     // maxCollision
        14159.0L,                  // periodicity(1)
        14106.0L,                  // periodicity(2)
        14352.0L,                  // periodicity(8)
        14098.0L,                  // periodicity(16)
        14240.0L,                  // periodicity(32)
        16000707087.0L,            // covariance(1)
        16001117960.0L,            // covariance(2)
        16000663025.0L,            // covariance(8)
        16000590634.0L,            // covariance(16)
        16000314179.0L,            // covariance(32)
        852701.0L                  // compression
    };

    ASSERT_EQ(result.initialStats.size(), expected_stats.size());
    for (size_t i = 0; i < expected_stats.size(); i++) {
        EXPECT_DOUBLE_EQ(result.initialStats[i], expected_stats[i]);
    }

    // Reference counts for C[i][0..2]
    const std::vector<std::array<int,3>> expected_counts = {
        {6,0,116}, {21,0,6}, {1,5,12}, {6,0,157}, {6,0,10}, {2,4,5}, {15,0,6}, {10,6,0}, {6,0,9}, {6,0,12}, 
        {6,0,396}, {9,0,6}, {6,0,117}, {6,0,17}, {6,0,24}, {6,0,9}, {6,0,20}, {6,0,12}, {7,0,6}
    };

    ASSERT_EQ(result.counts.size(), expected_counts.size());
    for (size_t i = 0; i < expected_counts.size(); i++) {
        EXPECT_EQ(result.counts[i][0], expected_counts[i][0]) << "Mismatch C[i][0] at stat " << i;
        EXPECT_EQ(result.counts[i][1], expected_counts[i][1]) << "Mismatch C[i][1] at stat " << i;
        EXPECT_EQ(result.counts[i][2], expected_counts[i][2]) << "Mismatch C[i][2] at stat " << i;
    }

    // Passed flag
    EXPECT_TRUE(result.passed);
}
