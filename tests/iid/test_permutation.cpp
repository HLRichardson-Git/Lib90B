#include <gtest/gtest.h>
#include <filesystem>

#include <lib90b/permutation.h>

using namespace lib90b;

inline std::filesystem::path normal_perm_filepath = "../../../tests/vectors/normal.bin";

TEST(PermutationTest, NormalBinLiteral) {
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
        14159.0L,                  // periodicity_1
        14106.0L,                  // periodicity_2
        14352.0L,                  // periodicity_8
        14098.0L,                  // periodicity_16
        14240.0L,                  // periodicity_32
        16000707087.0L,            // covariance_1
        16001117960.0L,            // covariance_2
        16000663025.0L,            // covariance_8
        16000590634.0L,            // covariance_16
        16000314179.0L,            // covariance_32
        852701.0L                  // compression
    };

    const std::vector<PermutationSubtestResult*> subtests = {
        &result.excursion,
        &result.numDirectionalRuns,
        &result.lenDirectionalRuns,
        &result.numIncreasesDecreases,
        &result.numRunsMedian,
        &result.lenRunsMedian,
        &result.avgCollision,
        &result.maxCollision,
        &result.periodicity_1,
        &result.periodicity_2,
        &result.periodicity_8,
        &result.periodicity_16,
        &result.periodicity_32,
        &result.covariance_1,
        &result.covariance_2,
        &result.covariance_8,
        &result.covariance_16,
        &result.covariance_32,
        &result.compression
    };

    // Reference counts C[i][0..2]
    const std::vector<std::array<int,3>> expected_counts = {
        {6,0,116}, {21,0,6}, {1,5,12}, {6,0,157}, {6,0,10}, {2,4,5}, {15,0,6}, {10,6,0}, 
        {6,0,9}, {6,0,12}, {6,0,396}, {9,0,6}, {6,0,117}, {6,0,17}, {6,0,24}, {6,0,9}, 
        {6,0,20}, {6,0,12}, {7,0,6}
    };

    ASSERT_EQ(subtests.size(), expected_stats.size());
    ASSERT_EQ(subtests.size(), expected_counts.size());

    for (size_t i = 0; i < subtests.size(); i++) {
        // Check initial statistic
        EXPECT_DOUBLE_EQ(subtests[i]->initialStat, expected_stats[i]) << "Mismatch in initialStat for subtest " << i;

        // Check counts
        EXPECT_EQ(subtests[i]->counts[0], expected_counts[i][0]) << "Mismatch C[i][0] for subtest " << i;
        EXPECT_EQ(subtests[i]->counts[1], expected_counts[i][1]) << "Mismatch C[i][1] for subtest " << i;
        EXPECT_EQ(subtests[i]->counts[2], expected_counts[i][2]) << "Mismatch C[i][2] for subtest " << i;

        // Check that each individual subtest passed
        EXPECT_TRUE(subtests[i]->passed) << "Subtest " << i << " did not pass";
    }

    // Check overall passed flag
    EXPECT_TRUE(result.passed);
}
