
#include <gtest/gtest.h>
#include <filesystem>

#include <lib90b/t_tuple.h>  // include the new t-tuple API

using namespace lib90b;

inline std::filesystem::path normal_filepath = "../../../tests/vectors/normal.bin";

TEST(TTupleTest, NormalBinLiteral) {

    TtupleResult literal_result = tTupleTest(normal_filepath);

    // Reference values from your vector
    const int expected_literal_t = 2;
    const double expected_literal_p_hat_max = 0.021283807295699071;
    const double expected_literal_p_u = 0.021655573872636183;
    const double expected_literal_min_entropy = 5.529117785448844;

    // Compare
    EXPECT_EQ(literal_result.t, expected_literal_t);
    EXPECT_DOUBLE_EQ(literal_result.p_hat_max, expected_literal_p_hat_max);
    EXPECT_DOUBLE_EQ(literal_result.p_u, expected_literal_p_u);
    ASSERT_TRUE(literal_result.h_original.has_value());
    EXPECT_DOUBLE_EQ(literal_result.h_original.value(), expected_literal_min_entropy);
}

TEST(TTupleTest, NormalBinBitstring) {

    TtupleResult bitstring_result = tTupleTest(normal_filepath, SymbolMode::Bitstring);

    // Reference values from your vector
    const int expected_bitstring_t = 23;
    const double expected_bitstring_p_hat_max = 0.58478886944022623;
    const double expected_bitstring_p_u = 0.58523762119077116;
    const double expected_bitstring_min_entropy = 0.77290558077529137;

    // Compare
    EXPECT_EQ(bitstring_result.t, expected_bitstring_t);
    EXPECT_DOUBLE_EQ(bitstring_result.p_hat_max, expected_bitstring_p_hat_max);
    EXPECT_DOUBLE_EQ(bitstring_result.p_u, expected_bitstring_p_u);
    ASSERT_TRUE(bitstring_result.h_bitstring.has_value());
    EXPECT_DOUBLE_EQ(bitstring_result.h_bitstring.value(), expected_bitstring_min_entropy);
}
