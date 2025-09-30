
#include <vector>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <algorithm>
#include <stdexcept>
#include <unordered_set>

#include <lib90b/lrs.h>
#include "shared/lrs/lrs_algs.h"
#include <lib90b/entropy_tests.h>
#include "utils.h"

namespace lib90b {

LrsResult lrsTest(const EntropyInputData& data, SymbolMode mode) {
    const auto& symbols = (mode == SymbolMode::Bitstring) ? data.getBits() : data.symbols;
    InternalLrsResult internal = SAalgs(symbols.data(), symbols.size(),
                                        (mode == SymbolMode::Bitstring) ? 2 : (1 << data.word_size));

    LrsResult result;
    result.u = internal.u;
    result.v = internal.v;
    result.p_hat = internal.p_hat;
    result.p_u = internal.p_u;
    result.lrs_res = internal.lrs_res;

    if (mode == SymbolMode::Original) result.h_original = internal.lrs_res;
    else result.h_bitstring = internal.lrs_res;

    return result;
}


LrsResult lrsTest(const std::filesystem::path& filepath, SymbolMode mode) {
    std::ifstream file(filepath, std::ios::binary);
    if (!file) {
        throw std::runtime_error("lrsTest: failed to open file " + filepath.string());
    }

    std::vector<uint8_t> buffer(
        (std::istreambuf_iterator<char>(file)),
        std::istreambuf_iterator<char>());

    EntropyInputData data{ buffer, 8 };
    return lrsTest(data, mode);
}

lenLrsResult lenLrsTest(const EntropyInputData& data) {
    lenLrsResult result{};
    
    int L = static_cast<int>(data.symbols.size());
    int k = data.alph_size > 0 ? data.alph_size : 256;
    
    // Calculate proportions and collision probability
    std::vector<double> p(k, 0.0);
    calc_proportions(data.symbols.data(), p, L);
    
    long double p_col = 0.0;
    calc_collision_proportion(p, p_col);
    
    assert(p_col >= 1.0L / ((long double)k));
    assert(p_col <= 1.0L);
    
    result.p_col = static_cast<double>(p_col);
    
    // Special case: if p_col == 1.0 (all same symbol)
    if (p_col > 1.0L - LDBL_EPSILON) {
        result.prob_x_ge_1 = 1.0;
        result.passed = true;
        return result;
    }
    
    assert(p_col < 1.0L);
    
    // Get the length of the longest repeated substring
    long int W;
    if (L < SAINDEX_MAX) {
        W = len_LRS32(data.symbols.data(), L);
    } else {
        W = len_LRS64(data.symbols.data(), L);
    }
    
    result.W = static_cast<int>(W);
    
    // Check that p_col^W is representable
    long double p_colPower = powl(p_col, (long double)W);
    assert(p_colPower >= LDBL_MIN);
    assert(p_colPower <= 1.0L - LDBL_EPSILON);
    
    // Calculate log probability of no collision per pair
    long double logProbNoColsPerPair = log1pl(-p_colPower);
    assert(logProbNoColsPerPair < 0.0L);
    
    // Number of pairs of W-length substrings
    uint64_t N = n_choose_2(static_cast<uint64_t>(L - W + 1));
    
    // Calculate probability of no collisions
    long double probNoCols = expl(((long double)N) * logProbNoColsPerPair);
    long double prob_x_ge_1 = 1.0L - probNoCols;
    
    result.prob_x_ge_1 = static_cast<double>(prob_x_ge_1);
    
    // Test passes if Pr(X >= 1) >= 1/1000
    // Equivalently: log(0.999) >= N * log1p(-p_col^W)
    result.passed = (logl(0.999L) >= ((long double)N) * logProbNoColsPerPair);
    
    return result;
}

lenLrsResult lenLrsTest(const std::filesystem::path& filepath) {
    std::ifstream file(filepath, std::ios::binary);
    if (!file) {
        throw std::runtime_error("lenLrsTest: failed to open file " + filepath.string());
    }

    std::vector<uint8_t> buffer(
        (std::istreambuf_iterator<char>(file)),
        std::istreambuf_iterator<char>());

    EntropyInputData data{buffer, 8};
    return lenLrsTest(data);
}

} // namespace lib90b