
#include <vector>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <algorithm> // for std::min
#include <stdexcept>

#include <lib90b/most_common.h>
#include <lib90b/entropy_tests.h>
#include "utils.h"

namespace lib90b {

// Section 6.3.1 - Most Common Value Estimate
MostCommonResult mostCommonTest(const EntropyInputData& data, SymbolMode mode) {

    const std::vector<uint8_t>& symbols_to_use = (mode == SymbolMode::Bitstring) ? data.getBits() : data.symbols;
    const long len = static_cast<long>(symbols_to_use.size());
    const int alph_size = (mode == SymbolMode::Bitstring) ? 2 : (1 << data.word_size);

    if (len <= 1) {
        throw std::runtime_error("runMostCommonTest: insufficient data length");
    }

    // Count frequencies
    std::vector<long> counts(alph_size, 0);
    for (long i = 0; i < len; i++) {
        counts[symbols_to_use[i]]++;
    }

    // Find mode
    long mode_count = 0;
    for (int i = 0; i < alph_size; i++) {
        if (counts[i] > mode_count) mode_count = counts[i];
    }

    double p_hat = mode_count / static_cast<double>(len);
    double p_u = std::min(
        1.0,
        p_hat + ZALPHA * std::sqrt(p_hat * (1.0 - p_hat) / (len - 1.0))
    );

    double entEst = -std::log2(p_u);

    MostCommonResult result;
    result.mode_count = static_cast<uint64_t>(mode_count);
    result.p_hat = p_hat;
    result.p_u = p_u;
    result.h_original = entEst;
    if (mode == SymbolMode::Bitstring) result.h_bitstring = entEst;

    return result;
}

MostCommonResult mostCommonTest(const std::filesystem::path& filepath, SymbolMode mode) {
    // Read file as raw bytes
    std::ifstream file(filepath, std::ios::binary);
    if (!file) throw std::runtime_error("mostCommonTest: failed to open file " + filepath.string());

    std::vector<uint8_t> buffer(
        (std::istreambuf_iterator<char>(file)),
        std::istreambuf_iterator<char>()
    );

    // Default word_size = 8 (byte-oriented)
    EntropyInputData data{ buffer, 8 };
    return mostCommonTest(data, mode);
}

} // namespace lib90b