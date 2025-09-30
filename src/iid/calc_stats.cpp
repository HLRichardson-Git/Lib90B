
#include <vector>
#include <algorithm>
#include <cstdint>
#include <stdexcept>
#include <fstream>
#include <iterator>
#include <unordered_set>
#include <unordered_map>

#include <lib90b/calc_stats.h>
#include <lib90b/entropy_tests.h>

namespace lib90b {

// Core function: works on prepared EntropyInputData
CalcStatsResult calcStats(const EntropyInputData& data) {
    if (data.symbols.empty()) {
        throw std::runtime_error("calcStats: no symbols provided");
    }

    // Mean: use raw symbols
    double sum = 0.0;
    for (uint8_t b : data.symbols) {
        sum += static_cast<double>(b);
    }
    double mean = sum / data.symbols.size();

    // Build translated symbols for median
    std::unordered_map<uint8_t, uint8_t> symbol_to_index;
    std::vector<uint8_t> sorted_unique(data.symbols.begin(), data.symbols.end());
    std::sort(sorted_unique.begin(), sorted_unique.end());
    sorted_unique.erase(std::unique(sorted_unique.begin(), sorted_unique.end()), sorted_unique.end());

    for (size_t i = 0; i < sorted_unique.size(); ++i) {
        symbol_to_index[sorted_unique[i]] = static_cast<uint8_t>(i);
    }

    std::vector<uint8_t> translated_symbols;
    translated_symbols.reserve(data.symbols.size());
    for (uint8_t b : data.symbols) {
        translated_symbols.push_back(symbol_to_index[b]);
    }

    // Median: translated symbols
    std::vector<uint8_t> sorted(translated_symbols.begin(), translated_symbols.end());
    std::sort(sorted.begin(), sorted.end());
    size_t n = sorted.size();
    double median;
    if (sorted_unique.size() == 2) {
        median = 0.5;  // special IID convention
    } else {
        if (n % 2 == 1) {
            median = static_cast<double>(sorted[n / 2]);
        } else {
            median = (static_cast<double>(sorted[n / 2]) +
                      static_cast<double>(sorted[n / 2 - 1])) / 2.0;
        }
    }

    return CalcStatsResult{mean, median};
}

// Convenience function: takes a file path
CalcStatsResult calcStats(const std::filesystem::path& filepath) {
    std::ifstream file(filepath, std::ios::binary);
    if (!file) {
        throw std::runtime_error("calcStats: failed to open file " + filepath.string());
    }

    std::vector<uint8_t> raw_symbols(
        (std::istreambuf_iterator<char>(file)),
        std::istreambuf_iterator<char>()
    );

    if (raw_symbols.empty()) {
        throw std::runtime_error("calcStats: file is empty " + filepath.string());
    }

    // Fill EntropyInputData: raw bytes in symbols
    EntropyInputData data;
    data.symbols = std::move(raw_symbols);
    data.alph_size = 0; // will be computed in the core function if needed

    return calcStats(data);
}

} // namespace lib90b
