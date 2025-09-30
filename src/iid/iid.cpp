
#include <fstream>
#include <filesystem>
#include <future>
#include <optional>
#include <semaphore>
#include <thread>
#include <unordered_set>
#include <vector>

#include <lib90b/entropy_tests.h>
#include <lib90b/most_common.h>
#include <lib90b/chi_square.h>
#include <lib90b/lrs.h>
#include <lib90b/permutation.h>
#include "utils.h"

namespace lib90b {

IidResult iidTestSuite(EntropyInputData& data) {
    IidResult result;

    // Track alphabet size if not set
    if (data.alph_size == 0) {
        std::unordered_set<uint8_t> unique_symbols(data.symbols.begin(), data.symbols.end());
        data.alph_size = unique_symbols.size();
    }

    // Limit concurrency
    const size_t max_concurrent = std::min(4u, std::thread::hardware_concurrency());
    auto sem = std::make_shared<std::counting_semaphore<10>>(max_concurrent);

    auto safe_async = [sem](auto&& func, const EntropyInputData& data_copy) {
        return std::async(std::launch::async, [sem, func, data_copy]() {
            sem->acquire();
            try {
                auto res = func(data_copy);
                sem->release();
                return res;
            } catch (...) {
                sem->release();
                throw;
            }
        });
    };

    // Launch tests
    auto fut_most_common = safe_async([](const auto& d) { return mostCommonTest(d); }, data);
    auto fut_chi_square  = safe_async([](const auto& d) { return chiSquareTest(d); }, data);
    auto fut_lrs         = safe_async([](const auto& d) { return lenLrsTest(d); }, data);
    auto fut_perm        = safe_async([](const auto& d) { return permutationTest(d); }, data);

    // Collect results
    result.most_common = fut_most_common.get();
    result.chi_square  = fut_chi_square.get();
    result.lrs         = fut_lrs.get();
    result.permutation = fut_perm.get();

    // Bitstring versions if necessary
    if (data.alph_size > 2) {
        auto fut_most_common_bit = safe_async([](const auto& d) { return mostCommonTest(d, SymbolMode::Bitstring); }, data);

        auto patch_bit = [](auto &orig, auto &fut) {
            auto bit_res = fut.get();
            if (bit_res.h_bitstring.has_value()) orig.h_bitstring = bit_res.h_bitstring;
        };

        patch_bit(result.most_common, fut_most_common_bit);
    }

    // Compute H_bitstring
    result.H_bitstring = std::min({
        result.most_common.h_bitstring.value_or(1.0)
    });

    // Compute H_original
    result.H_original = std::min({
        result.most_common.h_original.value_or((double)data.word_size)
    });

    // Overall min entropy
    double h_assessed = data.word_size;
    if (result.H_bitstring > 0.0) h_assessed = std::min(h_assessed, result.H_bitstring * data.word_size);
    if (result.H_original > 0.0) h_assessed = std::min(h_assessed, result.H_original);
    result.min_entropy = h_assessed;

    result.passed = (
        result.chi_square.passed &&
        result.lrs.passed &&
        result.permutation.passed
    );

    return result;
}

IidResult iidTestSuite(const std::filesystem::path& filepath) {
    std::ifstream file(filepath, std::ios::binary);
    if (!file) throw std::runtime_error("Failed to open file: " + filepath.string());

    std::vector<uint8_t> buffer((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
    EntropyInputData data{ buffer, 8 };

    IidResult result = iidTestSuite(data);
    result.filename = filepath;

    return result;
}

} // namespace lib90b