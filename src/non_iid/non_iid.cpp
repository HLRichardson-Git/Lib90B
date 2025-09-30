
#pragma once

#include <future>
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <string>
#include <optional>
#include <thread>
#include <semaphore>
#include <unordered_set>

#include <lib90b/non_iid.h>
#include <lib90b/entropy_tests.h>
#include <lib90b/most_common.h>
#include <lib90b/collision.h>
#include <lib90b/markov.h>
#include <lib90b/compression.h>
#include <lib90b/t_tuple.h>
#include <lib90b/lrs.h>
#include <lib90b/multi_mcw.h>
#include <lib90b/lag.h>
#include <lib90b/multi_mmc.h>
#include <lib90b/lz78y.h>
#include "utils.h"
#include "TestRunUtils.h"

namespace lib90b {

NonIidResult nonIidTestSuite(EntropyInputData& data) {
    NonIidResult result;
    result.word_size = data.word_size;

    // Track alphabet size if not set
    if (data.alph_size == 0) {
        std::unordered_set<uint8_t> unique_symbols(data.symbols.begin(), data.symbols.end());
        data.alph_size = unique_symbols.size();
    } else {
        data.alph_size = data.alph_size;
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
    auto fut_compression = safe_async([](const auto& d) { return compressionTest(d); }, data);
    auto fut_lz78y       = safe_async([](const auto& d) { return lz78yTest(d); }, data);
    auto fut_multi_mmc   = safe_async([](const auto& d) { return multiMmcTest(d); }, data);
    auto fut_lag         = safe_async([](const auto& d) { return lagTest(d); }, data);
    auto fut_lrs         = safe_async([](const auto& d) { return lrsTest(d); }, data);
    auto fut_multi_mcw   = safe_async([](const auto& d) { return multiMcwTest(d); }, data);
    auto fut_t_tuple     = safe_async([](const auto& d) { return tTupleTest(d); }, data);
    auto fut_markov      = safe_async([](const auto& d) { return markovTest(d); }, data);
    auto fut_most_common = safe_async([](const auto& d) { return mostCommonTest(d); }, data);
    auto fut_collision   = safe_async([](const auto& d) { return collisionTest(d); }, data);

    // Collect results
    result.compression = fut_compression.get();
    result.lz78y       = fut_lz78y.get();
    result.multi_mmc   = fut_multi_mmc.get();
    result.lag         = fut_lag.get();
    result.lrs         = fut_lrs.get();
    result.multi_mcw   = fut_multi_mcw.get();
    result.t_tuple     = fut_t_tuple.get();
    result.markov      = fut_markov.get();
    result.most_common = fut_most_common.get();
    result.collision   = fut_collision.get();

    // Run bitstring versions if necessary
    if (data.alph_size > 2) {
        auto fut_lz78y_bit       = safe_async([](const auto& d) { return lz78yTest(d, SymbolMode::Bitstring); }, data);
        auto fut_multi_mmc_bit   = safe_async([](const auto& d) { return multiMmcTest(d, SymbolMode::Bitstring); }, data);
        auto fut_lag_bit         = safe_async([](const auto& d) { return lagTest(d, SymbolMode::Bitstring); }, data);
        auto fut_lrs_bit         = safe_async([](const auto& d) { return lrsTest(d, SymbolMode::Bitstring); }, data);
        auto fut_multi_mcw_bit   = safe_async([](const auto& d) { return multiMcwTest(d, SymbolMode::Bitstring); }, data);
        auto fut_t_tuple_bit     = safe_async([](const auto& d) { return tTupleTest(d, SymbolMode::Bitstring); }, data);
        auto fut_most_common_bit = safe_async([](const auto& d) { return mostCommonTest(d, SymbolMode::Bitstring); }, data);

        auto patch = [](auto &orig, auto &fut) {
            auto bit_res = fut.get();
            if (bit_res.h_bitstring.has_value())
                orig.h_bitstring = bit_res.h_bitstring;
        };

        patch(result.lz78y,       fut_lz78y_bit);
        patch(result.multi_mmc,   fut_multi_mmc_bit);
        patch(result.lag,         fut_lag_bit);
        patch(result.lrs,         fut_lrs_bit);
        patch(result.multi_mcw,   fut_multi_mcw_bit);
        patch(result.t_tuple,     fut_t_tuple_bit);
        patch(result.most_common, fut_most_common_bit);
    }

    // Compute H_bitstring
    result.H_bitstring = std::min({
        result.most_common.h_bitstring.value_or(1.0),
        result.collision.h_bitstring.value_or(1.0),
        result.markov.h_bitstring.value_or(1.0),
        result.compression.h_bitstring.value_or(1.0),
        result.t_tuple.h_bitstring.value_or(1.0),
        result.lrs.h_bitstring.value_or(1.0),
        result.multi_mcw.h_bitstring.value_or(1.0),
        result.lag.h_bitstring.value_or(1.0),
        result.multi_mmc.h_bitstring.value_or(1.0),
        result.lz78y.h_bitstring.value_or(1.0)
    });

    // Compute H_original
    result.H_original = std::min({
        result.most_common.h_original.value_or((double)data.word_size),
        result.collision.h_original.value_or((double)data.word_size),
        result.markov.h_original.value_or((double)data.word_size),
        result.compression.h_original.value_or((double)data.word_size),
        result.t_tuple.h_original.value_or((double)data.word_size),
        result.lrs.h_original.value_or((double)data.word_size),
        result.multi_mcw.h_original.value_or((double)data.word_size),
        result.lag.h_original.value_or((double)data.word_size),
        result.multi_mmc.h_original.value_or((double)data.word_size),
        result.lz78y.h_original.value_or((double)data.word_size)
    });

    // Compute overall min_entropy
    double h_assessed = data.word_size;
    if (result.H_bitstring.has_value()) h_assessed = std::min(h_assessed, result.H_bitstring.value() * data.word_size);
    if (result.H_original.has_value()) h_assessed = std::min(h_assessed, result.H_original.value());
    result.min_entropy = h_assessed;

    return result;
}

NonIidResult nonIidTestSuite(const std::filesystem::path& filepath) {
    std::ifstream file(filepath, std::ios::binary);
    if (!file) {
        throw std::runtime_error("Failed to open file: " + filepath.string());
    }

    std::vector<uint8_t> buffer(
        (std::istreambuf_iterator<char>(file)),
        std::istreambuf_iterator<char>());

    EntropyInputData data{ buffer, 8 };
    NonIidResult result = nonIidTestSuite(data);
    result.filename = filepath;
    return result;
}

} // namespace lib90b
