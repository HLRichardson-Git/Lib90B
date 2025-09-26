#pragma once

#include <array>
#include <map>
#include <vector>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <cstring>
#include <stdexcept>
#include <cmath>
#include <algorithm>

#include <lib90b/entropy_tests.h>
#include "utils.h"

constexpr int D_MMC = 16;
constexpr long MAX_ENTRIES = 100000;

namespace lib90b {

static MultiMmcResult binaryMultiMmc(const std::vector<uint8_t>& S) {
    const long L = static_cast<long>(S.size());
    if (L < 3) throw std::runtime_error("binaryMultiMmc: insufficient data");

    long scoreboard[D_MMC] = {0};
    long winner = 0;
    long curRun = 0, maxRun = 0, correctCount = 0;
    long dictElems[D_MMC] = {0};
    uint32_t curPattern = 0;

    // Allocate binary dictionaries
    std::array<std::vector<long>, D_MMC> binaryDict;
    for (int j = 0; j < D_MMC; ++j) {
        binaryDict[j].resize(1U << (j + 2), 0);
    }

    // Initialize dictionaries with first D_MMC symbols
    for (int d = 0; d < D_MMC && d < L - 1; ++d) {
        curPattern = (curPattern << 1) | (S[d] & 1);
        binaryDict[d][curPattern * 2 + (S[d + 1] & 1)] = 1;
        dictElems[d] = 1;
    }

    // Prediction loop
    for (long i = 2; i < L; ++i) {
        bool found_x = false;
        long curWinner = winner;
        curPattern = 0;

        for (int d = 0; d < D_MMC && d <= i - 2; ++d) {
            // Build pattern for this d value
            curPattern |= static_cast<uint32_t>(S[i - d - 1] & 1) << d;
            
            auto& dictEntry = binaryDict[d];
            long* pairCounts = &dictEntry[curPattern * 2];
            
            uint8_t curPrediction = 2; // Invalid prediction initially
            long curCount = 0;

            // Check if we should make a prediction
            if (d == 0 || found_x) {
                // Make prediction based on counts
                if (pairCounts[0] > pairCounts[1]) {
                    curPrediction = 0;
                    curCount = pairCounts[0];
                } else {
                    curPrediction = 1;
                    curCount = pairCounts[1];
                }

                found_x = (curCount > 0);
            }

            if (found_x) {
                // Check if prediction is correct
                if (curPrediction == S[i]) {
                    ++scoreboard[d];
                    if (scoreboard[d] >= scoreboard[winner]) winner = d;

                    if (d == curWinner) {
                        ++correctCount;
                        if (++curRun > maxRun) maxRun = curRun;
                    }
                } else if (d == curWinner) {
                    curRun = 0;
                }

                // Update dictionary
                if (pairCounts[S[i] & 1] != 0) {
                    ++pairCounts[S[i] & 1];
                } else if (dictElems[d] < MAX_ENTRIES) {
                    pairCounts[S[i] & 1] = 1;
                    ++dictElems[d];
                }
            } else if (dictElems[d] < MAX_ENTRIES) {
                // Pattern not found, so add it
                pairCounts[S[i] & 1] = 1;
                ++dictElems[d];
            }
        }
    }

    PredictionEstimateResult pred = predictionEstimate(correctCount, L - 2, maxRun, 2);
    MultiMmcResult result;
    static_cast<PredictionEstimateResult&>(result) = pred;
    result.h_bitstring = -std::log2(pred.curMax);
    return result;
}

// Section 6.3.9 - MultiMMC Prediction Estimate
/* This implementation of the MultiMMC test is a based on NIST's really cleaver implementation,
 * which interleaves the predictions and updates. This makes optimization much easier.
 * It's opaque why it works correctly (in particular, the first few symbols are added in a different
 * order than in the reference implementation), but once the initialization is performed, the rest of the
 * operations are done in the correct order.
 * The general observations that explains why this approach works are that
 * 1) each prediction that could succeed (i.e., ignoring some of the early predictions that must fail due
 *    to lack of strings of the queried length) must occur only after all the correct (x,y) tuples for that
 *    length have been processed. One is free to reorder otherwise.
 * 2) If there is a distinct string of length n, then this induces corresponding unique strings of all
 *    lengths greater than n. We track all string lengths independently (thus conceptually, we
 *    could run out of a short-string length prior to a long string length, thus erroneously not add
 *    some long string to the dictionary after no longer looking for a string to the dictionary when
 *    we should have), this can't happen in practice because we add strings from shortest to longest.
 */
MultiMmcResult multiMmcTest(const EntropyInputData& data, SymbolMode mode) {
    const std::vector<uint8_t>& symbols = (mode == SymbolMode::Bitstring) ? data.getBits() : data.symbols;

    if (symbols.size() < 3) throw std::runtime_error("multiMmcTest: insufficient data");

    // Binary case
    if (data.word_size == 1 || mode == SymbolMode::Bitstring) {
        return binaryMultiMmc(symbols);
    }

    // Multi-symbol case
    const long L = static_cast<long>(symbols.size());
    long scoreboard[D_MMC] = {0};
    long winner = 0;
    long curRun = 0, maxRun = 0, correctCount = 0;
    long entries[D_MMC] = {0};
    std::array<uint8_t, D_MMC> x{};
    std::array<std::map<std::array<uint8_t, D_MMC>, PostfixDictionary>, D_MMC> M;

    // Initialize first D_MMC prefixes
    std::fill(x.begin(), x.end(), 0);
    for (int d = 0; d < D_MMC && d < L - 1; ++d) {
        std::copy_n(symbols.begin(), d + 1, x.begin());
        M[d][x].incrementPostfix(symbols[d + 1], true);
        entries[d] = 1;
    }

    // Prediction loop
    for (long i = 2; i < L; ++i) {
        bool found_x = false;
        long curWinner = winner;
        std::fill(x.begin(), x.end(), 0);

        for (int d = 0; d < D_MMC && d <= i - 2; ++d) {
            // Build the pattern for this d value
            std::copy_n(symbols.begin() + i - d - 1, d + 1, x.begin());

            // Check if we should make a prediction
            if (d == 0 || found_x) {
                auto curp = M[d].find(x);
                found_x = (curp != M[d].end());

                if (found_x) {
                    long predictCount;
                    uint8_t pred = curp->second.predict(predictCount);

                    // Check if prediction is correct
                    if (pred == symbols[i]) {
                        ++scoreboard[d];
                        if (scoreboard[d] >= scoreboard[winner]) winner = d;

                        if (d == curWinner) {
                            ++correctCount;
                            if (++curRun > maxRun) maxRun = curRun;
                        }
                    } else if (d == curWinner) {
                        curRun = 0;
                    }

                    // Update dictionary
                    if (curp->second.incrementPostfix(symbols[i], entries[d] < MAX_ENTRIES)) {
                        ++entries[d];
                    }
                }
            }

            if (!found_x && entries[d] < MAX_ENTRIES) {
                // Pattern not found, so add it
                M[d][x].incrementPostfix(symbols[i], true);
                ++entries[d];
            }
        }
    }

    PredictionEstimateResult pred = predictionEstimate(correctCount, L - 2, maxRun, 1U << data.word_size);
    MultiMmcResult result;
    static_cast<PredictionEstimateResult&>(result) = pred;
    result.h_original = -std::log2(pred.curMax);

    return result;
}

// File overload
MultiMmcResult multiMmcTest(const std::filesystem::path& filepath, SymbolMode mode) {
    std::ifstream file(filepath, std::ios::binary);
    if (!file) throw std::runtime_error("multiMmcTest: failed to open " + filepath.string());

    std::vector<uint8_t> buffer(
        (std::istreambuf_iterator<char>(file)),
        std::istreambuf_iterator<char>()
    );

    EntropyInputData data{buffer, 8};
    return multiMmcTest(data, mode);
}

} // namespace lib90b
