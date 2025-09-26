
#include <vector>
#include <map>
#include <array>
#include <cstring>
#include <fstream>
#include <stdexcept>
#include <algorithm>

#include <lib90b/lz78y.h>
#include "utils.h"

#define B_len 16
#define MAX_DICTIONARY_SIZE 65536

PredictionEstimateResult binaryLZ78YPredictionEstimate(const std::vector<uint8_t>& bits) {
    const long L = static_cast<long>(bits.size());

    if (L <= B_len + 2) {
        throw std::runtime_error("LZ78Y: not enough samples (need > B_len+2)");
    }

    long* binaryDict[B_len];
    long curRunOfCorrects = 0;
    long maxRunOfCorrects = 0;
    long correctCount = 0;
    uint32_t curPattern = 0;
    long dictElems = 0;

    // Initialize
    for (int j = 0; j < B_len; j++) {
        binaryDict[j] = new long[1U << (j + 2)];
        std::memset(binaryDict[j], 0, sizeof(long) * (1U << (j + 2)));
    }

    // Initialize B counts
    for (int j = 0; j < B_len; j++) {
        curPattern |= ((bits[B_len - j - 1] & 1U) << j);
        auto entry = &binaryDict[j][(curPattern << 1)];
        entry[bits[B_len] & 1U] = 1;
        dictElems++;
    }

    // Process stream
    for (long i = B_len + 1; i < L; i++) {
        curPattern = compressedBitSymbols(bits.data() + i - B_len, B_len);

        bool havePrediction = false;
        uint8_t curPrediction = 2;
        long maxCount = 0;

        for (int j = B_len; j > 0; j--) {
            curPattern &= ((1U << j) - 1);

            long* binaryDictEntry = &binaryDict[j - 1][(curPattern << 1)];

            int roundPrediction;
            long curCount;
            if (binaryDictEntry[0] > binaryDictEntry[1]) {
                roundPrediction = 0;
                curCount = binaryDictEntry[0];
            } else {
                roundPrediction = 1;
                curCount = binaryDictEntry[1];
            }

            if (curCount > 0) {
                if (curCount > maxCount) {
                    maxCount = curCount;
                    curPrediction = roundPrediction;
                    havePrediction = true;
                }
                binaryDictEntry[bits[i] & 1U]++;
            } else if (dictElems < MAX_DICTIONARY_SIZE) {
                binaryDictEntry[bits[i] & 1U] = 1;
                dictElems++;
            }
        }

        if (havePrediction && (curPrediction == bits[i])) {
            correctCount++;
            if (++curRunOfCorrects > maxRunOfCorrects) {
                maxRunOfCorrects = curRunOfCorrects;
            }
        } else {
            curRunOfCorrects = 0;
        }
    }

    for (int j = 0; j < B_len; j++) {
        delete[] binaryDict[j];
    }

    return predictionEstimate(correctCount, L - B_len - 1, maxRunOfCorrects, 2);
}

namespace lib90b {

// Section 6.3.10 - LZ78Y Prediction Estimate
Lz78yResult lz78yTest(const EntropyInputData& data, SymbolMode mode) {
    const std::vector<uint8_t>& symbols_to_use = (mode == SymbolMode::Bitstring) ? data.getBits() : data.symbols;
    const long len = static_cast<long>(symbols_to_use.size());
    const int alph_size = (mode == SymbolMode::Bitstring) ? 2 : (1 << data.word_size);

    if (alph_size == 2) {
        PredictionEstimateResult pred = binaryLZ78YPredictionEstimate(symbols_to_use);

        Lz78yResult result;
        static_cast<PredictionEstimateResult&>(result) = pred;
        result.h_original = -std::log2(pred.curMax);
        if (mode == SymbolMode::Bitstring) result.h_bitstring = result.h_original;
        return result;
    }

    // General dictionary case
    if (len < B_len + 2) {
        throw std::runtime_error("LZ78Y: not enough samples to run test");
    }

    using Key = std::array<uint8_t, B_len>;
    std::array<std::map<Key, PostfixDictionary>, B_len> D;
    int dict_size = 0;
    long C = 0;
    long run_len = 0;
    long max_run_len = 0;
    Key x{};

    // Initialize
    for (int j = 1; j <= B_len; j++) {
        std::memcpy(x.data(), symbols_to_use.data() + B_len - j, j);
        D[j - 1][x].incrementPostfix(symbols_to_use[B_len], true);
        dict_size++;
    }

    // Process stream
    for (long i = B_len + 1; i < len; i++) {
        bool havePrediction = false;
        uint8_t prediction = 0;
        long max_count = 0;

        for (int j = B_len; j > 0; j--) {
            std::memset(x.data(), 0, B_len);
            std::memcpy(x.data(), symbols_to_use.data() + i - j, j);

            auto curp = D[j - 1].find(x);
            if (curp != D[j - 1].end()) {
                long count;
                uint8_t y = curp->second.predict(count);

                if (count > max_count) {
                    max_count = count;
                    prediction = y;
                    havePrediction = true;
                }

                curp->second.incrementPostfix(symbols_to_use[i], true);
            } else if (dict_size < MAX_DICTIONARY_SIZE) {
                D[j - 1][x].incrementPostfix(symbols_to_use[i], true);
                dict_size++;
            }
        }

        if (havePrediction && prediction == symbols_to_use[i]) {
            C++;
            if (++run_len > max_run_len) max_run_len = run_len;
        } else {
            run_len = 0;
        }
    }

    PredictionEstimateResult pred = predictionEstimate(C, len - B_len - 1, max_run_len, alph_size);

    Lz78yResult result;
    static_cast<PredictionEstimateResult&>(result) = pred;
    result.h_original = -std::log2(pred.curMax);
    if (mode == SymbolMode::Bitstring) result.h_bitstring = result.h_original;
    return result;
}


// Public API: file-based
Lz78yResult lz78yTest(const std::filesystem::path& filepath, SymbolMode mode) {
    std::ifstream file(filepath, std::ios::binary);
    if (!file) throw std::runtime_error("lz78yTest: failed to open file " + filepath.string());

    std::vector<uint8_t> buffer(
        (std::istreambuf_iterator<char>(file)),
        std::istreambuf_iterator<char>());

    EntropyInputData data{buffer, 8};
    return lz78yTest(data, mode);
}

} // namespace lib90b
