
#include <vector>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <algorithm>
#include <stdexcept>

#include <lib90b/entropy_tests.h>
#include "utils.h"

#define D_LAG 128U
#define LAGMASK (D_LAG-1)

/*Convention:
 * if start==end, the buffer is empty. Start and end values are not ANDed, so range from 0...255
 * This extended range allows for use of all the entries of the buffer
 * See https://www.snellman.net/blog/archive/2016-12-13-ring-buffers/
 * The actual data locations range 0...127
 * start&LAGMASK is the index of the oldest data
 * end&LAGMASK is the index where the *next* data goes.
 */

struct lagBuf {
	uint8_t start;
	uint8_t end;
	long buf[D_LAG];
};

namespace lib90b {

LagResult lagTest(const EntropyInputData& data, SymbolMode mode) {
    // Choose symbol representation
    const std::vector<uint8_t>& symbols_to_use = (mode == SymbolMode::Bitstring) ? data.getBits() : data.symbols;
    const long L = static_cast<long>(symbols_to_use.size());
    const int k = (mode == SymbolMode::Bitstring) ? 2 : (1 << data.word_size);

    if (L <= 2) {
        throw std::runtime_error("lagTest: insufficient data length");
    }
    if (k < 2) {
        throw std::runtime_error("lagTest: invalid alphabet size");
    }

    // --- original lag_test implementation ---
    long scoreboard[D_LAG] = {0};
    int winner = 0;
    long curRunOfCorrects = 0;
    long maxRunOfCorrects = 0;
    long correctCount = 0;

    std::vector<lagBuf> ringBuffers(k);

    // Init first symbol
    ringBuffers[symbols_to_use[0]].buf[0] = 0;
    ringBuffers[symbols_to_use[0]].start = 0;
    ringBuffers[symbols_to_use[0]].end = 1;

    for (long i = 1; i < L; i++) {
        const uint8_t curSymbol = symbols_to_use[i];
        lagBuf* const curRingBuffer = &ringBuffers[curSymbol];

        // Check prediction
        if (curSymbol == symbols_to_use[i - winner - 1]) {
            correctCount++;
            curRunOfCorrects++;
            if (curRunOfCorrects > maxRunOfCorrects) {
                maxRunOfCorrects = curRunOfCorrects;
            }
        } else {
            curRunOfCorrects = 0;
        }

        // Update counters
        if (curRingBuffer->start != curRingBuffer->end) {
            uint8_t counterIndex = curRingBuffer->end;
            const long cutoff = (i >= D_LAG) ? (i - D_LAG) : 0;

            do {
                counterIndex--;
                if (curRingBuffer->buf[counterIndex & LAGMASK] >= cutoff) {
                    long curOffset = i - curRingBuffer->buf[counterIndex & LAGMASK] - 1;
                    assert(curOffset < D_LAG);
                    long curScore = ++scoreboard[curOffset];

                    if (curScore >= scoreboard[winner]) {
                        winner = curOffset;
                    }
                } else {
                    curRingBuffer->start = static_cast<uint8_t>(counterIndex + 1U);
                    break;
                }
            } while (counterIndex != curRingBuffer->start);
        }

        // Add symbol
        if ((uint8_t)(curRingBuffer->end - curRingBuffer->start) == D_LAG) {
            curRingBuffer->start++;
        }
        curRingBuffer->buf[(curRingBuffer->end) & LAGMASK] = i;
        curRingBuffer->end++;
        assert((uint8_t)(curRingBuffer->end - curRingBuffer->start) <= D_LAG);
    }

    // Use shared prediction estimate
    PredictionEstimateResult pred = predictionEstimate(correctCount, L - 1, maxRunOfCorrects, k);

    LagResult result;
    static_cast<PredictionEstimateResult&>(result) = pred;

    double entEst = -std::log2(pred.curMax);
    if (mode == SymbolMode::Original) result.h_original = entEst;
    else result.h_bitstring = entEst;

    return result;
}

LagResult lagTest(const std::filesystem::path& filepath, SymbolMode mode) {
    std::ifstream file(filepath, std::ios::binary);
    if (!file) throw std::runtime_error("lagTest: failed to open file " + filepath.string());

    std::vector<uint8_t> buffer(
        (std::istreambuf_iterator<char>(file)),
        std::istreambuf_iterator<char>()
    );

    EntropyInputData data{ buffer, 8 };
    return lagTest(data, mode);
}

} // namespace lib90b

