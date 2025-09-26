// multi_mcw.cpp
#include <fstream>
#include <stdexcept>
#include <vector>
#include <algorithm>
#include <array>
#include <limits>

#include <lib90b/multi_mcw.h>
#include "utils.h"

static constexpr int NUM_WINS = 4;
static constexpr std::array<int, NUM_WINS> W = {63, 255, 1023, 4095};

namespace lib90b {

MultiMcwResult multiMcwTest(const EntropyInputData& data, SymbolMode mode) {
    // Select symbols (bitstring or original symbols)
    const std::vector<uint8_t>& symbols_to_use = (mode == SymbolMode::Bitstring) ? data.getBits() : data.symbols;
    const long len = static_cast<long>(symbols_to_use.size());
    const int alph_size = (mode == SymbolMode::Bitstring) ? 2 : (1 << data.word_size);

    if (len < W.back() + 1) {
        throw std::runtime_error("multiMcwTest: insufficient data length");
    }

    long N = len - W[0];
    long C = 0;
    long run_len = 0;
    long max_run_len = 0;

    int winner = 0;
    std::array<long, NUM_WINS> scoreboard = {0, 0, 0, 0};
    std::array<long, NUM_WINS> max_cnts = {0, 0, 0, 0};
    std::array<uint8_t, NUM_WINS> frequent = {0, 0, 0, 0};

    // Window counts and positions
    std::vector<std::vector<long>> win_cnts(NUM_WINS, std::vector<long>(alph_size, 0));
    std::vector<std::vector<long>> win_poses(NUM_WINS, std::vector<long>(alph_size, 0));

    // Initialize window counts
    for (long i = 0; i < W.back(); i++) {
        for (int j = 0; j < NUM_WINS; j++) {
            if (i < W[j]) {
                if (++win_cnts[j][symbols_to_use[i]] >= max_cnts[j]) {
                    max_cnts[j] = win_cnts[j][symbols_to_use[i]];
                    frequent[j] = symbols_to_use[i];
                }
                win_poses[j][symbols_to_use[i]] = i;
            }
        }
    }

    // Perform predictions
    for (long i = W[0]; i < len; i++) {
        // Test prediction of winner
        if (frequent[winner] == symbols_to_use[i]) {
            C++;
            if (++run_len > max_run_len) max_run_len = run_len;
        } else {
            run_len = 0;
        }

        // Update scoreboard and select new winner
        for (int j = 0; j < NUM_WINS; j++) {
            if (i >= W[j] && frequent[j] == symbols_to_use[i]) {
                if (++scoreboard[j] >= scoreboard[winner]) {
                    winner = j;
                }
            }
        }

        // Update window counts and select new frequent values
        for (int j = 0; j < NUM_WINS; j++) {
            if (i >= W[j]) {
                win_cnts[j][symbols_to_use[i - W[j]]]--;
                win_cnts[j][symbols_to_use[i]]++;
                win_poses[j][symbols_to_use[i]] = i;

                if ((symbols_to_use[i - W[j]] != frequent[j]) &&
                    (win_cnts[j][symbols_to_use[i]] >= max_cnts[j])) {
                    max_cnts[j] = win_cnts[j][symbols_to_use[i]];
                    frequent[j] = symbols_to_use[i];
                } else if (symbols_to_use[i - W[j]] == frequent[j]) {
                    max_cnts[j]--;
                    long max_pos = i - W[j];
                    for (int k = 0; k < alph_size; k++) {
                        if ((win_cnts[j][k] > max_cnts[j]) ||
                            ((win_cnts[j][k] == max_cnts[j]) && (win_poses[j][k] >= max_pos))) {
                            max_cnts[j] = win_cnts[j][k];
                            frequent[j] = static_cast<uint8_t>(k);
                            max_pos = win_poses[j][k];
                        }
                    }
                }
            }
        }
    }

    // Run prediction estimate
    PredictionEstimateResult predRes = predictionEstimate(C, N, max_run_len, alph_size);

    MultiMcwResult result;
    static_cast<PredictionEstimateResult&>(result) = predRes;

    double entEst = -std::log2(predRes.curMax);
    if (mode == SymbolMode::Original) result.h_original = entEst;
    else result.h_bitstring = entEst;

    return result;
}

MultiMcwResult multiMcwTest(const std::filesystem::path& filepath, SymbolMode mode) {
    std::ifstream file(filepath, std::ios::binary);
    if (!file) throw std::runtime_error("multiMcwTest: failed to open file " + filepath.string());

    std::vector<uint8_t> buffer(
        (std::istreambuf_iterator<char>(file)),
        std::istreambuf_iterator<char>());

    EntropyInputData data{buffer, 8};
    return multiMcwTest(data, mode);
}

} // namespace lib90b
