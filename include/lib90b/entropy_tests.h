
#pragma once

#include <optional>
#include <cstdint>

#include "utils.h"

namespace lib90b {

enum class SymbolMode {
    Original,   // Use symbols as-is
    Bitstring   // Expand symbols to bits
};

struct EntropyInputData {
    std::vector<uint8_t> symbols;
    int word_size;
    size_t alph_size = 0;

    // bits generated on-demand
    mutable std::optional<std::vector<uint8_t>> bits;

    const std::vector<uint8_t>& getBits() const {
        if (!bits) {
            bits = expandSymbolsToBits(symbols, word_size);
        }
        return *bits;
    }
};

// Base struct
struct EntropyTestResult {
    std::optional<double> h_original;   // Entropy per symbol
    std::optional<double> h_bitstring;  // Entropy per bit
    virtual ~EntropyTestResult() = default;
};

// structs for each statistical test
struct MostCommonResult : public EntropyTestResult {
    uint64_t mode_count = 0;
    double p_hat = 0.0;
    double p_u = 0.0;
};

struct LrsResult : public EntropyTestResult {
    int u = -1;                 // smallest substring length considered
    int v = -1;                 // longest repeated substring length
    long double p_hat = 0.0;    // estimated probability
    long double p_u = 0.0;      // upper bound probability
    long double lrs_res = -1.0; // min-entropy estimate
};

struct CollisionResult : public EntropyTestResult {
    uint64_t v = 0;
    uint64_t sum_ti = 0;
    double x_bar = 0.0;
    double sigma_hat = 0.0;
    double x_bar_prime = 0.0;
    double p = 0.0;
    bool used_lower_bound = false;
};

struct CompressionResult : public EntropyTestResult {
    double x_bar = 0.0;
    double sigma_hat = 0.0;
    double x_bar_prime = 0.0;
    double p = 0.0;
    bool found_p = false;
    bool success = true;
};

struct LagResult : public EntropyTestResult, public PredictionEstimateResult {};

struct Lz78yResult : public EntropyTestResult, public PredictionEstimateResult {};

struct MarkovResult : public EntropyTestResult {
    double p_0 = 0.0;
    double p_1 = 0.0;
    double p_00 = 0.0;
    double p_01 = 0.0;
    double p_10 = 0.0;
    double p_11 = 0.0;
    double p_hat_max = 0.0;
};

struct MultiMcwResult : public EntropyTestResult, public PredictionEstimateResult {};

struct MultiMmcResult : public EntropyTestResult, public PredictionEstimateResult {};

struct TtupleResult : public EntropyTestResult {
    uint64_t t = 0;
    double p_hat_max = 0.0;
    double p_u = 0.0;
    double t_tuple_res = -1.0;
};

} // namespace lib90b
