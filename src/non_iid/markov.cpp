// markov.cpp
#include <fstream>
#include <stdexcept>
#include <cmath>
#include <algorithm>

#include <lib90b/markov.h>
#include "utils.h"

namespace lib90b {

MarkovResult markovTest(const EntropyInputData& data) {
    const std::vector<uint8_t>& bits = data.getBits();
    const long len = static_cast<long>(bits.size());

    if (len <= 1) {
        throw std::runtime_error("markovTest: insufficient data length (need > 1 bit)");
    }

    long C_0 = 0;
    long C_00 = 0;
    long C_10 = 0;

    // Count unconditional and transitions
    for (long i = 0; i < len - 1; i++) {
        if (bits[i] == 0) {
            C_0++;
            if (bits[i + 1] == 0) C_00++;
        } else if (bits[i + 1] == 0) {
            C_10++;
        }
    }

    long C_1 = (len - 1) - C_0;

    double P_00 = (C_0 > 0) ? static_cast<double>(C_00) / C_0 : 0.0;
    double P_01 = (C_0 > 0) ? 1.0 - P_00 : 0.0;

    double P_10 = (C_1 > 0) ? static_cast<double>(C_10) / C_1 : 0.0;
    double P_11 = (C_1 > 0) ? 1.0 - P_10 : 0.0;

    // account for last symbol
    if (bits[len - 1] == 0) C_0++;

    double P_0 = static_cast<double>(C_0) / len;
    double P_1 = 1.0 - P_0;

    // Calculate min-entropy
    double H_min = 128.0;

    if (P_00 > 0.0) {
        H_min = std::min(H_min, -std::log2(P_0) - 127.0 * std::log2(P_00));
    }
    if (P_01 > 0.0 && P_10 > 0.0) {
        H_min = std::min(H_min, -std::log2(P_0) - 64.0 * std::log2(P_01) - 63.0 * std::log2(P_10));
    }
    if (P_01 > 0.0 && P_11 > 0.0) {
        H_min = std::min(H_min, -std::log2(P_0) - std::log2(P_01) - 126.0 * std::log2(P_11));
    }
    if (P_10 > 0.0 && P_00 > 0.0) {
        H_min = std::min(H_min, -std::log2(P_1) - std::log2(P_10) - 126.0 * std::log2(P_00));
    }
    if (P_10 > 0.0 && P_01 > 0.0) {
        H_min = std::min(H_min, -std::log2(P_1) - 64.0 * std::log2(P_10) - 63.0 * std::log2(P_01));
    }
    if (P_11 > 0.0) {
        H_min = std::min(H_min, -std::log2(P_1) - 127.0 * std::log2(P_11));
    }

    double p_hat_max = std::pow(2.0, -H_min);
    double entEst = std::min(H_min / 128.0, 1.0);

    // Fill result
    MarkovResult result;
    result.p_0 = P_0;
    result.p_1 = P_1;
    result.p_00 = P_00;
    result.p_01 = P_01;
    result.p_10 = P_10;
    result.p_11 = P_11;
    result.p_hat_max = p_hat_max;

    result.h_bitstring = entEst; // always bitstring
    return result;
}

MarkovResult markovTest(const std::filesystem::path& filepath) {
    std::ifstream file(filepath, std::ios::binary);
    if (!file) throw std::runtime_error("markovTest: failed to open file " + filepath.string());

    std::vector<uint8_t> buffer(
        (std::istreambuf_iterator<char>(file)),
        std::istreambuf_iterator<char>());

    EntropyInputData data{buffer, 8};
    return markovTest(data);
}

} // namespace lib90b
