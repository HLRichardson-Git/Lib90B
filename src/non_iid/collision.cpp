
#include <fstream>
#include <filesystem>

#include <lib90B/collision.h>
#include "utils.h"

namespace lib90b {

// Section 6.3.2 - Collision Estimate
// Collision test is defined only for bitstring input
CollisionResult collisionTest(const EntropyInputData& data) {
    // Always expand to bits, ignore literal mode
    const std::vector<uint8_t>& symbols_to_use = data.getBits();
    const long len = static_cast<long>(symbols_to_use.size());
    
    if (len < 2) {
        throw std::runtime_error("collisionTest: insufficient data length");
    }

    CollisionResult result;

    long i = 0;
    long v = 0;
    double s = 0.0;
    int t_v = 0;

    // Compute wait times until collisions
    while (i < len - 1) {
        if (symbols_to_use[i] == symbols_to_use[i + 1]) t_v = 2;        // 00 or 11
        else if (i < len - 2) t_v = 3;                                  // patterns like 101, 011, etc.
        else break;

        v++;
        s += t_v * t_v;
        i += t_v;
    }

    result.v = static_cast<uint64_t>(v);
    result.sum_ti = static_cast<uint64_t>(i);

    // X is mean of t_v's, s is sample std deviation
    double X = i / static_cast<double>(v);
    result.x_bar = X;
    double sigma_hat = std::sqrt((s - i * X) / (v - 1));
    result.sigma_hat = sigma_hat;

    // Adjust X by confidence interval
    double X_prime = X - ZALPHA * sigma_hat / std::sqrt(v);
    if (X_prime < 2.0) X_prime = 2.0;
    result.x_bar_prime = X_prime;

    double p = 0.5;
    double entEst = 1.0;
    bool used_lower_bound = false;

    if (X_prime < 2.5) {
        p = 0.5 + std::sqrt(1.25 - 0.5 * X_prime);
        entEst = -std::log2(p);
        used_lower_bound = true;
    }

    result.p = p;
    result.h_bitstring = entEst;  // only bitstring is meaningful
    result.used_lower_bound = used_lower_bound;

    return result;
}

CollisionResult collisionTest(const std::filesystem::path& filepath) {
    // Read file as raw bytes
    std::ifstream file(filepath, std::ios::binary);
    if (!file) throw std::runtime_error("collisionTest: failed to open file " + filepath.string());

    std::vector<uint8_t> buffer((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
    EntropyInputData data{ buffer, 8 }; // default word size = 8
    return collisionTest(data);
}

} // namespace lib90b
