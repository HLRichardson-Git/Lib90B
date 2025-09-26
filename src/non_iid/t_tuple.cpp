#pragma once

#include <vector>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <stdexcept>

#include <lib90b/entropy_tests.h>
#include "shared/lrs/lrs_algs.h"
#include "utils.h"

namespace lib90b {

// --- t-Tuple API using InternalLrsResult ---
TtupleResult tTupleTest(const EntropyInputData& data, SymbolMode mode) {
    const auto& symbols = (mode == SymbolMode::Bitstring) ? data.getBits() : data.symbols;

    InternalLrsResult internal = SAalgs(symbols.data(), symbols.size(),
                                        (mode == SymbolMode::Bitstring) ? 2 : (1 << data.word_size));

    TtupleResult result;
    result.t = internal.t;
    result.p_hat_max = internal.t_p_hat_max;
    result.p_u = internal.t_p_u;
    result.t_tuple_res = internal.t_tuple_res;

    if (mode == SymbolMode::Original) result.h_original = internal.t_tuple_res;
    else result.h_bitstring = internal.t_tuple_res;

    return result;
}

TtupleResult tTupleTest(const std::filesystem::path& filepath, SymbolMode mode) {
    std::ifstream file(filepath, std::ios::binary);
    if (!file) {
        throw std::runtime_error("tTupleTest: failed to open file " + filepath.string());
    }

    std::vector<uint8_t> buffer((std::istreambuf_iterator<char>(file)),
                                std::istreambuf_iterator<char>());

    EntropyInputData data{ buffer, 8 };
    return tTupleTest(data, mode);
}

} // namespace lib90b
