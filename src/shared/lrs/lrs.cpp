#include <vector>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <algorithm> // for std::min
#include <stdexcept>

#include <lib90b/lrs.h>
#include "shared/lrs/lrs_algs.h"
#include <lib90b/entropy_tests.h>
#include "utils.h"

namespace lib90b {

LrsResult lrsTest(const EntropyInputData& data, SymbolMode mode) {
    const auto& symbols = (mode == SymbolMode::Bitstring) ? data.getBits() : data.symbols;
    InternalLrsResult internal = SAalgs(symbols.data(), symbols.size(),
                                        (mode == SymbolMode::Bitstring) ? 2 : (1 << data.word_size));

    LrsResult result;
    result.u = internal.u;
    result.v = internal.v;
    result.p_hat = internal.p_hat;
    result.p_u = internal.p_u;
    result.lrs_res = internal.lrs_res;

    if (mode == SymbolMode::Original) result.h_original = internal.lrs_res;
    else result.h_bitstring = internal.lrs_res;

    return result;
}


LrsResult lrsTest(const std::filesystem::path& filepath, SymbolMode mode) {
    std::ifstream file(filepath, std::ios::binary);
    if (!file) {
        throw std::runtime_error("lrsTest: failed to open file " + filepath.string());
    }

    std::vector<uint8_t> buffer(
        (std::istreambuf_iterator<char>(file)),
        std::istreambuf_iterator<char>());

    EntropyInputData data{ buffer, 8 };
    return lrsTest(data, mode);
}

} // namespace lib90b