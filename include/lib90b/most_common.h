
#pragma once

#include <string>
#include <filesystem>

#include "entropy_tests.h"

namespace lib90b {

MostCommonResult mostCommonTest(const EntropyInputData& data, SymbolMode mode = SymbolMode::Original);
MostCommonResult mostCommonTest(const std::filesystem::path& filepath, SymbolMode mode = SymbolMode::Original);

} // namespace lib90b
