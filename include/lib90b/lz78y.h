
#pragma once

#include <string>
#include <filesystem>

#include "entropy_tests.h"

namespace lib90b {

Lz78yResult lz78yTest(const EntropyInputData& data, SymbolMode mode = SymbolMode::Original);
Lz78yResult lz78yTest(const std::filesystem::path& filepath, SymbolMode mode = SymbolMode::Original);

} // namespace lib90b
