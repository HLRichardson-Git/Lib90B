
#pragma once

#include <string>
#include <filesystem>

#include "entropy_tests.h"

namespace lib90b {

MultiMmcResult multiMmcTest(const EntropyInputData& data, SymbolMode mode = SymbolMode::Original);
MultiMmcResult multiMmcTest(const std::filesystem::path& filepath, SymbolMode mode = SymbolMode::Original);

} // namespace lib90b
