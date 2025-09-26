
#pragma once

#include <string>
#include <filesystem>

#include "entropy_tests.h"

namespace lib90b {

LagResult lagTest(const EntropyInputData& data, SymbolMode mode = SymbolMode::Original);
LagResult lagTest(const std::filesystem::path& filepath, SymbolMode mode = SymbolMode::Original);

} // namespace lib90b
