
#pragma once

#include <string>
#include <filesystem>

#include "entropy_tests.h"

namespace lib90b {

LrsResult lrsTest(const EntropyInputData& data, SymbolMode mode = SymbolMode::Original);
LrsResult lrsTest(const std::filesystem::path& filepath, SymbolMode mode = SymbolMode::Original);

lenLrsResult lenLrsTest(const EntropyInputData& data);
lenLrsResult lenLrsTest(const std::filesystem::path& filepath);

} // namespace lib90b