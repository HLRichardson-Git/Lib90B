
#pragma once

#include <string>
#include <filesystem>

#include <lib90b/entropy_tests.h>

namespace lib90b {

TtupleResult tTupleTest(const EntropyInputData& data, SymbolMode mode = SymbolMode::Original);
TtupleResult tTupleTest(const std::filesystem::path& filepath, SymbolMode mode = SymbolMode::Original);

} // namespace lib90b