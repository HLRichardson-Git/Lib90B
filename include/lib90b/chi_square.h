
#pragma once

#include <string>
#include <filesystem>

#include <lib90b/entropy_tests.h>

namespace lib90b {

ChiSquareResult chiSquareTest(const EntropyInputData& data);
ChiSquareResult chiSquareTest(const std::filesystem::path& filepath);

} // namespace lib90b