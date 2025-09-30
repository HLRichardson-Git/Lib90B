
#pragma once

#include <string>
#include <filesystem>

#include <lib90b/entropy_tests.h>

namespace lib90b {

PermutationTestResult permutationTest(const EntropyInputData& data, std::optional<uint64_t> fixedSeed = std::nullopt);
PermutationTestResult permutationTest(const std::filesystem::path& filepath, std::optional<uint64_t> fixedSeed = std::nullopt);

} // namespace lib90b