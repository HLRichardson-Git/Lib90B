
#pragma once

#include <string>
#include <filesystem>

#include <lib90b/entropy_tests.h>

namespace lib90b {

CalcStatsResult calcStats(const EntropyInputData& data);
CalcStatsResult calcStats(const std::filesystem::path& filepath);

} // namespace lib90b
