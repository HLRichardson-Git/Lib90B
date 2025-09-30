
#pragma once

#include <filesystem>
#include <optional>
#include <vector>

#include <lib90b/entropy_tests.h>

namespace lib90b {

IidResult iidTestSuite(const EntropyInputData& data);
IidResult iidTestSuite(const std::filesystem::path& filepath);

} // namespace lib90b
