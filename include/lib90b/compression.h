
#pragma once

#include <string>
#include <filesystem>

#include "entropy_tests.h"

namespace lib90b {

CompressionResult compressionTest(const EntropyInputData& data);
CompressionResult compressionTest(const std::filesystem::path& filepath);

} // namespace lib90b