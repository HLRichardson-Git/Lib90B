
#pragma once

#include <string>
#include <filesystem>

#include "entropy_tests.h"

namespace lib90b {

MarkovResult markovTest(const EntropyInputData& data);
MarkovResult markovTest(const std::filesystem::path& filepath);

} // namespace lib90b