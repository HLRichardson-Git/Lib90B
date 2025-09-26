
#pragma once

#include <string>
#include <filesystem>

#include "entropy_tests.h"

namespace lib90b {

CollisionResult collisionTest(const EntropyInputData& data);
CollisionResult collisionTest(const std::filesystem::path& filepath);

} // namespace lib90b
