
#pragma once

#include <string>
#include <filesystem>

#include <lib90b/entropy_tests.h>

namespace lib90b {

struct NonIidResult {
    std::filesystem::path filename;
    std::string sha256;
    int word_size = 0;
    bool is_initial_entropy = true;

    std::optional<double> min_entropy;
    std::optional<double> H_original;
    std::optional<double> H_bitstring;

    MostCommonResult most_common;
    CollisionResult collision;
    MarkovResult markov;
    CompressionResult compression;
    TtupleResult t_tuple;
    LrsResult lrs;
    MultiMcwResult multi_mcw;
    LagResult lag;
    MultiMmcResult multi_mmc;
    Lz78yResult lz78y;
};

NonIidResult nonIidTestSuite(const std::filesystem::path& filepath);

} // namespace lib90b
