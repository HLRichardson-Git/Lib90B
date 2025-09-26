
#pragma once

#include <cstdint>

struct InternalLrsResult {
    // t-Tuple
    int t = -1;
    long double t_p_hat_max = 0.0;
    long double t_p_u = 0.0;
    long double t_tuple_res = -1.0;

    // LRS
    int u = -1;
    int v = -1;
    long double p_hat = 0.0;
    long double p_u = 0.0;
    long double lrs_res = -1.0;
};

InternalLrsResult SAalgs(const uint8_t text[], long n, int k);