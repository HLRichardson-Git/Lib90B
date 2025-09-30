
#pragma once

#include <vector>
#include <climits>

#include <divsufsort.h>
#include <divsufsort64.h>

#include "lrs_algs.h"
#include "utils.h"

#define SAINDEX_MAX INT32_MAX
#define SAINDEX64_MAX INT64_MAX

//Using the Kasai (et al.) O(n) time "13n space" algorithm.
//"Linear-Time Longest-Common-Prefix Computation in Suffix Arrays and Its Applications", by Kasai, Lee, Arimura, Arikawa, and Park
//https://doi.org/10.1007/3-540-48194-X_17
//http://web.cs.iastate.edu/~cs548/references/linear_lcp.pdf
//The default implementation uses 4 byte indexes
static void sa2lcp32(const uint8_t text[], long int n, const std::vector<saidx_t> &sa, std::vector<saidx_t> &lcp) {
	saidx_t h;
	std::vector<saidx_t> rank(n+1,-1);

	assert(n>1);

	lcp[0] = -1;
	lcp[1] = 0;

	// compute rank = sa^{-1}
	for(saidx_t i=0; i<=(saidx_t)n; i++) {
		rank[sa[i]] = i;
	}

	// traverse suffixes in rank order
	h=0;

	for(saidx_t i=0; i<(saidx_t)n; i++) {
		saidx_t k = rank[i]; // rank of s[i ... n-1]
		if(k>1) {
			saidx_t j = sa[k-1]; // predecessor of s[i ... n-1]
			while((i+h<(saidx_t)n) && (j+h<(saidx_t)n) && (text[i+h]==text[j+h])) {
				h++;
			}

			lcp[k] = h;
		}
		if(h>0) {
			h--;
		}
	}
}

//Using the Kasai (et al.) O(n) time "25n space" algorithm (with 64-bit indicies)
static void sa2lcp64(const uint8_t text[], long int n, const std::vector<saidx64_t> &sa, std::vector<saidx64_t> &lcp) {
	saidx64_t h;
	std::vector<saidx64_t> rank(n+1,-1);

	assert(n>1);

	lcp[0] = -1;
	lcp[1] = 0;

	// compute rank = sa^{-1}
	for(saidx64_t i=0; i<=(saidx64_t)n; i++) {
		rank[sa[i]] = i;
	}

	// traverse suffixes in rank order
	h=0;

	for(saidx64_t i=0; i<(saidx64_t)n; i++) {
		saidx64_t k = rank[i]; // rank of s[i ... n-1]
		if(k>1) {
			saidx64_t j = sa[k-1]; // predecessor of s[i ... n-1]
			while((i+h<(saidx64_t)n) && (j+h<(saidx64_t)n) && (text[i+h]==text[j+h])) {
				h++;
			}

			lcp[k] = h;
		}
		if(h>0) {
			h--;
		}
	}
}


static void calcSALCP32(const uint8_t text[], long int n, std::vector<saidx_t> &sa, std::vector<saidx_t> &lcp) {
	int32_t res;

	assert(n < SAINDEX_MAX);
	assert(n > 0); 
	assert(sa.size() == (size_t)(n+1));
	assert(sa.size() == (size_t)(n+1));

	sa[0] = (saidx_t)n;

	res=divsufsort((const sauchar_t *)text, (saidx_t *)(sa.data()+1), (saidx_t)n);
	assert(res==0);
   	sa2lcp32(text, n, sa, lcp);
}

static void calcSALCP64(const uint8_t text[], long int n, std::vector<saidx64_t> &sa, std::vector<saidx64_t> &lcp) {
	int32_t res;

	assert(n < SAINDEX64_MAX);
	assert(n > 0);
	assert(sa.size() == (size_t)(n+1));
	assert(sa.size() == (size_t)(n+1));

	sa[0] = (saidx64_t)n;

	res=divsufsort64((const sauchar_t *)text, (saidx64_t *)(sa.data()+1), (saidx64_t)n);
	assert(res==0);
   	sa2lcp64(text, n, sa, lcp);
}
/* Based on the algorithm outlined by Aaron Kaufer
 * This is described here:
 * http://www.untruth.org/~josh/sp80090b/Kaufer%20Further%20Improvements%20for%20SP%20800-90B%20Tuple%20Counts.pdf
 */
static InternalLrsResult SAalgs32(const uint8_t text[], long n, int k) {
    InternalLrsResult result;

    std::vector<saidx_t> sa(n+1, -1); 
    std::vector<saidx_t> L(n+2, -1); 

    long int u; 
    long int v; 
    long int c; 
    long int j; 
    saidx_t t;   

    long double Pmax;
    long double pu;

    assert(n > 0);
    assert(k > 0);
    assert(n < SAINDEX_MAX);
    assert((UINT64_MAX / (uint64_t)n) >= ((uint64_t)n+1U));

    calcSALCP32(text, n, sa, L);

    // adjust indexing to conform with Kaufer
    L.erase(L.begin());
    L[n] = 0;
    assert(L[0] == 0);

    // Find v (longest repeated substring length)
    v = 0;
    for (long int i = 0; i < n; i++) {
        if (L[i] > v) v = L[i];
    }
    assert((v > 0) && (v < n));

    std::vector<saidx_t> Q(v+1, 1);
    std::vector<saidx_t> A(v+2, 0);
    std::vector<saidx_t> I(v+3, 0);

    j = 0;
    for (long int i = 1; i <= n; i++) {
        c = 0;
        assert(L[i] >= 0);

        if (L[i] < L[i-1]) {
            t = L[i-1];
            assert(j > 0);
            j--;
            assert(j <= v);

            while (t > L[i]) {
                assert((t > 0) && (t <= v));
                if ((j > 0) && (I[j] == t)) {
                    A[I[j]] += A[I[j+1]];
                    A[I[j+1]] = 0;
                    j--;
                }

                if (Q[t] >= A[I[j+1]]+1) {
                    if (j > 0) {
                        t = I[j];
                    } else {
                        t = L[i];
                    }
                } else {
                    Q[t--] = A[I[j+1]]+1;
                }
            }
            c = A[I[j+1]];
            A[I[j+1]] = 0;
        }

        if (L[i] > 0) {
            if ((j < 1) || (I[j] < L[i])) {
                assert(j < v);
                I[++j] = L[i];
            }
            A[I[j]] += c+1;
        }
    }

    // Calculate u
    for (u = 1; (u <= v) && (Q[u] >= 35); u++);
    assert(u > 0);
    assert(((u == v+1) || ((u <= v) && (Q[u] < 35))));
	assert(((u == 1) || (Q[u-1] >= 35)));

    result.u = (int)u;
    result.v = (int)v;

    // --- t-Tuple estimate ---
    Pmax = -1.0;
    for (long int i = 1; i < u; i++) {
        long double curP = ((long double)(Q[i])) / ((long double)(n-i+1));
        long double curPMax = powl(curP, 1.0L/(long double)i);
        if (curPMax > Pmax) Pmax = curPMax;
    }

    if (Pmax > 0.0L) {
        pu = Pmax + ZALPHA_L * sqrtl(Pmax * (1.0L - Pmax) / ((long double)(n - 1)));
        if (pu > 1.0L) pu = 1.0L;

        result.t = (int)(u-1);
        result.t_p_hat_max = Pmax;
        result.t_p_u = pu;
        result.t_tuple_res = (double)-log2l(pu);
    }

    // --- LRS estimate ---
    if (v >= u) {
        std::vector<uint64_t> S(v+1, 0);
        memset(A.data(), 0, sizeof(saidx_t) * ((size_t)v+2));

        for (long int i = 1; i <= n; i++) {
            if ((L[i-1] >= u) && (L[i] < L[i-1])) {
                saidx_t b = L[i];
                if (b < u) b = u-1;

                for (t = L[i-1]; t > b; t--) {
                    uint64_t priorS = S[t];
                    uint64_t choices = ((((uint64_t)(A[t]+1) * (uint64_t)(A[t])))) >> 1;
                    A[t] += A[t+1];
                    A[t+1] = 0;
                    assert(A[t] >= 0);

                    priorS = S[t];
					choices = ((((uint64_t)(A[t]+1) * (uint64_t)(A[t]))))>>1;
                    S[t] = priorS + choices;
                    assert(S[t] >= priorS);
                }
                if (b >= u) A[b] += A[b+1];
                A[b+1] = 0;
            }
            if (L[i] >= u) A[L[i]]++;
        }

        Pmax = 0.0;
        for (long int i = u; i <= v; i++) {
            uint64_t choices = (((uint64_t)n-(uint64_t)i)*((uint64_t)n-(uint64_t)i+1U)) >> 1;
            long double curP = ((long double)S[i]) / (long double)choices;
            long double curPMax = pow(curP, 1.0/((long double)i));
            if (Pmax < curPMax) Pmax = curPMax;
        }

        pu = Pmax + ZALPHA_L * sqrtl(Pmax * (1.0L - Pmax) / ((long double)(n - 1)));
        if (pu > 1.0L) pu = 1.0L;

        result.p_hat = Pmax;
        result.p_u = pu;
        result.lrs_res = (double)-log2l(pu);
    } else {
        result.lrs_res = -1.0;
    }

    return result;
}

static InternalLrsResult SAalgs64(const uint8_t text[], long n, int k) {
    InternalLrsResult result;

    std::vector<saidx64_t> sa(n+1, -1);
    std::vector<saidx64_t> L(n+2, -1);

    long int u;
    long int v;
    long int c;
    long int j;
    saidx64_t t;

    long double Pmax;
    long double pu;

    assert(n > 0);
    assert(k > 0);
    assert(n <= SAINDEX64_MAX - 1);
    assert((UINT128_MAX / (uint128_t)n) >= ((uint128_t)n+1U));

    calcSALCP64(text, n, sa, L);

    // Adjust indexing to conform with Kaufer
    L.erase(L.begin());
    L[n] = 0;
    assert(L[0] == 0);

    // Find v (length of longest repeated substring)
    v = 0;
    for (long int i = 0; i < n; i++) {
        if (L[i] > v) v = L[i];
    }
    assert((v > 0) && (v < n));

    std::vector<saidx64_t> Q(v+1, 1);
    std::vector<saidx64_t> A(v+2, 0);
    std::vector<saidx64_t> I(v+3, 0);

    j = 0;
    for (long int i = 1; i <= n; i++) {
        c = 0;
        assert(L[i] >= 0);

        if (L[i] < L[i-1]) {
            t = L[i-1];
            assert(j > 0);
            j--;
            assert(j <= v);

            while (t > L[i]) {
                assert((t > 0) && (t <= v));
                if ((j > 0) && (I[j] == t)) {
                    A[I[j]] += A[I[j+1]];
                    A[I[j+1]] = 0;
                    j--;
                }

                if (Q[t] >= A[I[j+1]]+1) {
                    if (j > 0) {
                        t = I[j];
                    } else {
                        t = L[i];
                    }
                } else {
                    Q[t--] = A[I[j+1]]+1;
                }
            }
            c = A[I[j+1]];
            A[I[j+1]] = 0;
        }

        if (L[i] > 0) {
            if ((j < 1) || (I[j] < L[i])) {
                assert(j < v);
                I[++j] = L[i];
            }
            A[I[j]] += c+1;
        }
    }

    // Calculate u
    for (u = 1; (u <= v) && (Q[u] >= 35); u++);
    assert(u > 0);
    assert(((u == v+1) || ((u <= v) && (Q[u] < 35))));
	assert(((u == 1) || (Q[u-1] >= 35)));

    result.u = (int)u;
    result.v = (int)v;

    // --- t-Tuple estimate ---
    Pmax = -1.0;
    for (long int i = 1; i < u; i++) {
        long double curP = ((long double)(Q[i])) / ((long double)(n-i+1));
        long double curPMax = powl(curP, 1.0L/(long double)i);
        if (curPMax > Pmax) Pmax = curPMax;
    }

    if (Pmax > 0.0L) {
        pu = Pmax + ZALPHA_L * sqrtl(Pmax * (1.0L - Pmax) / ((long double)(n - 1)));
        if (pu > 1.0L) pu = 1.0L;

        result.t = (int)(u-1);
        result.t_p_hat_max = Pmax;
        result.t_p_u = pu;
        result.t_tuple_res = (double)-log2l(pu);
    }

    // --- LRS estimate ---
    if (v >= u) {
        std::vector<uint128_t> S(v+1, 0);
        memset(A.data(), 0, sizeof(saidx64_t) * ((size_t)v+2));

        for (long int i = 1; i <= n; i++) {
            if ((L[i-1] >= u) && (L[i] < L[i-1])) {
                saidx64_t b = L[i];
                if (b < u) b = u-1;

                for (t = L[i-1]; t > b; t--) {
                    uint128_t priorS = S[t];
                    uint128_t choices = ((((uint128_t)(A[t]+1) * (uint128_t)(A[t])))) >> 1;
                    A[t] += A[t+1];
                    A[t+1] = 0;
                    assert(A[t] >= 0);

                    priorS = S[t];
					choices = ((((uint128_t)(A[t]+1) * (uint128_t)(A[t]))))>>1;
                    S[t] = priorS + choices;
                    assert(S[t] >= priorS);
                }
                if (b >= u) A[b] += A[b+1];
                A[b+1] = 0;
            }
            if (L[i] >= u) A[L[i]]++;
        }

        Pmax = 0.0;
        for (long int i = u; i <= v; i++) {
            uint128_t choices = (((uint128_t)n-(uint128_t)i) * ((uint128_t)n-(uint128_t)i+1U)) >> 1;
            long double curP = ((long double)S[i]) / (long double)choices;
            long double curPMax = powl(curP, 1.0L/((long double)i));
            if (Pmax < curPMax) Pmax = curPMax;
        }

        pu = Pmax + ZALPHA_L * sqrtl(Pmax * (1.0L - Pmax) / ((long double)(n - 1)));
        if (pu > 1.0L) pu = 1.0L;

        result.p_hat = Pmax;
        result.p_u = pu;
        result.lrs_res = (double)-log2l(pu);
    } else {
        result.lrs_res = -1.0;
    }

    return result;
}

InternalLrsResult SAalgs(const uint8_t text[], long n, int k) {
	if(n<SAINDEX_MAX) {
        return SAalgs32(text, n, k);
	} else {
        return SAalgs64(text, n, k);
	}
}

// Helper function to calculate collision proportion
void calc_collision_proportion(const std::vector<double> &p, long double &p_col) {
    p_col = 0.0L;
    for (size_t i = 0; i < p.size(); i++) {
        p_col += powl((long double)(p[i]), 2.0L);
    }
}

// Helper to get the length of the LRS
long int len_LRS32(const uint8_t text[], int sample_size) {
    std::vector<saidx_t> sa(sample_size + 1, -1);
    std::vector<saidx_t> lcp(sample_size + 1, -1);
    saidx_t lrs_len = -1;

    calcSALCP32(text, sample_size, sa, lcp);

    for (saidx_t j = 0; j <= sample_size; j++) {
        if (lcp[j] > lrs_len) lrs_len = lcp[j];
    }

    return lrs_len;
}

long int len_LRS64(const uint8_t text[], int sample_size) {
    std::vector<saidx64_t> sa(sample_size + 1, -1);
    std::vector<saidx64_t> lcp(sample_size + 1, -1);
    saidx64_t lrs_len = -1;

    calcSALCP64(text, sample_size, sa, lcp);

    for (saidx64_t j = 0; j <= sample_size; j++) {
        if (lcp[j] > lrs_len) lrs_len = lcp[j];
    }

    return lrs_len;
}