
#pragma once

#include <vector>
#include <filesystem>
#include <fstream>

#include <lib90b/compression.h>
#include "utils.h"

static void kahan_add(double &sum, double &comp, double in){
	double y, t; 

	y = in - comp;
	t = sum + y;
	comp = (t - sum) - y;
	sum = t;
}

//There is some cleverness associated with this calculation of G; in particular,
//one doesn't need to calculate all the terms independently (they are inter-related!)
//See UL's implementation comments here: https://bit.ly/UL90BCOM 
//Look in the section "Compression Estimate G Function Calculation"
static double G(double z, int d, long num_blocks){
	double Ai=0.0, Ai_comp=0.0;
	double firstSum=0.0, firstSum_comp=0.0;
	long v = num_blocks - d;
	double Ad1;

	long double Bi;
	long double Bterm;
	long double ai;
	long double aiScaled;
	bool underflowTruncate;

	assert(d>0);
	assert(num_blocks>d);

	//i=2
	Bterm = (1.0L-(long double)z);
	//Note: B_1 isn't needed, as a_1 = 0
	//B_2
	Bi = Bterm;

	//Calculate A_{d+1}
	for(int i=2; i<=d; i++) {
		//calculate the a_i term
		kahan_add(Ai, Ai_comp, log2l((long double)i)*Bi);

		//Calculate B_{i+1}
		Bi *= Bterm;
	}

	//Store A_{d+1}
	Ad1 = Ai;

	underflowTruncate = false;
	//Now calculate A_{num_blocks} and the sum of sums term (firstsum)
	for(long i=d+1; i<=num_blocks-1; i++) {
		//calculate the a_i term
		ai = log2l((long double)i)*Bi;

		//Calculate A_{i+1}
		kahan_add(Ai, Ai_comp, (double)ai);
		//Sum in A_{i+1} into the firstSum

		//Calculate the tail of the sum of sums term (firstsum)
		aiScaled = (long double)(num_blocks-i) * ai;
		if((double)aiScaled > 0.0) {
			kahan_add(firstSum, firstSum_comp, (double)aiScaled);
		} else {
			underflowTruncate = true;
			break;
		}

		//Calculate B_{i+1}
		Bi *= Bterm;
	}

	//Ai now contains A_{num_blocks} and firstsum contains the tail
	//finalize the calculation of firstsum
	kahan_add(firstSum, firstSum_comp, ((double)(num_blocks-d))*Ad1);

	//Calculate A_{num_blocks+1}
	if(!underflowTruncate) {
		ai = log2l((long double)num_blocks)*Bi;
		kahan_add(Ai, Ai_comp, (double)ai);
	}

	return 1/(double)v * z*(z*firstSum + (Ai - Ad1));
}

static double com_exp(double p, unsigned int alph_size, int d, long num_blocks){
	double q = (1.0-p)/((double)alph_size-1.0);
        return G(p, d, num_blocks) + ((double)alph_size-1.0) * G(q, d, num_blocks);
}

namespace lib90b {

// Section 6.3.4 - Compression Estimate
CompressionResult compressionTest(const EntropyInputData& data) {
    const std::vector<uint8_t>& symbols_to_use = data.getBits();
    const long len = static_cast<long>(symbols_to_use.size());

    const int b = 6;
    const unsigned int alph_size = 1 << b;
    const int d = 1000;
    const long num_blocks = len / b;

    if (num_blocks <= d) {
        throw std::runtime_error("compressionTest: insufficient data length (need > 1000 blocks)");
    }

    std::vector<unsigned int> dict(alph_size, 0);
    double X = 0.0, X_comp = 0.0;
    double sigma = 0.0, sigma_comp = 0.0;

    // Initialize dictionary
    for (long i = 0; i < d; i++) {
        unsigned int block = 0;
        for (int j = 0; j < b; j++) {
            block |= (symbols_to_use[i * b + j] & 0x1) << (b - j - 1);
        }
        dict[block] = i + 1;
    }

    long v = num_blocks - d;

    // Test data against dictionary with long double for log
    for (long i = d; i < num_blocks; i++) {
        unsigned int block = 0;
        for (int j = 0; j < b; j++) {
            block |= (symbols_to_use[i * b + j] & 0x1) << (b - j - 1);
        }

        long double log_val = log2l((long double)(i + 1 - dict[block]));
        kahan_add(X, X_comp, (double)log_val);
        kahan_add(sigma, sigma_comp, (double)(log_val * log_val));
        dict[block] = i + 1;
    }

    // Compute mean and std
    X /= v;
    sigma = 0.5907 * std::sqrt(sigma / (v - 1.0) - X * X);

    // Adjust mean for confidence interval
    double X_prime = X - ZALPHA * sigma / std::sqrt(v);

    double p = -1.0;
    bool found_p = false;

    if (com_exp(1.0 / (double)alph_size, alph_size, d, num_blocks) > X_prime) {
        double ldomain = 1.0 / (double)alph_size;
        double hdomain = 1.0;
        double lbound = ldomain, hbound = hdomain;
        double lvalue = DBL_INFINITY, hvalue = -DBL_INFINITY;
        double lastP;

        p = (lbound + hbound) / 2.0;
        double pVal = com_exp(p, alph_size, d, num_blocks);

        for (int iter = 0; iter < ITERMAX; iter++) {
            if (relEpsilonEqual(pVal, X_prime, ABSEPSILON, RELEPSILON, 4)) break;

            if (X_prime < pVal) {
                lbound = p;
                lvalue = pVal;
            } else {
                hbound = p;
                hvalue = pVal;
            }

            if (lbound >= hbound ||
                !INCLOSEDINTERVAL(X_prime, lvalue, hvalue)) {
                p = lbound;
                break;
            }

            lastP = p;
            p = (lbound + hbound) / 2.0;

            if (!INOPENINTERVAL(p, lbound, hbound)) {
                p = hbound;
                break;
            }

            if (lastP == p) {
                p = hbound;
                break;
            }

            pVal = com_exp(p, alph_size, d, num_blocks);

            if (!INCLOSEDINTERVAL(pVal, lvalue, hvalue)) {
                p = hbound;
                break;
            }
        }
    }

    double entEst = 1.0;
    if (p > 1.0 / (double)alph_size) {
        entEst = -std::log2(p) / b;
        found_p = true;
    } else {
        p = 1.0 / (double)alph_size;
        entEst = 1.0;
    }
    
    CompressionResult result;
    result.x_bar = X;
    result.sigma_hat = sigma;
    result.x_bar_prime = X_prime;
    result.p = p;
    result.found_p = found_p;
    result.h_bitstring = entEst;

    return result;
}

CompressionResult compressionTest(const std::filesystem::path& filepath) {
    std::ifstream file(filepath, std::ios::binary);
    if (!file) throw std::runtime_error("compressionTest: failed to open file " + filepath.string());

    std::vector<uint8_t> buffer((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
    EntropyInputData data{ buffer, 8 }; // default word_size = 8
    return compressionTest(data);
}

} // namespace lib90b
