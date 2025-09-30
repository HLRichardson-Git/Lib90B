#include <fstream>
#include <stdexcept>
#include <random>
#include <algorithm>
#include <iostream>
#include <vector>
#include <numeric>
#include <unordered_set>
#include <unordered_map>
#include <bzlib.h>

#include <lib90b/permutation.h>
#include <lib90b/calc_stats.h>
#include "utils.h"

// Number of permutations to run (should match PERMS constant from original)
#ifndef PERMS
#define PERMS 10000
#endif

/*
 * ---------------------------------------------
 * 	  TASKS FOR PERMUTATION TESTS
 * ---------------------------------------------
 */

// 5.1 Conversion I
static std::vector<uint8_t> conversion1(const uint8_t data[], const int sample_size){
	std::vector<uint8_t> ret((sample_size / 8) + ((sample_size%8==0)?0:1), 0);
	for(int i = 0; i < sample_size; ++i){
		ret[i/8] += data[i];
	}
	return ret;
}

// 5.1 Conversion II
static std::vector<uint8_t> conversion2(const uint8_t data[], const int sample_size){
	std::vector<uint8_t> ret((sample_size / 8) + ((sample_size%8==0)?0:1), 0);
	for(int i = 0; i < sample_size; ++i) {
		ret[i/8] += data[i] << (7 - i%8);
	}
	return ret;
}

// 5.1.1 Excursion Test
static double excursion(const uint8_t data[], const double rawmean, const int sample_size){
	double d_i = 0;
	double max = 0;
	double running_sum = 0;

	for(int i = 0; i < sample_size; ++i){
		running_sum += data[i];
		d_i = abs(running_sum - ((i+1) * rawmean));
		if(d_i > max){
			max = d_i;
		}
	}
	return max;
}

// Helper for 5.1.2, 5.1.3, and 5.1.4
static std::vector<int> alt_sequence1(const uint8_t data[], const int sample_size){
	std::vector<int> ret(sample_size-1, 0);
	for(int i = 0; i < sample_size-1; ++i){
		ret[i] = ((data[i] > data[i+1]) ? -1 : 1);
	}
	return ret;
}

// Helper for 5.1.5 and 5.1.6
static std::vector<int> alt_sequence2(const uint8_t data[], const double median, const int sample_size){
	std::vector<int> ret(sample_size, 0);
	for(int i = 0; i < sample_size; ++i){
		ret[i] = ((data[i] < median) ? -1 : 1);
	}
	return ret;
}

// 5.1.2 & 5.1.5 Number of Directional Runs
static unsigned int num_directional_runs(const std::vector<int> &alt_seq){
	unsigned int num_runs = 0;
	if(alt_seq.size() > 0) num_runs++;
	
	for(unsigned int i = 1; i < alt_seq.size(); ++i){
		if(alt_seq[i] != alt_seq[i-1]){
			++num_runs;
		}
	}
	return num_runs;
}

// 5.1.3 & 5.1.6 Length of Directional Runs
static unsigned int len_directional_runs(const std::vector<int> &alt_seq){
	unsigned int max_run = 0;
	unsigned int run = 1;

	for(unsigned int i = 1; i < alt_seq.size(); ++i){
		if(alt_seq[i] == alt_seq[i-1]){
			++run;
		}else{
			if(run > max_run){
				max_run = run;
			}
			run = 1;
		}
	}
	if(run > max_run){
		max_run = run;
	}
	return max_run;
}

// 5.1.4 Number of Increases and Decreases
static unsigned int num_increases_decreases(const std::vector<int> &alt_seq){
	unsigned int pos = 0;
	for(unsigned int i = 0; i < alt_seq.size(); ++i){
		if(alt_seq[i] == 1) ++pos;
	}
	unsigned int reverse_pos = alt_seq.size() - pos;
	return max(pos, reverse_pos);
}

// Helper function for collision tests
static std::vector<unsigned int> find_collisions(const uint8_t data[], const unsigned int n, const unsigned int k){
	std::vector<unsigned int> ret;
	std::vector<bool> dups(k, false);

	unsigned long int i=0;
	unsigned long int j=0;

	while(i + j < n){
		for(unsigned int l=0; l<k; l++) dups[l] = false;

		while(i + j < n) {
			if(dups[data[i+j]]) {
				ret.push_back(j+1);
				i += j;
				j=0;
				break;
			} else {
				dups[data[i+j]]=true;
				++j;
			}
		}
		++i;
	}
	return ret;
}

// 5.1.7 Average Collision Test
static double avg_collision(const std::vector<unsigned int> &col_seq){
	return divide(sum(col_seq), col_seq.size());
}

// 5.1.8 Maximum Collision Test
static unsigned int max_collision(const std::vector<unsigned int> &col_seq){
	unsigned int max = 0;
	for(unsigned int i = 0; i < col_seq.size(); ++i){
		if(max < col_seq[i]) max = col_seq[i];
	}
	return max;
}

// 5.1.9 Periodicity Test
static unsigned int periodicity(const uint8_t data[], const unsigned int p, const unsigned int n){
	unsigned int T = 0;
	assert(n>=p);
	for(unsigned int i = 0; i < n-p; ++i){
		if(data[i] == data[i+p]){
			++T;
		}
	}
	return T;
}

// 5.1.10 Covariance Test
static uint64_t covariance(const uint8_t data[], const unsigned int p, const unsigned int n) {
    uint64_t T = 0;
    for (unsigned int i = 0; i < n - p; ++i) {
        T += static_cast<uint64_t>(data[i]) * data[i + p];
    }
    return T;
}

// 5.1.11 Compression Test
static unsigned int compression(const uint8_t data[], const int sample_size, const uint8_t max_symbol){
	char *msg;
	unsigned int curlen = 0;
	char *curmsg;

	assert(max_symbol > 0);

	msg = new char[(size_t)(floor(log10(max_symbol))+2.0)*sample_size+1];
	msg[0] = '\0';
	curmsg = msg;

	for(int i = 0; i < sample_size; ++i) {
		int res;
		res = sprintf(curmsg, "%u ", data[i]);
		assert(res >= 2);
		curlen += res;
		curmsg += res;
	}

	if(curlen > 0) {
		assert(curmsg > msg);
		curmsg--;
		*curmsg = '\0';
		curlen--;
	}

	unsigned int dest_len = ceil(1.01*curlen) + 600;
	char* dest = new char[dest_len];

	int rc = BZ2_bzBuffToBuffCompress(dest, &dest_len, msg, curlen, 5, 0, 0);

	delete[](dest);
	delete[](msg);

	if(rc == BZ_OK){
		return dest_len;
	}else{
		return 0;
	}
}

/*
 * ---------------------------------------------
 * 	  HELPERS FOR PERMUTATION TEST ITERATION
 * ---------------------------------------------
 */

static void excursion_test(const uint8_t data[], const double rawmean, const int sample_size, 
                           long double* stats, const bool *test_status){
	if(test_status[0]) stats[0] = excursion(data, rawmean, sample_size);
}

static void directional_tests(const uint8_t data[], const int alphabet_size, const int sample_size, 
                              long double *stats, const bool *test_status){
	std::vector<int> alt_seq;

	if(test_status[1] || test_status[2] || test_status[3]) {
		if(alphabet_size == 2){
			std::vector<uint8_t> cs1 = conversion1(data, sample_size);
			alt_seq = alt_sequence1(cs1.data(), cs1.size());
		}else{
			alt_seq = alt_sequence1(data, sample_size);
		}

		if(test_status[1]) stats[1] = num_directional_runs(alt_seq);
		if(test_status[2]) stats[2] = len_directional_runs(alt_seq);
		if(test_status[3]) stats[3] = num_increases_decreases(alt_seq);
	}
}

static void consecutive_runs_tests(const uint8_t data[], const double median, const int alphabet_size, 
                                   const int sample_size, long double *stats, const bool *test_status){
	std::vector<int> alt_seq;

	if(test_status[4] || test_status[5]) {
		if(alphabet_size == 2){
			alt_seq = alt_sequence2(data, 0.5, sample_size);
		}else{
			alt_seq = alt_sequence2(data, median, sample_size);
		}

		if(test_status[4]) stats[4] = num_directional_runs(alt_seq);
		if(test_status[5]) stats[5] = len_directional_runs(alt_seq);
	}
}

static void collision_tests(const uint8_t data[], const int alphabet_size, const int sample_size, 
                           long double *stats, const bool *test_status){
	std::vector<unsigned int> col_seq;

	if(test_status[7] || test_status[6]) {
		if(alphabet_size == 2){
			std::vector<uint8_t> cs2 = conversion2(data, sample_size);
			col_seq = find_collisions(cs2.data(), cs2.size(), 256);
		}else{
			col_seq = find_collisions(data, sample_size, alphabet_size);
		}

		if(test_status[6]) stats[6] = avg_collision(col_seq);
		if(test_status[7]) stats[7] = max_collision(col_seq);
	}
}

static void periodicity_tests(const uint8_t data[], const int alphabet_size, const int sample_size,
                             long double *stats, const bool *test_status){
	if(test_status[8] || test_status[9] || test_status[10] || test_status[11] || test_status[12]) {
		if(alphabet_size == 2){
			std::vector<uint8_t> cs1 = conversion1(data, sample_size);
			if(test_status[8]) stats[8] = periodicity(cs1.data(), 1, cs1.size());
			if(test_status[9]) stats[9] = periodicity(cs1.data(), 2, cs1.size());
			if(test_status[10]) stats[10] = periodicity(cs1.data(), 8, cs1.size());
			if(test_status[11]) stats[11] = periodicity(cs1.data(), 16, cs1.size());
			if(test_status[12]) stats[12] = periodicity(cs1.data(), 32, cs1.size());
		}else{
			if(test_status[8]) stats[8] = periodicity(data, 1, sample_size);
			if(test_status[9]) stats[9] = periodicity(data, 2, sample_size);
			if(test_status[10]) stats[10] = periodicity(data, 8, sample_size);
			if(test_status[11]) stats[11] = periodicity(data, 16, sample_size);
			if(test_status[12]) stats[12] = periodicity(data, 32, sample_size);
		}
	}
}

static void covariance_tests(const uint8_t data[], const int alphabet_size, const int sample_size, 
                            long double *stats, const bool *test_status){
	if(test_status[13] || test_status[14] || test_status[15] || test_status[16] || test_status[17]) {
		if(alphabet_size == 2){
			std::vector<uint8_t> cs1 = conversion1(data, sample_size);
			if(test_status[13]) stats[13] = covariance(cs1.data(), 1, cs1.size());
			if(test_status[14]) stats[14] = covariance(cs1.data(), 2, cs1.size());
			if(test_status[15]) stats[15] = covariance(cs1.data(), 8, cs1.size());
			if(test_status[16]) stats[16] = covariance(cs1.data(), 16, cs1.size());
			if(test_status[17]) stats[17] = covariance(cs1.data(), 32, cs1.size());
		}else{
			if(test_status[13]) stats[13] = covariance(data, 1, sample_size);
			if(test_status[14]) stats[14] = covariance(data, 2, sample_size);
			if(test_status[15]) stats[15] = covariance(data, 8, sample_size);
			if(test_status[16]) stats[16] = covariance(data, 16, sample_size);
			if(test_status[17]) stats[17] = covariance(data, 32, sample_size);
		}
	}
}

static void compression_test(const uint8_t data[], const int sample_size, long double *stats, 
                            const uint8_t max_symbol, const bool *test_status){
	if(test_status[18]) stats[18] = compression(data, sample_size, max_symbol);
}

static void run_tests(const uint8_t remapped_data[], const uint8_t raw_data[], 
                     const double rawmean, const double median, 
                     int alphabet_size, int sample_size, uint8_t max_symbol,
                     long double *stats, const bool *test_status){
	
	excursion_test(raw_data, rawmean, sample_size, stats, test_status);
	directional_tests(remapped_data, alphabet_size, sample_size, stats, test_status);
	consecutive_runs_tests(remapped_data, median, alphabet_size, sample_size, stats, test_status);
	collision_tests(remapped_data, alphabet_size, sample_size, stats, test_status);
	periodicity_tests(remapped_data, alphabet_size, sample_size, stats, test_status);
	
	if(alphabet_size == 2) {
		covariance_tests(remapped_data, alphabet_size, sample_size, stats, test_status);
	} else {
		covariance_tests(raw_data, alphabet_size, sample_size, stats, test_status);
	}
	compression_test(raw_data, sample_size, stats, max_symbol, test_status);
}

namespace lib90b {

PermutationTestResult permutationTest(const EntropyInputData& input, std::optional<uint64_t> fixedSeed)
{
    const auto& raw_symbols = input.symbols;
    const int sample_size = static_cast<int>(raw_symbols.size());

    // Determine alphabet size
    int alphabet_size = input.alph_size > 0 ? input.alph_size
                                            : static_cast<int>(std::unordered_set<uint8_t>(raw_symbols.begin(), raw_symbols.end()).size());

    // Map symbols
    std::unordered_set<uint8_t> unique_symbols(raw_symbols.begin(), raw_symbols.end());
    std::vector<uint8_t> sorted_symbols(unique_symbols.begin(), unique_symbols.end());
    std::sort(sorted_symbols.begin(), sorted_symbols.end());

    std::unordered_map<uint8_t, int> symbol_map;
    for (size_t i = 0; i < sorted_symbols.size(); i++)
        symbol_map[sorted_symbols[i]] = static_cast<int>(i);

    std::vector<uint8_t> remapped_data(sample_size);
    for (int i = 0; i < sample_size; i++)
        remapped_data[i] = static_cast<uint8_t>(symbol_map.at(raw_symbols[i]));

    // Calculate initial statistics
    double rawmean = std::accumulate(raw_symbols.begin(), raw_symbols.end(), 0.0) / sample_size;
    std::vector<uint8_t> sorted_remapped = remapped_data;
    std::sort(sorted_remapped.begin(), sorted_remapped.end());
    double median = (sample_size % 2 == 0)
                        ? (sorted_remapped[sample_size/2 - 1] + sorted_remapped[sample_size/2]) / 2.0
                        : sorted_remapped[sample_size/2];

    uint8_t max_symbol = *std::max_element(raw_symbols.begin(), raw_symbols.end());

    // Initialize RNG
    uint64_t seed = fixedSeed.value_or(std::random_device{}());
    std::mt19937_64 rng(seed);

    PermutationTestResult result{};
    result.initialStats.resize(19, 0);
    result.counts.resize(19, {0, 0, 0});

    bool test_status[19];
    std::fill_n(test_status, 19, true);

    run_tests(remapped_data.data(), raw_symbols.data(), rawmean, median,
              alphabet_size, sample_size, max_symbol,
              result.initialStats.data(), test_status);

    // Permutation loop
    std::vector<uint8_t> perm_remapped = remapped_data;
    std::vector<uint8_t> perm_raw = raw_symbols;
    int passed_count = 0;

    for (int p = 0; p < PERMS; ++p) {
        if (passed_count < 19) {
            std::vector<size_t> indices(sample_size);
            std::iota(indices.begin(), indices.end(), 0);
            std::shuffle(indices.begin(), indices.end(), rng);

            for (int i = 0; i < sample_size; i++) {
                perm_remapped[i] = remapped_data[indices[i]];
                perm_raw[i] = raw_symbols[indices[i]];
            }

            long double permuted_stats[19] = {0};
            run_tests(perm_remapped.data(), perm_raw.data(), rawmean, median,
                      alphabet_size, sample_size, max_symbol,
                      permuted_stats, test_status);

            for (int j = 0; j < 19; ++j) {
                if (test_status[j]) {
                    if (permuted_stats[j] > result.initialStats[j]) result.counts[j][0]++;
                    else if (permuted_stats[j] == result.initialStats[j]) result.counts[j][1]++;
                    else result.counts[j][2]++;

                    if ((result.counts[j][0] + result.counts[j][1] > 5) &&
                        (result.counts[j][1] + result.counts[j][2] > 5))
                        test_status[j] = false;
                }
            }

            passed_count = std::count(std::begin(test_status), std::end(test_status), false);
        }
    }

    result.passed = true;
    for (int i = 0; i < 19; ++i) {
        if ((result.counts[i][0] + result.counts[i][1] <= 5) ||
            (result.counts[i][1] + result.counts[i][2] <= 5)) {
            result.passed = false;
            break;
        }
    }

    return result;
}

PermutationTestResult permutationTest(const std::filesystem::path& filepath, std::optional<uint64_t> fixedSeed) {
    std::ifstream file(filepath, std::ios::binary);
    if (!file) {
        throw std::runtime_error("Failed to open file: " + filepath.string());
    }

    std::vector<uint8_t> buffer(
        (std::istreambuf_iterator<char>(file)),
        std::istreambuf_iterator<char>()
    );

    EntropyInputData data{buffer, 8};
    return permutationTest(data, fixedSeed);
}

} // namespace lib90b