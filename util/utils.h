
#pragma once

#include <iostream>
#include <string>
#include <map>
#include <set>
#include <iomanip>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <ctime>
#include <algorithm>
#include <cmath>
#include <array>
#include <omp.h>
#include <bitset>
#include <mutex>
#include <cassert>
#include <cfloat>
#include <cstdint>
#include <stdexcept>
#include <limits>

// Version
#define VERSION "1.1.8"

// Macros
#define SWAP(x, y) do { int s = x; x = y; y = s; } while(0)
#define INOPENINTERVAL(x, a, b) (((a)>(b))?(((x)>(b))&&((x)<(a))):(((x)>(a))&&((x)<(b))))
#define INCLOSEDINTERVAL(x, a, b) (((a)>(b))?(((x)>=(b))&&((x)<=(a))):(((x)>=(a))&&((x)<=(b))))

#define MIN_SIZE 1000000
#define PERMS 10000

#define RELEPSILON DBL_EPSILON
#define ABSEPSILON DBL_MIN
#define DBL_INFINITY std::numeric_limits<double>::infinity()
#define ITERMAX 1076
#define ZALPHA 2.5758293035489008
#define ZALPHA_L 2.575829303548900384158L

// Forward declarations
struct data_t {
    int word_size;
    int alph_size;
    uint8_t maxsymbol;
    uint8_t *rawsymbols;
    uint8_t *symbols;
    uint8_t *bsymbols;
    long len;
    long blen;
};

struct uint128_t {
    uint64_t low;
    uint64_t high;

    // Constructors
    uint128_t() : low(0), high(0) {}
    uint128_t(uint64_t v) : low(v), high(0) {}
    uint128_t(uint64_t low_, uint64_t high_) : low(low_), high(high_) {}

    // Conversions
    explicit operator uint64_t() const;
    explicit operator long double() const;
    explicit operator double() const;

    // Arithmetic
    uint128_t operator+(const uint128_t& other) const;
    uint128_t& operator+=(const uint128_t& other);
    uint128_t operator-(const uint128_t& other) const;
    uint128_t& operator-=(const uint128_t& other);
    uint128_t operator*(const uint128_t& other) const;
    uint128_t operator/(const uint128_t& other) const;
    uint128_t operator%(const uint128_t& other) const;

    // Bitwise
    uint128_t operator|(const uint128_t& other) const;
    uint128_t operator&(const uint128_t& other) const;
    uint128_t operator^(const uint128_t& other) const;
    uint128_t operator~() const;
    uint128_t& operator|=(const uint128_t& other);
    uint128_t& operator&=(const uint128_t& other);
    uint128_t& operator^=(const uint128_t& other);

    // Shifts
    uint128_t operator<<(unsigned int shift) const;
    uint128_t operator>>(unsigned int shift) const;

    // Comparisons
    bool operator<(const uint128_t& other) const;
    bool operator>=(const uint128_t& other) const;
    bool operator==(const uint128_t& other) const;
    bool operator!=(const uint128_t& other) const;
};
const uint128_t UINT128_MAX = uint128_t{0xFFFFFFFFFFFFFFFFULL, 0xFFFFFFFFFFFFFFFFULL};

std::vector<uint8_t> fileToBitstring(const std::filesystem::path& filepath);

// Function prototypes
std::vector<uint8_t> expandSymbolsToBits(const std::vector<uint8_t>& symbols, int word_size);

bool relEpsilonEqual(double A, double B, double maxAbsFactor, double maxRelFactor, uint32_t maxULP);
void free_data(data_t *dp);

void xoshiro_jump(unsigned int jump_count, uint64_t *xoshiro256starstarState);
void seed(uint64_t *xoshiro256starstarState);
uint64_t randomRange64(uint64_t s, uint64_t *xoshiro256starstarState);
double randomUnit(uint64_t *xoshiro256starstarState);
void FYshuffle(uint8_t data[], uint8_t rawdata[], int sample_size, uint64_t *xoshiro256starstarState);

long int sum(const uint8_t arr[], int sample_size);
template<size_t LENGTH>
int sum(const std::array<int, LENGTH>& arr);
template<typename T>
T sum(const std::vector<T>& v) {
    T total = 0;
    for (auto x : v) total += x;
    return total;
}

void map_init(std::map<uint8_t, int>& m);
void map_init(std::map<uint8_t, double>& m);
void map_init(std::map<std::pair<uint8_t, uint8_t>, int>& m);

void calc_proportions(const uint8_t data[], std::vector<double>& p, int sample_size);
void calc_counts(const uint8_t data[], std::vector<int>& c, int sample_size);
double std_dev(const std::vector<int> x, double x_mean);
uint64_t n_choose_2(uint64_t n);

std::vector<uint8_t> substr(const uint8_t text[], int pos, int len, int sample_size);
std::array<uint8_t, 16> fast_substr(const uint8_t text[], int pos, int len);

template<typename T>
T max_vector(const std::vector<T>& vals);
template<typename T>
T max_arr(const T* vals, unsigned int k);

double divide(int a, int b);
double prediction_estimate_function(long double p, long r, long N);
double calc_p_local(long max_run_len, long N, double ldomain);

struct PredictionEstimateResult {
    uint64_t C = 0;       // Correct prediction count
    uint64_t r = 0;       // Max run length of correct predictions
    uint64_t N = 0;       // Total number of predictions
    double p_global = 0.0;
    double p_global_prime = 0.0;
    double p_local = 0.0;
    double curMax = 0.0;  // the "max probability" used for entropy estimate
};

PredictionEstimateResult predictionEstimate(long C, long N, long max_run_len, long k);

uint32_t compressedBitSymbols(const uint8_t *S, long length);
void printVersion(std::string name);
std::string recreateCommandLine(int argc, char* argv[]);

class PostfixDictionary {
    std::map<uint8_t, long> postfixes;
    long curBest;
    uint8_t curPrediction;
    
public:
    PostfixDictionary() { curBest = 0; curPrediction = 0;}

    uint8_t predict(long &count);
    bool incrementPostfix(uint8_t in, bool makeNew);
};
