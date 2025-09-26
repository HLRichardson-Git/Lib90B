
#pragma once

#include <iostream>		// std::cout
#include <string>		// std::string
#include <map>			// std::map
#include <set>			// std::set
#include <string.h>		// strlen
#include <iomanip>		// setw / setfill
#include <stdio.h>
#include <cstdlib>
#include <vector>		// std::vector
#include <time.h>		// time
#include <algorithm>	// std::sort
#include <cmath>		// pow, log2
#include <array>		// std::array
#include <omp.h>		// openmp 4.0 with gcc 4.9
#include <bitset>
#include <mutex>		// std::mutex
#include <assert.h>
#include <cfloat>
#include <math.h>
#include <cstdint>
#include <stdexcept>
#include <cfloat>
#include <limits>

#include "utils.h"

#define SWAP(x, y) do { int s = x; x = y; y = s; } while(0)
#define INOPENINTERVAL(x, a, b) (((a)>(b))?(((x)>(b))&&((x)<(a))):(((x)>(a))&&((x)<(b))))
#define INCLOSEDINTERVAL(x, a, b) (((a)>(b))?(((x)>=(b))&&((x)<=(a))):(((x)>=(a))&&((x)<=(b))))

#define MIN_SIZE 1000000
#define PERMS 10000

//This is the smallest practical value (one can't do better with the double type)
#define RELEPSILON DBL_EPSILON
//This is clearly overkill, but it's difficult to do better without a view into the monotonic function
#define ABSEPSILON DBL_MIN
//#define DBL_INFINITY __builtin_inf ()
#define DBL_INFINITY std::numeric_limits<double>::infinity()
#define ITERMAX 1076
#define ZALPHA 2.5758293035489008
#define ZALPHA_L 2.575829303548900384158L

//Make uint128_t a supported type (standard as of C23)
/*#ifdef __SIZEOF_INT128__
typedef unsigned __int128 uint128_t;
typedef unsigned __int128 uint_least128_t;
# define UINT128_MAX         ((uint128_t)-1)
# define UINT128_WIDTH       128
# define UINT_LEAST128_WIDTH 128
# define UINT_LEAST128_MAX   UINT128_MAX
# define UINT128_C(N)        ((uint_least128_t)+N ## WBU)
#endif*/

typedef struct data_t data_t;

using namespace std;

// Conversions
uint128_t::operator uint64_t() const {
    return low;
}

uint128_t::operator long double() const {
    const long double factor = 18446744073709551616.0L; // 2^64
    return (long double)high * factor + (long double)low;
}

uint128_t::operator double() const {
    const double factor = 18446744073709551616.0; // 2^64
    return (double)high * factor + (double)low;
}

// Addition
uint128_t uint128_t::operator+(const uint128_t& other) const {
    uint128_t result;
    result.low = low + other.low;
    result.high = high + other.high + (result.low < low ? 1 : 0); // carry
    return result;
}
uint128_t& uint128_t::operator+=(const uint128_t& other) {
    *this = *this + other;
    return *this;
}

// Subtraction
uint128_t uint128_t::operator-(const uint128_t& other) const {
    uint128_t result;
    result.low = low - other.low;
    result.high = high - other.high - (low < other.low ? 1 : 0); // borrow
    return result;
}
uint128_t& uint128_t::operator-=(const uint128_t& other) {
    *this = *this - other;
    return *this;
}

// Bitwise
uint128_t uint128_t::operator|(const uint128_t& other) const {
    return uint128_t{low | other.low, high | other.high};
}
uint128_t uint128_t::operator&(const uint128_t& other) const {
    return uint128_t{low & other.low, high & other.high};
}
uint128_t uint128_t::operator^(const uint128_t& other) const {
    return uint128_t{low ^ other.low, high ^ other.high};
}
uint128_t uint128_t::operator~() const {
    return uint128_t{~low, ~high};
}
uint128_t& uint128_t::operator|=(const uint128_t& other) {
    low |= other.low;
    high |= other.high;
    return *this;
}
uint128_t& uint128_t::operator&=(const uint128_t& other) {
    low &= other.low;
    high &= other.high;
    return *this;
}
uint128_t& uint128_t::operator^=(const uint128_t& other) {
    low ^= other.low;
    high ^= other.high;
    return *this;
}

// Shifts
uint128_t uint128_t::operator<<(unsigned int shift) const {
    if (shift >= 128) return uint128_t(0);
    if (shift == 0) return *this;
    uint128_t result;
    if (shift < 64) {
        result.high = (high << shift) | (low >> (64 - shift));
        result.low  = low << shift;
    } else {
        result.high = low << (shift - 64);
        result.low  = 0;
    }
    return result;
}
uint128_t uint128_t::operator>>(unsigned int shift) const {
    if (shift >= 128) return uint128_t(0);
    if (shift == 0) return *this;
    uint128_t result;
    if (shift < 64) {
        result.low  = (low >> shift) | (high << (64 - shift));
        result.high = high >> shift;
    } else {
        result.low  = high >> (shift - 64);
        result.high = 0;
    }
    return result;
}

// Multiplication
uint128_t uint128_t::operator*(const uint128_t& other) const {
    uint64_t a0 = low & 0xFFFFFFFFULL;
    uint64_t a1 = low >> 32;
    uint64_t a2 = high & 0xFFFFFFFFULL;
    uint64_t a3 = high >> 32;

    uint64_t b0 = other.low & 0xFFFFFFFFULL;
    uint64_t b1 = other.low >> 32;
    uint64_t b2 = other.high & 0xFFFFFFFFULL;
    uint64_t b3 = other.high >> 32;

    uint128_t result;
    uint64_t r0 = a0 * b0;
    uint64_t r1 = a0 * b1 + a1 * b0;
    uint64_t r2 = a0 * b2 + a1 * b1 + a2 * b0;
    uint64_t r3 = a0 * b3 + a1 * b2 + a2 * b1 + a3 * b0;
    uint64_t r4 = a1 * b3 + a2 * b2 + a3 * b1;
    uint64_t r5 = a2 * b3 + a3 * b2;
    uint64_t r6 = a3 * b3;

    uint64_t carry = 0;
    result.low = r0 + (r1 << 32);
    if (result.low < r0) carry++;
    result.high = r3 + r4 + r5 + r6 + carry;

    return result;
}

// Division
uint128_t uint128_t::operator/(const uint128_t& other) const {
    if (other == uint128_t(0)) throw std::runtime_error("Division by zero");
    uint128_t quotient = 0;
    uint128_t remainder = 0;

    for (int i = 127; i >= 0; --i) {
        remainder = (remainder << 1);
        remainder.low |= ((*this >> i).low & 1);
        if (remainder >= other) {
            remainder -= other;
            quotient = quotient | (uint128_t(1) << i);
        }
    }
    return quotient;
}

// Modulo
uint128_t uint128_t::operator%(const uint128_t& other) const {
    if (other == uint128_t(0)) throw std::runtime_error("Modulo by zero");
    uint128_t quotient = 0;
    uint128_t remainder = 0;

    for (int i = 127; i >= 0; --i) {
        remainder = (remainder << 1);
        remainder.low |= ((*this >> i).low & 1);
        if (remainder >= other) {
            remainder -= other;
            quotient = quotient | (uint128_t(1) << i);
        }
    }
    return remainder;
}

// Comparisons
bool uint128_t::operator<(const uint128_t& other) const {
    return high < other.high || (high == other.high && low < other.low);
}
bool uint128_t::operator>=(const uint128_t& other) const {
    return !(*this < other);
}
bool uint128_t::operator==(const uint128_t& other) const {
    return high == other.high && low == other.low;
}
bool uint128_t::operator!=(const uint128_t& other) const {
    return !(*this == other);
}



std::vector<uint8_t> expandSymbolsToBits(const std::vector<uint8_t>& symbols, int word_size) {
    std::vector<uint8_t> bits;
    bits.reserve(symbols.size() * word_size);

    for (uint8_t sym : symbols) {
        // Extract bits from MSB to LSB
        for (int i = word_size - 1; i >= 0; --i) {
            uint8_t bit = (sym >> i) & 0x01;
            bits.push_back(bit);
        }
    }
    return bits;
}

//This generally performs a check for relative closeness, but (if that check would be nonsense)
//it can check for an absolute separation, using either the distance between the numbers, or
//the number of ULPs that separate the two numbers.
//See the following for details and discussion of this approach:
//https://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/
//https://floating-point-gui.de/errors/comparison/
//https://www.boost.org/doc/libs/1_62_0/libs/test/doc/html/boost_test/testing_tools/extended_comparison/floating_point/floating_points_comparison_theory.html
//Knuth AoCP vol II (section 4.2.2)
//Tested using modified test cases from https://floating-point-gui.de/errors/NearlyEqualsTest.java
bool relEpsilonEqual(double A, double B, double maxAbsFactor, double maxRelFactor, uint32_t maxULP)
{
   double diff;
   double absA, absB;
   uint64_t Aint;
   uint64_t Bint;

   assert(sizeof(uint64_t) == sizeof(double));
   assert(maxAbsFactor >= 0.0);
   assert(maxRelFactor >= 0.0);

   ///NaN is by definition not equal to anything (including itself)
   if(std::isnan(A) || std::isnan(B)) {
      return false;
   }

   //Deals with equal infinities, and the corner case where they are actually copies
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"
   if(A==B) {
      return true;
   }
#pragma GCC diagnostic pop

   //If either is infinity, but they are not equal, then they aren't close.
   if(std::isinf(A) || std::isinf(B)) {
      return false;
   }

   absA = fabs(A);
   absB = fabs(B);
   //Make sure that A is the closest to 0.
   if(absA > absB) {
      double tmp;

      //Swap A and B
      tmp = B;
      B=A;
      A=tmp;

      //Swap absA and absB
      tmp = absB;
      absB = absA;
      absA = tmp;
   }

   //Capture the difference of the largest magnitude from the smallest magnitude
   diff=fabs(B-A);

   //Is absA, diff, or absB * maxRelFactor subnormal?
   //Did diff overflow?
   //if absA is subnormal (effectively 0) or 0, then relative difference isn't meaningful, as fabs(B-A)/B≈1 for all values of B
   //In the instance of overflows, the resulting relative comparison will be nonsense.
   if((absA < DBL_MIN) || (diff < DBL_MIN) || std::isinf(diff) || (absB * maxRelFactor < DBL_MIN)) {
      //Yes. Relative closeness is going to be nonsense
      return diff <= maxAbsFactor;
   } else {
      //No. Using relative closeness is probably the right thing to do.
      //Proceeding roughly as per Knuth AoCP vol II (section 4.2.2)
      if(diff <= absB * maxRelFactor) {
         //These are relatively close
         return true;
      } 
   }

   //Neither A or B is subnormal, and they aren't close in the conventional sense, 
   //but perhaps that's just due to IEEE representation. Check to see if the value is within maxULP ULPs.

   //We can't meaningfully compare non-zero values with 0.0 in this way,
   //but absA >= DBL_MIN if we're here, so neither value is 0.0.

   //if they aren't the same sign, then these can't be only a few ULPs away from each other
   if(signbit(A) != signbit(B)) {
      return false;
   }

   //Note, casting from one type to another is undefined behavior, but memcpy will necessarily work
   memcpy(&Aint, &absA, sizeof(double));
   memcpy(&Bint, &absB, sizeof(double));
   //This should be true by the construction of IEEE doubles
   assert(Bint > Aint);

   return (Bint - Aint <= maxULP);
}


void free_data(data_t *dp){
	if(dp->symbols != NULL) free(dp->symbols);
	if(dp->rawsymbols != NULL) free(dp->rawsymbols);
	if((dp->word_size > 1) && (dp->bsymbols != NULL)) free(dp->bsymbols);
} 

/* This is xoshiro256** 1.0*/
/*This implementation is derived from David Blackman and Sebastiano Vigna, which they placed into
the public domain. See http://xoshiro.di.unimi.it/xoshiro256starstar.c
*/
static inline uint64_t rotl(const uint64_t x, int k) {
	return (x << k) | (x >> (64 - k));
}

static inline uint64_t xoshiro256starstar(uint64_t *xoshiro256starstarState)
{
	const uint64_t result_starstar = rotl(xoshiro256starstarState[1] * 5, 7) * 9;
	const uint64_t t = xoshiro256starstarState[1] << 17;

	xoshiro256starstarState[2] ^= xoshiro256starstarState[0];
	xoshiro256starstarState[3] ^= xoshiro256starstarState[1];
	xoshiro256starstarState[1] ^= xoshiro256starstarState[2];
	xoshiro256starstarState[0] ^= xoshiro256starstarState[3];

	xoshiro256starstarState[2] ^= t;

	xoshiro256starstarState[3] = rotl(xoshiro256starstarState[3], 45);

   return result_starstar;
}

/* This is the jump function for the generator. It is equivalent
   to 2^128 calls to xoshiro256starstar(); it can be used to generate 2^128
   non-overlapping subsequences for parallel computations. */
void xoshiro_jump(unsigned int jump_count, uint64_t *xoshiro256starstarState) {
	static const uint64_t JUMP[] = { 0x180ec6d33cfd0aba, 0xd5a61266f0c9392c, 0xa9582618e03fc9aa, 0x39abdc4529b1661c };
	uint64_t s0 = 0;
	uint64_t s1 = 0;
	uint64_t s2 = 0;
	uint64_t s3 = 0;

	for(unsigned int j=0; j < jump_count; j++) {
		for(unsigned int i = 0; i < sizeof JUMP / sizeof *JUMP; i++)
			for(unsigned int b = 0; b < 64; b++) {
				if (JUMP[i] & ((uint64_t)1) << b) {
					s0 ^= xoshiro256starstarState[0];
					s1 ^= xoshiro256starstarState[1];
					s2 ^= xoshiro256starstarState[2];
					s3 ^= xoshiro256starstarState[3];
				}
				xoshiro256starstar(xoshiro256starstarState);	
			}
			
		xoshiro256starstarState[0] = s0;
		xoshiro256starstarState[1] = s1;
		xoshiro256starstarState[2] = s2;
		xoshiro256starstarState[3] = s3;
	}
}

//This seeds using an external source
//We use /dev/urandom here. 
//We could alternately use the RdRand (or some other OS or HW source of pseudo-random numbers)
void seed(uint64_t *xoshiro256starstarState){
	FILE *infp;

	if((infp=fopen("/dev/urandom", "rb"))==NULL) {
		perror("Can't open random source. Reverting to a deterministic seed.");
		exit(-1);
	} 

	if(fread(xoshiro256starstarState, sizeof(uint64_t), 4, infp)!=4) {
		perror("Can't read random seed");
		exit(-1);
	}

	if(fclose(infp)!=0) {
		perror("Couldn't close random source");
		exit(-1);
	}
}

/*Return an integer in the range [0, high], without modular bias*/
/*This is a slight modification of Lemire's approach (as we want [0,s] rather than [0,s)*/
/*See "Fast Random Integer Generation in an Interval" by Lemire (2018) (https://arxiv.org/abs/1805.10941) */
 /* The relevant text explaining the central factor underlying this opaque approach is:
  * "Given an integer x ∈ [0, 2^L), we have that (x × s) ÷ 2^L ∈ [0, s). By multiplying by s, we take
  * integer values in the range [0, 2^L) and map them to multiples of s in [0, s × 2^L). By dividing by 2^L,
  * we map all multiples of s in [0, 2^L) to 0, all multiples of s in [2^L, 2 × 2^L) to one, and so forth. The
  * (i + 1)th interval is [i × 2^L, (i + 1) × 2^L). By Lemma 2.1, there are exactly floor(2^L/s) multiples of s in
  * intervals [i × 2^L + (2^L mod s), (i + 1) × 2^L) since s divides the size of the interval (2^L − (2^L mod s)).
  * Thus if we reject the multiples of s that appear in [i × 2^L, i × 2^L + (2^L mod s)), we get that all
  * intervals have exactly floor(2^L/s) multiples of s."
  *
  * This approach allows us to avoid _any_ modular reductions with high probability, and at worst case one
  * reduction. It's an opaque approach, but lovely.
  */
uint64_t randomRange64(uint64_t s, uint64_t *xoshiro256starstarState){
	uint64_t x;
	uint128_t m;
	uint64_t l;

	x = xoshiro256starstar(xoshiro256starstarState);

	if(UINT64_MAX == s) {
		return x;
	} else {
		s++; // We want an integer in the range [0,s], not [0,s)
		m = (uint128_t)x * (uint128_t)s;
		l = (uint64_t)m; //This is m mod 2^64

		if(l<s) {
			uint64_t t = ((uint64_t)(-s)) % s; //t = (2^64 - s) mod s (by definition of unsigned arithmetic in C)
			while(l < t) {
				x = xoshiro256starstar(xoshiro256starstarState);
				m = (uint128_t)x * (uint128_t)s;
				l = (uint64_t)m; //This is m mod 2^64
			}
		}

		return (uint64_t)(m >> 64U); //return floor(m/2^64)
	}
}

/*
 * This function produces a double that is uniformly distributed in the interval [0, 1).
 * Note that 2^53 is the largest integer that can be represented in a 64 bit IEEE 754 double, such that all 
 * smaller positive integers can also be represented. Shifting the initial random 64-bit value right by 11 
 * bits makes the result only in the lower 53 bits, so the resulting integer is in the range [0, 2^53 - 1].
 * 1.1102230246251565e-16 (0x1.0p-53) is 2^(-53). Multiplying by this value just effects the exponent of the 
 * resulting double, not the significand. We get a double uniformly distributed in the range [0, 1).  
 * The delta between adjacent values is 2^(-53).
 */
double randomUnit(uint64_t *xoshiro256starstarState) {
	return((xoshiro256starstar(xoshiro256starstarState) >> 11) * 1.1102230246251565e-16);
}

// Fisher-Yates Fast (in place) shuffle algorithm
void FYshuffle(uint8_t data[], uint8_t rawdata[], const int sample_size, uint64_t *xoshiro256starstarState) {
	long int r;
	static mutex shuffle_mutex;
	unique_lock<mutex> lock(shuffle_mutex);

	for (long int i = sample_size - 1; i > 0; --i) {
		r = (long int)randomRange64((uint64_t)i, xoshiro256starstarState);
		SWAP(data[r], data[i]);
		SWAP(rawdata[r], rawdata[i]);
	}
}

// Quick sum array  // TODO
long int sum(const uint8_t arr[], const int sample_size) {
	long int sum = 0;
	for (long int i = 0; i < sample_size; ++i) {
		sum += arr[i];
	}

	return sum;
}

// Quick sum std::array // TODO
template<size_t LENGTH>
int sum(const array<int, LENGTH> &arr) {
	int sum = 0;
	for (int i = 0; i < LENGTH; ++i) {
		sum += arr[i];
	}

	return sum;
}

// Quick sum vector
template<typename T>
T sum(const vector<T> &v) {
	T sum = 0;
	for (unsigned long int i = 0; i < v.size(); ++i) {
		sum += v[i];
	}

	return sum;
}

// Calculate baseline statistics
// Finds mean, median, and whether or not the data is binary
void calc_stats(const data_t *dp, double &rawmean, double &median) {

	// Calculate mean
	rawmean = sum(dp->rawsymbols, dp->len) / (double)dp->len;

	// Sort in a vector for median/min/max
	vector<uint8_t> v(dp->symbols, dp->symbols + dp->len);
	sort(v.begin(), v.end());

	long int half = dp->len / 2;
	if(dp->alph_size == 2) {
		//This isn't necessarily true, but we are supposed to set it this way.
		//See 5.1.5, 5.1.6.
		median = 0.5;
	} else {
		if((dp->len & 1) == 1) {
			//the length is odd
			median = v[half];
		} else {
			//the length is even
			median = (v[half] + v[half - 1]) / 2.0;
		}
	}
}


// Map initialization for integers
void map_init(map<uint8_t, int> &m) {
	for (int i = 0; i < 256; i++) {
		m[i] = 0;
	}
}

// Map initialization for doubles
void map_init(map<uint8_t, double> &m) {
	for (int i = 0; i < 256; i++) {
		m[i] = 0.0;
	}
}

// Map initialization for pair<uint8_t, uint8_t> to int
void map_init(map<pair<uint8_t, uint8_t>, int> &m) {
	for (int i = 0; i < 256; i++) {
		for (int j = 0; j < 256; j++) {
			m[pair<uint8_t, uint8_t>(i, j)] = 0;
		}
	}
}

// Calculates proportions of each value as an index
void calc_proportions(const uint8_t data[], vector<double> &p, const int sample_size) {
	for (int i = 0; i < sample_size; i++) {
		p[data[i]] += (1.0 / sample_size);
	}
}

// Calculates proportions of each value as an index
void calc_counts(const uint8_t data[], vector<int> &c, const int sample_size) {
	for (int i = 0; i < sample_size; i++) {
		c[data[i]] ++;
	}
}

// Determines the standard deviation of a dataset
double std_dev(const vector<int> x, const double x_mean) {
	double sum = 0.0;

	for (unsigned int i = 0; i < x.size(); i++) {
		sum += pow(x[i] - x_mean, 2);
	}

	return sqrt(sum / x.size());
}

// Quick formula for n choose 2 (which can be simplified to [n^2 - n] / 2)
long int n_choose_2(const long int n) {
	return ((n*n) - n) / 2;
}

vector<uint8_t> substr(const uint8_t text[], const int pos, const int len, const int sample_size) {
	int substr_len = len;

	if (pos + len > sample_size) {
		substr_len = sample_size - pos;
	}

	vector<uint8_t> substring;

	for (int i = 0; i < substr_len; i++) {
		substring.push_back(text[pos + i]);
	}

	return substring;
}

// Fast substring with no bounds checking
array<uint8_t, 16> fast_substr(const uint8_t text[], const int pos, const int len) {
	array<uint8_t, 16> substring = { 0 };

	for (int i = 0; i < len; i++) {
		substring[i] = text[pos + i];
	}

	return substring;
}

template<typename T>
T max_vector(const vector<T> &vals) {
	T max = vals[0];
	for (unsigned int i = 0; i < vals.size(); i++) {
		if (vals[i] > max) {
			max = vals[i];
		}
	}

	return max;
}

template<typename T>
T max_arr(const T* vals, const unsigned int k){
	T max = vals[0];
	for (unsigned int i = 0; i < k; i++){
		if (vals[i] > max) {
			max = vals[i];
		}
	}

	return max;
}

double divide(const int a, const int b) {
	return ((double)a / (double)b);
}

double prediction_estimate_function(long double p, long r, long N){
	long double q, x, xlast=0.0L;
	assert(p > 0.0L); //In fact, it is >= 1/k
	assert(p < 1.0L);

	q = 1.0L-p;
	x = 1.0L;

	//We know that 
	// * x is in the interval [1,1/p] which is a subset of [1,k]
	// * the sequence is monotonic up (so x >= xlast).
	//As such we don't need much fancyness for looking for "equality"
	for(int i = 0; (i <= 65) && ((x - xlast) > (LDBL_EPSILON*x)); i++) {
		xlast = x;
		x = 1.0L + q*powl(p, r)*powl(x, r+1.0L);
		//We expect this convergence to be monotonic up.
		assert(x >= xlast);
		//We expect x<=1/p
		assert(p*x <= 1.0L);
	}

	return((double)(logl(1.0L-p*x) - logl((r+1.0L-r*x)*q) - (N+1.0L)*logl(x)));
}

double calc_p_local(long max_run_len, long N, double ldomain){
	int j;
	double p, log_alpha;
	double lastP, pVal;
	double lvalue, hvalue;
	double hbound, lbound;
	double hdomain;

	// binary search for p_local
	log_alpha = log(0.99);
	
	hdomain = 1.0;

	lbound = ldomain;
	hbound = hdomain;

	lvalue = DBL_INFINITY;
	hvalue = -DBL_INFINITY;

	//Note that the bounds are in [0,1], so overflows aren't an issue
	//But underflows are.
	p = (lbound + hbound) / 2.0;
	pVal = prediction_estimate_function(p, max_run_len+1, N);

	//We don't need the initial pVal invariant, as our initial bounds are infinite.
	//We don't need the initial bounds, as they are set to the domain bounds
	for(j=0; j<ITERMAX; j++) {
		//Have we reached "equality"?
		if(relEpsilonEqual(pVal, log_alpha, ABSEPSILON, RELEPSILON, 4)) break;

		//Now update based on the found pVal
		if(log_alpha < pVal) {
			lbound = p;
			lvalue = pVal;
		} else {
			hbound = p;
			hvalue = pVal;
		}

		//We now verify that ldomain <= lbound < p < hbound <= hdomain
		//and that target in [ lvalue, hvalue ]
		if(lbound >= hbound) {
			p = fmin(fmax(lbound, hbound),hdomain);
			break;
		}

		//invariant. If this isn't true, then we can't evaluate here.
		if(!(INCLOSEDINTERVAL(lbound, ldomain, hdomain) && INCLOSEDINTERVAL(hbound,  ldomain, hdomain))) {
			p = hdomain;
			break;
		}

		//invariant. If this isn't true, then seeking the value within this interval doesn't make sense.
		if(!INCLOSEDINTERVAL(log_alpha, lvalue, hvalue)) {
			p = hdomain;
			break;
		}

		//Update p
		lastP = p;
		p = (lbound + hbound) / 2.0;

		//invariant. If this isn't true, then further calculation isn't really meaningful.
		if(!INOPENINTERVAL(p,  lbound, hbound)) {
			p = hbound;
			break;
		}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"
		//Look for a cycle
		if(lastP == p) {
			p = hbound;
			break;
		}
#pragma GCC diagnostic pop

		pVal = prediction_estimate_function(p, max_run_len+1, N);

		//invariant: If this isn't true, then this isn't loosely monotonic
		if(!INCLOSEDINTERVAL(pVal, lvalue, hvalue)) {
			p = hbound;
			break;
		}
	}//for loop

	return p;
}

PredictionEstimateResult predictionEstimate(long C, long N, long max_run_len, long k) {
    PredictionEstimateResult res;
    res.C = C;
    res.N = N;
    res.r = max_run_len + 1;

    double curMax = 1.0 / static_cast<double>(k);

    res.p_global = static_cast<double>(C) / static_cast<double>(N);

    if (res.p_global > 0) {
        res.p_global_prime = std::min(
            1.0,
            res.p_global + ZALPHA * sqrt((res.p_global * (1.0 - res.p_global)) / (static_cast<double>(N) - 1.0))
        );
    } else {
        res.p_global_prime = 1.0 - pow(0.01, 1.0 / static_cast<double>(N));
    }

    curMax = std::max(curMax, res.p_global_prime);

    double p_local = -1.0;
    if ((curMax < 1.0) && (prediction_estimate_function(curMax, max_run_len + 1, N) > log(0.99))) {
        p_local = calc_p_local(max_run_len, N, curMax);
        curMax = std::max(curMax, p_local);
    }
    res.p_local = (p_local > 0.0) ? p_local : 0.0;

    res.curMax = curMax;
    return res;
}

//The idea here is that we've given an array of pointers (binaryDict). 
//We are trying to produce the address of the length-2 array associated with the length-d prefix "b".
// array The dth index is d-1, so we first find the start of the address space (binaryDict[(d)-1])
//We take the least significant d bits from "b": this is the expression "(b) & ((1U << (d)) - 1)"
//We then multiply this by 2 (as each pattern is associated with a length-2 array) by left shifting by 1.
#define BINARYDICTLOC(d, b) (binaryDict[(d)-1] + (((b) & ((1U << (d)) - 1))<<1))

uint32_t compressedBitSymbols(const uint8_t *S, long length) {
   uint32_t retPattern;
   long j;

   assert(length<=32);

   retPattern = 0;

   for(j=0; j<length; j++) {
      assert(S[j] <= 1);
      retPattern = (retPattern << 1) | S[j];
   }

   return retPattern;
}

static void printVersion(string name) {
    cout << name << " " << VERSION << "\n\n";
    cout << "Disclaimer: ";
    cout << "NIST-developed software is provided by NIST as a public service. You may use, copy, and distribute copies of the software in any medium, provided that you keep intact this entire notice. You may improve, modify, and create derivative works of the software or any portion of the software, and you may copy and distribute such modifications or works. Modified works should carry a notice stating that you changed the software and should note the date and nature of any such change. Please explicitly acknowledge the National Institute of Standards and Technology as the source of the software.";
    cout << "\n\n";
    cout << "NIST-developed software is expressly provided \"AS IS.\" NIST MAKES NO WARRANTY OF ANY KIND, EXPRESS, IMPLIED, IN FACT, OR ARISING BY OPERATION OF LAW, INCLUDING, WITHOUT LIMITATION, THE IMPLIED WARRANTY OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NON-INFRINGEMENT, AND DATA ACCURACY. NIST NEITHER REPRESENTS NOR WARRANTS THAT THE OPERATION OF THE SOFTWARE WILL BE UNINTERRUPTED OR ERROR-FREE, OR THAT ANY DEFECTS WILL BE CORRECTED. NIST DOES NOT WARRANT OR MAKE ANY REPRESENTATIONS REGARDING THE USE OF THE SOFTWARE OR THE RESULTS THEREOF, INCLUDING BUT NOT LIMITED TO THE CORRECTNESS, ACCURACY, RELIABILITY, OR USEFULNESS OF THE SOFTWARE.";
    cout << "\n\n";
    cout << "You are solely responsible for determining the appropriateness of using and distributing the software and you assume all risks associated with its use, including but not limited to the risks and costs of program errors, compliance with applicable laws, damage to or loss of data, programs or equipment, and the unavailability or interruption of operation. This software is not intended to be used in any situation where a failure could cause risk of injury or damage to property. The software developed by NIST employees is not subject to copyright protection within the United States.";
    cout << "\n\n";
}

static string recreateCommandLine(int argc, char* argv[]) {
    string commandLine = "";
    for(int i = 0; i < argc; ++i) {
        commandLine.append(argv[i]);
        if((i + 1) != argc) {
            commandLine.append(" ");
        }
    }
    return commandLine;
}

uint8_t PostfixDictionary::predict(long &count) {assert(curBest > 0); count = curBest; return curPrediction;}

bool PostfixDictionary::incrementPostfix(uint8_t in, bool makeNew) {
	map<uint8_t, long>::iterator curp = postfixes.find(in);
	long curCount;
	bool newEntry=false;

	if(curp != postfixes.end()) {
		//The entry is already there. We always increment in this case.
		curCount = ++(curp->second);
	} else if(makeNew) {
		//The entry is not here, but we are allowed to create a new entry
		newEntry = true;
		curCount = postfixes[in] = 1;
	} else {
		//The entry is not here, we are not allowed to create a new entry
		return false;
	}

	//Only instances where curCount is set and an increment was performed get here
	if((curCount > curBest) || ((curCount == curBest) && (in > curPrediction))) { 
		curPrediction = in; 
		curBest = curCount; 
	} 

	return newEntry;
}
