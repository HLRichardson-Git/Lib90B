#include <vector>
#include <algorithm>
#include <cstdint>
#include <stdexcept>
#include <fstream>
#include <iterator>
#include <unordered_set>
#include <unordered_map>
#include <cmath>
#include <cassert>

#include <lib90b/chi_square.h>
#include <lib90b/entropy_tests.h>
#include "utils.h"

/* -----------------------------
   Cephes Gamma & Chi-Square Helpers
----------------------------- */

static double MACHEP = 1.11022302462515654042E-16;     // 2**-53
static double MAXLOG = 7.09782712893383996732224E2;    // log(MAXNUM)
static double MAXNUM = 1.7976931348623158E308;         // 2**1024*(1-MACHEP)
static double PI     = 3.14159265358979323846;         // pi

static double big = 4.503599627370496e15;
static double biginv = 2.22044604925031308085e-16;

static int sgngam = 0;

static double A[] = {
   8.11614167470508450300E-4,
   -5.95061904284301438324E-4,
   7.93650340457716943945E-4,
   -2.77777777730099687205E-3,
   8.33333333333331927722E-2
};
static double B[] = {
   -1.37825152569120859100E3,
   -3.88016315134637840924E4,
   -3.31612992738871184744E5,
   -1.16237097492762307383E6,
   -1.72173700820839662146E6,
   -8.53555664245765465627E5
};
static double C[] = {
   -3.51815701436523470549E2,
   -1.70642106651881159223E4,
   -2.20528590553854454839E5,
   -1.13933444367982507207E6,
   -2.53252307177582951285E6,
   -2.01889141433532773231E6
};

#define MAXLGM 2.556348e305

double cephes_igamc(double a, double x);

static double cephes_polevl(double x, double *coef, int N)
{
   double ans;
   int i;
   double *p;

   p = coef;
   ans = *p++;
   i = N;

   do {
      ans = ans * x + *p++;
   } while (--i);

   return ans;
}

static double cephes_p1evl(double x, double *coef, int N)
{
   double ans;
   double *p;
   int i;

   p = coef;
   ans = x + *p++;
   i = N - 1;

   do {
      ans = ans * x + *p++;
   } while (--i);

   return ans;
}

static double cephes_lgam(double x)
{
   double p, q, u, w, z;
   int i;

   sgngam = 1;

   if (x < -34.0) {
      q = -x;
      w = cephes_lgam(q);
      p = floor(q);

      if (relEpsilonEqual(p, q, DBL_EPSILON, DBL_EPSILON, 4)) {
         goto loverf;
      }

      i = (int)p;

      if ((i & 1) == 0) {
         sgngam = -1;
      } else {
         sgngam = 1;
      }

      z = q - p;

      if (z > 0.5) {
         p += 1.0;
         z = p - q;
      }

      z = q * sin(PI * z);

      if (relEpsilonEqual(z, 0.0, DBL_EPSILON, DBL_EPSILON, 4)) {
         goto loverf;
      }

      z = log(PI) - log(z) - w;
      return z;
   }

   if (x < 13.0) {
      z = 1.0;
      p = 0.0;
      u = x;

      while (u >= 3.0) {
         p -= 1.0;
         u = x + p;
         z *= u;
      }

      while (u < 2.0) {
         if (relEpsilonEqual(u, 0.0, DBL_EPSILON, DBL_EPSILON, 4)) {
            goto loverf;
         }

         z /= u;
         p += 1.0;
         u = x + p;
      }

      if (z < 0.0) {
         sgngam = -1;
         z = -z;
      } else {
         sgngam = 1;
      }

      if (relEpsilonEqual(u, 2.0, DBL_EPSILON, DBL_EPSILON, 4)) {
         return (log(z));
      }

      p -= 2.0;
      x = x + p;
      p = x * cephes_polevl(x, B, 5) / cephes_p1evl(x, (double *)C, 6);

      return log(z) + p;
   }

   if (x > MAXLGM) {
loverf:
      fprintf(stderr, "lgam: OVERFLOW\n");
      return sgngam * MAXNUM;
   }

   q = (x - 0.5) * log(x) - x + log(sqrt(2 * PI));

   if (x > 1.0e8) {
      return q;
   }

   p = 1.0 / (x * x);

   if (x >= 1000.0)
      q += ((7.9365079365079365079365e-4 * p - 2.7777777777777777777778e-3) * p + 0.0833333333333333333333) / x;
   else {
      q += cephes_polevl(p, A, 4) / x;
   }

   return q;
}

static double cephes_igam(double a, double x)
{
   double ans, ax, c, r;

   if ((x <= 0) || (a <= 0)) {
      return 0.0;
   }

   if ((x > 1.0) && (x > a)) {
      return 1.e0 - cephes_igamc(a, x);
   }

   ax = a * log(x) - x - cephes_lgam(a);

   if (ax < -MAXLOG) {
      fprintf(stderr, "igam: UNDERFLOW\n");
      return 0.0;
   }

   ax = exp(ax);

   r = a;
   c = 1.0;
   ans = 1.0;

   do {
      r += 1.0;
      c *= x / r;
      ans += c;
   } while (c / ans > MACHEP);

   return ans * ax / a;
}

double cephes_igamc(double a, double x)
{
   double ans, ax, c, yc, r, t, y, z;
   double pk, pkm1, pkm2, qk, qkm1, qkm2;

   if ((x <= 0) || (a <= 0)) {
      return (1.0);
   }

   if ((x < 1.0) || (x < a)) {
      return (1.e0 - cephes_igam(a, x));
   }

   ax = a * log(x) - x - cephes_lgam(a);

   if (ax < -MAXLOG) {
      fprintf(stderr, "igamc: UNDERFLOW\n");
      return 0.0;
   }

   ax = exp(ax);

   y = 1.0 - a;
   z = x + y + 1.0;
   c = 0.0;
   pkm2 = 1.0;
   qkm2 = x;
   pkm1 = x + 1.0;
   qkm1 = z * x;
   ans = pkm1 / qkm1;

   do {
      c += 1.0;
      y += 1.0;
      z += 2.0;
      yc = y * c;
      pk = pkm1 * z - pkm2 * yc;
      qk = qkm1 * z - qkm2 * yc;

      if (!relEpsilonEqual(qk, 0.0, DBL_EPSILON, DBL_EPSILON, 4)) {
         r = pk / qk;
         t = fabs((ans - r) / r);
         ans = r;
      } else {
         t = 1.0;
      }

      pkm2 = pkm1;
      pkm1 = pk;
      qkm2 = qkm1;
      qkm1 = qk;

      if (fabs(pk) > big) {
         pkm2 *= biginv;
         pkm1 *= biginv;
         qkm2 *= biginv;
         qkm1 *= biginv;
      }
   } while (t > MACHEP);

   return ans * ax;
}

static double chi_square_pvalue(double x, double k) {
    return cephes_igamc(k / 2.0, x / 2.0);
}

/* -----------------------------
   Tuple Structure for Chi-Square
----------------------------- */
struct tupleTranslateEntry {
    uint16_t tuple;
    double expectation;
    int bin;
};

static bool expectationOrder(const tupleTranslateEntry &a, const tupleTranslateEntry &b) {
    if (a.expectation != b.expectation) return a.expectation < b.expectation;
    else return a.tuple < b.tuple;
}

static bool tupleOrder(const tupleTranslateEntry &a, const tupleTranslateEntry &b) {
    return a.tuple < b.tuple;
}

/* -----------------------------
   Independence Test Helpers
----------------------------- */
static void independence_calc_expectations(const std::vector<double> &p,
                                           std::vector<tupleTranslateEntry> &e,
                                           int sample_size)
{
    size_t alphabet_size = p.size();
    assert(alphabet_size <= UINT8_MAX + 1);
    for (size_t i = 0; i < alphabet_size; i++) {
        for (size_t j = 0; j < alphabet_size; j++) {
            uint16_t index = static_cast<uint16_t>(i * alphabet_size + j);
            e[index].tuple = index;
            e[index].expectation = p[i] * p[j] * floor(sample_size * 0.5);
            e[index].bin = -1;
        }
    }
}

static void allocate_bins(std::vector<tupleTranslateEntry> &e, std::vector<double> &bin_exp)
{
    int current_bin = 0;
    double current_expectation = 0.0;

    for (size_t i = 0; i < e.size(); i++) {
        if (current_expectation >= 5.0) {
            bin_exp.push_back(current_expectation);
            current_bin++;
            current_expectation = 0.0;
        }

        e[i].bin = current_bin;
        current_expectation += e[i].expectation;
    }

    // If the current_bin is 0, we can't combine anything
    if ((current_bin != 0) && (current_expectation < 5.0)) {
        // Combine the last two bins
        for (size_t i = e.size() - 1; e[i].bin == current_bin; i--) {
            e[i].bin = current_bin - 1;
        }
        bin_exp[current_bin - 1] += current_expectation;
    } else {
        bin_exp.push_back(current_expectation);
    }
}

static void independence_calc_observed(const uint8_t data[],
                                      const std::vector<tupleTranslateEntry> &e,
                                      std::vector<int> &o,
                                      int sample_size,
                                      int alphabet_size,
                                      const std::unordered_map<uint8_t, int>& symbol_map)
{
    for (int j = 0; j < sample_size - 1; j += 2) {
        // Map symbols to consecutive indices
        auto it1 = symbol_map.find(data[j]);
        auto it2 = symbol_map.find(data[j + 1]);
        
        if (it1 == symbol_map.end() || it2 == symbol_map.end()) {
            throw std::runtime_error("independence_calc_observed: unmapped symbol found");
        }
        
        int val1 = it1->second;
        int val2 = it2->second;
        uint16_t index = static_cast<uint16_t>((val1 * alphabet_size) + val2);
        
        if (index >= e.size()) {
            throw std::runtime_error("independence_calc_observed: index out of range");
        }
        
        int bin = e[index].bin;
        if (bin < 0 || static_cast<size_t>(bin) >= o.size()) {
            throw std::runtime_error("independence_calc_observed: bin index out of range");
        }
        
        o[bin]++;
    }
}

static double calc_T(const std::vector<double> &bin_expectations, const std::vector<int> &o)
{
    double T = 0.0;
    assert(bin_expectations.size() == o.size());

    for (size_t i = 0; i < bin_expectations.size(); i++) {
        T += pow((o[i] - bin_expectations[i]), 2) / bin_expectations[i];
    }

    return T;
}

/* -----------------------------
   Goodness-of-Fit Helpers
----------------------------- */
static void goodness_of_fit_calc_observed(const uint8_t data[],
                                         const std::vector<tupleTranslateEntry> &e,
                                         std::vector<int> &o,
                                         int sample_size,
                                         const std::unordered_map<uint8_t, int>& symbol_map)
{
    for (int j = 0; j < sample_size; j++) {
        // Map symbol to consecutive index
        auto it = symbol_map.find(data[j]);
        if (it == symbol_map.end()) {
            throw std::runtime_error("goodness_of_fit_calc_observed: unmapped symbol found");
        }
        
        int idx = it->second;
        int bin = e[idx].bin;
        
        if (bin < 0 || static_cast<size_t>(bin) >= o.size()) {
            throw std::runtime_error("goodness_of_fit_calc_observed: bin index out of range");
        }
        
        o[bin]++;
    }
}

/* -----------------------------
   Chi-Square Test Implementations
----------------------------- */
void chi_square_independence(const uint8_t data[],
                             double &score,
                             int &df,
                             int sample_size,
                             int alphabet_size,
                             const std::unordered_map<uint8_t, int>& symbol_map)
{
    if (sample_size < 2 || alphabet_size < 2) {
        score = 0.0;
        df = 0;
        return;
    }

    // Build remapped data
    std::vector<uint8_t> remapped_data(sample_size);
    for (int i = 0; i < sample_size; i++) {
        auto it = symbol_map.find(data[i]);
        if (it == symbol_map.end()) {
            throw std::runtime_error("chi_square_independence: unmapped symbol found");
        }
        remapped_data[i] = static_cast<uint8_t>(it->second);
    }

    std::vector<double> p(alphabet_size, 0.0);
    calc_proportions(remapped_data.data(), p, sample_size);

    size_t tuple_size = static_cast<size_t>(alphabet_size) * alphabet_size;
    std::vector<tupleTranslateEntry> e(tuple_size);
    independence_calc_expectations(p, e, sample_size);
    std::sort(e.begin(), e.end(), expectationOrder);

    std::vector<double> bin_exp;
    allocate_bins(e, bin_exp);
    
    if (bin_exp.empty()) {
        score = 0.0;
        df = 0;
        return;
    }

    std::sort(e.begin(), e.end(), tupleOrder);

    std::vector<int> o(bin_exp.size(), 0);
    // Use remapped data directly since it's already in [0, alphabet_size)
    for (int j = 0; j < sample_size - 1; j += 2) {
        uint16_t index = static_cast<uint16_t>((remapped_data[j] * alphabet_size) + remapped_data[j + 1]);
        o[e[index].bin]++;
    }

    score = calc_T(bin_exp, o);
    df = static_cast<int>(bin_exp.size()) - alphabet_size;
}

void goodness_of_fit(const uint8_t data[],
                     double &score,
                     int &df,
                     int sample_size,
                     int alphabet_size,
                     const std::unordered_map<uint8_t, int>& symbol_map)
{
    if (sample_size < 10 || alphabet_size < 2) {
        score = 0.0;
        df = 0;
        return;
    }

    // Build remapped data
    std::vector<uint8_t> remapped_data(sample_size);
    for (int i = 0; i < sample_size; i++) {
        auto it = symbol_map.find(data[i]);
        if (it == symbol_map.end()) {
            throw std::runtime_error("goodness_of_fit: unmapped symbol found");
        }
        remapped_data[i] = static_cast<uint8_t>(it->second);
    }

    std::vector<double> p(alphabet_size, 0.0);
    calc_proportions(remapped_data.data(), p, sample_size);

    std::vector<tupleTranslateEntry> e(alphabet_size);
    for (int j = 0; j < alphabet_size; j++) {
        e[j].tuple = static_cast<uint16_t>(j);
        e[j].expectation = p[j] * floor(static_cast<double>(sample_size) / 10.0);
        e[j].bin = -1;
    }

    std::sort(e.begin(), e.end(), expectationOrder);

    std::vector<double> bin_exp;
    allocate_bins(e, bin_exp);
    
    if (bin_exp.empty()) {
        score = 0.0;
        df = 0;
        return;
    }

    std::sort(e.begin(), e.end(), tupleOrder);

    int block_size = sample_size / 10;
    double T = 0.0;
    std::vector<int> o(bin_exp.size());

    for (int j = 0; j < 10; j++) {
        std::fill(o.begin(), o.end(), 0);
        // Use remapped data directly
        for (int k = 0; k < block_size; k++) {
            o[e[remapped_data[j * block_size + k]].bin]++;
        }
        T += calc_T(bin_exp, o);
    }

    score = T;
    df = 9 * (static_cast<int>(bin_exp.size()) - 1);
}

/* -----------------------------
   Binary Special Cases
----------------------------- */
void binary_chi_square_independence(const uint8_t data[],
                                    double &score,
                                    int &df,
                                    int sample_size)
{
    double p0 = 0.0, p1 = 0.0;

    for (int i = 0; i < sample_size; i++) {
        p1 += data[i];
    }

    p1 /= sample_size;
    p0 = 1.0 - p1;

    double min_p = std::min(p0, p1);
    int m = 11;
    int threshhold = 5;
    
    while (m > 1) {
        if (pow(min_p, m) * (sample_size / m) >= threshhold) {
            break;
        } else {
            m--;
        }
    }

    unsigned int tuple_count = 1 << m;

    if (m < 2) {
        score = 0.0;
        df = 0;
        return;
    }

    double T = 0;
    std::vector<int> occ(tuple_count, 0);
    int block_count = sample_size / m;

    for (int i = 0; i < block_count; i++) {
        int symbol = 0;
        for (int j = 0; j < m; j++) {
            symbol = (symbol << 1) | data[i * m + j];
        }
        occ[symbol]++;
    }

    for (unsigned int i = 0; i < occ.size(); i++) {
#ifdef _MSC_VER
        int w = __popcnt(i);
#else
        int w = __builtin_popcount(i);
#endif
        double e = pow(p1, w) * pow(p0, m - w) * block_count;
        T += pow(occ[i] - e, 2) / e;
    }

    score = T;
    df = (1 << m) - 2;
}

void binary_goodness_of_fit(const uint8_t data[],
                            double &score,
                            int &df,
                            int sample_size)
{
    int sublength = sample_size / 10;
    int ones = 0;

    for (int i = 0; i < sample_size; i++) {
        ones += data[i];
    }

    double p = divide(ones, sample_size);
    double T = 0;

    double e0 = (1.0 - p) * sublength;
    double e1 = p * sublength;

    for (int i = 0; i < 10; i++) {
        int o0 = 0, o1 = 0;

        for (int j = 0; j < sublength; j++) {
            o1 += data[i * sublength + j];
        }

        o0 = sublength - o1;

        T += (pow(o0 - e0, 2) / e0) + (pow(o1 - e1, 2) / e1);
    }

    score = T;
    df = 9;
}

namespace lib90b {

ChiSquareResult chiSquareTest(const EntropyInputData& data) {
    ChiSquareResult result{};
    if (data.symbols.empty()) {
        throw std::runtime_error("chiSquareTest: no symbols provided");
    }

    int sample_size = static_cast<int>(data.symbols.size());
    
    // Build symbol mapping: map each unique symbol to consecutive indices [0, alphabet_size)
    std::unordered_set<uint8_t> unique_symbols(data.symbols.begin(), data.symbols.end());
    std::vector<uint8_t> sorted_symbols(unique_symbols.begin(), unique_symbols.end());
    std::sort(sorted_symbols.begin(), sorted_symbols.end());
    
    std::unordered_map<uint8_t, int> symbol_map;
    for (size_t i = 0; i < sorted_symbols.size(); i++) {
        symbol_map[sorted_symbols[i]] = static_cast<int>(i);
    }
    
    int alphabet_size;
    if (data.alph_size > 0) {
        alphabet_size = data.alph_size;
    } else {
        alphabet_size = static_cast<int>(unique_symbols.size());
    }

    double score = 0.0;
    int df = 0;

    // Chi-square independence test
    if (alphabet_size == 2) {
        // Binary version works directly on data
        binary_chi_square_independence(data.symbols.data(), score, df, sample_size);
    } else {
        chi_square_independence(data.symbols.data(), score, df, sample_size, alphabet_size, symbol_map);
    }

    result.independence_score = score;
    result.independence_df = df;
    result.independence_pvalue = chi_square_pvalue(score, df);
    result.passed = result.independence_pvalue >= 0.001;

    // Reset for goodness-of-fit
    score = 0.0;
    df = 0;

    // Chi-square goodness-of-fit test
    if (alphabet_size == 2) {
        binary_goodness_of_fit(data.symbols.data(), score, df, sample_size);
    } else {
        goodness_of_fit(data.symbols.data(), score, df, sample_size, alphabet_size, symbol_map);
    }

    result.goodness_of_fit_score = score;
    result.goodness_of_fit_df = df;
    result.goodness_of_fit_pvalue = chi_square_pvalue(score, df);
    result.passed = result.passed && (result.goodness_of_fit_pvalue >= 0.001);

    return result;
}

ChiSquareResult chiSquareTest(const std::filesystem::path& filepath) {
    std::ifstream file(filepath, std::ios::binary);
    if (!file) {
        throw std::runtime_error("chiSquareTest: failed to open file " + filepath.string());
    }

    std::vector<uint8_t> buffer((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
    
    // Don't set alphabet size - let it be calculated from unique symbols
    EntropyInputData data{buffer, 8};

    return chiSquareTest(data);
}

} // namespace lib90b