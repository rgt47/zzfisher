// net_v6_opt.cpp
// Version 6: v5 + constrained suffix_max DP + hypergeometric recurrence
//
// Over v5, adds:
//   (4) Constrained suffix_max: precompute exact maximum of
//       prod(choose(r_i, y_i)) subject to sum(y_i)=c1 via DP.
//       Tighter than v5's unconstrained product of individual maxima.
//   (5) Hypergeometric recurrence: compute choose(r_k, y_k)
//       incrementally via ratio (r_k - y_k)/(y_k + 1).
//
// Pure C++ math (std::lgamma) throughout — no R API in hot paths.

#include <Rcpp.h>
#include <cmath>
#include <algorithm>
#include <vector>

using namespace Rcpp;

#define MAX_ROWS 20

static inline double cpp_choose(int n, int k) {
    if (k < 0 || k > n) return 0.0;
    if (k == 0 || k == n) return 1.0;
    return std::round(std::exp(
        std::lgamma(n + 1) - std::lgamma(k + 1) - std::lgamma(n - k + 1)
    ));
}

struct FisherStateV6 {
    int m;
    int n;
    int r[MAX_ROWS];
    int c[2];
    double tol;
    int suffix_r[MAX_ROWS + 1];
    double threshold;

    std::vector<double> cmax;
    int c0p1;
};

static inline double get_cmax(int k, int c1, const FisherStateV6& s) {
    return s.cmax[k * s.c0p1 + c1];
}

static double enumerate_pruned(int k, int c1, double prefix_factor,
                               FisherStateV6& s) {

    if (k >= s.m) {
        return (prefix_factor <= s.threshold) ? 1.0 : 0.0;
    }

    if (k == s.m - 1) {
        if (c1 < 0 || c1 > s.r[k]) return 0.0;
        double f = cpp_choose(s.r[k], c1);
        return (prefix_factor * f <= s.threshold) ? f : 0.0;
    }

    int y_lo = std::max(0, c1 - s.suffix_r[k + 1]);
    int y_hi = std::min(s.r[k], c1);
    if (y_lo > y_hi) return 0.0;

    if (prefix_factor * get_cmax(k, c1, s) <= s.threshold) {
        return cpp_choose(s.suffix_r[k], c1);
    }

    if (prefix_factor > s.threshold) {
        return 0.0;
    }

    double sum = 0.0;
    double factor_k = cpp_choose(s.r[k], y_lo);
    int rk = s.r[k];

    for (int y_k = y_lo; y_k <= y_hi; y_k++) {
        double new_prefix = prefix_factor * factor_k;
        int new_c1 = c1 - y_k;

        if (new_prefix * get_cmax(k + 1, new_c1, s) <= s.threshold) {
            sum += factor_k * cpp_choose(s.suffix_r[k + 1], new_c1);
        } else if (new_prefix <= s.threshold) {
            sum += factor_k * enumerate_pruned(k + 1, new_c1, new_prefix, s);
        }

        if (y_k < y_hi) {
            factor_k = factor_k * (rk - y_k) / (y_k + 1);
        }
    }

    return sum;
}

static void precompute_cmax(FisherStateV6& s) {
    int c0p1 = s.c[0] + 1;
    s.c0p1 = c0p1;
    s.cmax.assign((s.m + 1) * c0p1, 0.0);

    s.cmax[s.m * c0p1 + 0] = 1.0;

    for (int k = s.m - 1; k >= 0; k--) {
        for (int c1 = 0; c1 < c0p1; c1++) {
            int y_lo = std::max(0, c1 - s.suffix_r[k + 1]);
            int y_hi = std::min(s.r[k], c1);
            if (y_lo > y_hi) continue;

            double best = 0.0;
            double binom = cpp_choose(s.r[k], y_lo);
            int rk = s.r[k];

            for (int y_k = y_lo; y_k <= y_hi; y_k++) {
                double val = binom * s.cmax[(k + 1) * c0p1 + (c1 - y_k)];
                if (val > best) best = val;

                if (y_k < y_hi) {
                    binom = binom * (rk - y_k) / (y_k + 1);
                }
            }
            s.cmax[k * c0p1 + c1] = best;
        }
    }
}

// [[Rcpp::export(name = ".net_dp_cpp")]]
double net_dp_cpp(IntegerMatrix dat) {
    FisherStateV6 s;
    s.m = dat.nrow();
    if (s.m > MAX_ROWS) Rcpp::stop("Too many rows (max %d)", MAX_ROWS);

    s.tol = 3.45254e-7;
    s.n = 0;
    s.c[0] = s.c[1] = 0;

    int r_unsorted[MAX_ROWS];
    for (int i = 0; i < s.m; i++) {
        r_unsorted[i] = dat(i, 0) + dat(i, 1);
        s.n += r_unsorted[i];
        s.c[0] += dat(i, 0);
        s.c[1] += dat(i, 1);
    }

    int order[MAX_ROWS];
    for (int i = 0; i < s.m; i++) order[i] = i;
    std::sort(order, order + s.m, [&](int a, int b) {
        return r_unsorted[a] > r_unsorted[b];
    });

    int y_obs[MAX_ROWS];
    bool flipped = s.c[0] > s.c[1];
    if (flipped) std::swap(s.c[0], s.c[1]);

    for (int i = 0; i < s.m; i++) {
        int orig = order[i];
        s.r[i] = r_unsorted[orig];
        y_obs[i] = flipped ? dat(orig, 1) : dat(orig, 0);
    }

    s.suffix_r[s.m] = 0;
    for (int i = s.m - 1; i >= 0; i--)
        s.suffix_r[i] = s.suffix_r[i + 1] + s.r[i];

    precompute_cmax(s);

    double log_D = std::lgamma(s.c[0] + 1) + std::lgamma(s.c[1] + 1)
                 - std::lgamma(s.n + 1);
    double D = std::exp(log_D);

    double log_p_obs = log_D;
    for (int i = 0; i < s.m; i++) {
        log_p_obs += std::lgamma(s.r[i] + 1)
                   - std::lgamma(y_obs[i] + 1)
                   - std::lgamma(s.r[i] - y_obs[i] + 1);
    }
    double p_obs = std::exp(log_p_obs);
    s.threshold = p_obs * (1 + s.tol) / D;

    double factor_sum = enumerate_pruned(0, s.c[0], 1.0, s);

    return D * factor_sum;
}
