// net_v5_opt.cpp
// Version 5: Row ordering + suffix enumeration + branch pruning
//
// Combines:
//   (1) Sort rows by decreasing margin (exp3)
//   (2) Enumerate suffix completions as products of choose(r_i, y_i),
//       threshold-filter, and sum (exp4)
//   (3) Prune during enumeration using suffix bounds:
//       - All-below: if max possible factor for a branch <= threshold,
//         bulk-sum via Vandermonde identity instead of enumerating
//       - All-above: if min possible factor > threshold, skip branch
//
// The threshold is tau = p_obs*(1+tol)/D where D = c1!*c2!/n!.
// (tau, not t: the manuscripts reserve t for the rxc/rxck path index.)
// P(table) = D * prod(choose(r_i, y_i)), so factor = prod(choose(r_i, y_i)).
//
// prefix_factor tracks the accumulated product of choose values for rows
// already assigned. This lets us apply pruning at each branch point:
//   full_factor = prefix_factor * suffix_factor
// We need full_factor <= t for the table to contribute to the p-value.
//
// No caching: pruned enumeration depends on prefix_factor, so (k, c1)
// subproblems are not reusable across different prefixes. The pruning
// savings (skipping or bulk-summing entire subtrees) should compensate.

#include <Rcpp.h>
#include <cmath>
#include <algorithm>
#include <vector>

using namespace Rcpp;

#define MAX_ROWS 20

namespace {

struct FisherState {
    int m;
    int n;
    int r[MAX_ROWS];
    int c[2];
    double tol;
    int suffix_r[MAX_ROWS + 1];
    double suffix_max[MAX_ROWS + 1];
    double threshold;
};

// Pruned enumeration for rows k..m-1 with c1 remaining.
// prefix_factor = product of choose(r[i], y[i]) for rows 0..k-1.
// Returns the sum of suffix factors (rows k..m-1) for all completions
// whose FULL factor (prefix * suffix) is <= threshold.
static double enumerate_pruned(int k, int c1, double prefix_factor,
                               FisherState& s) {

    if (k >= s.m) {
        return (prefix_factor <= s.threshold) ? 1.0 : 0.0;
    }

    if (k == s.m - 1) {
        if (c1 < 0 || c1 > s.r[k]) return 0.0;
        double f = R::choose(s.r[k], c1);
        return (prefix_factor * f <= s.threshold) ? f : 0.0;
    }

    int y_lo = std::max(0, c1 - s.suffix_r[k + 1]);
    int y_hi = std::min(s.r[k], c1);
    if (y_lo > y_hi) return 0.0;

    // All-below: every factor from rows k..m-1 is at most suffix_max[k].
    // If prefix * suffix_max[k] <= threshold, every completion qualifies.
    // Sum of all choose products = choose(suffix_r[k], c1) by Vandermonde.
    if (prefix_factor * s.suffix_max[k] <= s.threshold) {
        return R::choose(s.suffix_r[k], c1);
    }

    // All-above: every factor from rows k..m-1 is at least 1 (choose(r,0)=1).
    // If prefix * 1 > threshold, no completion qualifies.
    if (prefix_factor > s.threshold) {
        return 0.0;
    }

    double sum = 0.0;

    for (int y_k = y_lo; y_k <= y_hi; y_k++) {
        double factor_k = R::choose(s.r[k], y_k);
        double new_prefix = prefix_factor * factor_k;

        // Per-branch all-below: if new_prefix * suffix_max[k+1] <= threshold,
        // every completion through this y_k qualifies.
        if (new_prefix * s.suffix_max[k + 1] <= s.threshold) {
            sum += factor_k * R::choose(s.suffix_r[k + 1], c1 - y_k);
            continue;
        }

        // Per-branch all-above: if new_prefix > threshold, even the minimum
        // suffix factor (= 1) makes the full factor exceed threshold.
        if (new_prefix > s.threshold) {
            continue;
        }

        sum += factor_k * enumerate_pruned(k + 1, c1 - y_k, new_prefix, s);
    }

    return sum;
}

}  // anonymous namespace

// [[Rcpp::export(name = ".rx2_net_vander_cpp")]]
double rx2_net_vander_cpp(IntegerMatrix dat) {
    FisherState s;
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

    s.suffix_max[s.m] = 1.0;
    for (int i = s.m - 1; i >= 0; i--)
        s.suffix_max[i] = s.suffix_max[i + 1] * R::choose(s.r[i], s.r[i] / 2);

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
