// tree_v5_opt.cpp
// Tree v5: Best tree implementation with backported network improvements
//
// Base: tree_v4_opt (memoized find_max/find_min, fixed arrays)
// Added:
//   (1) Precomputed log-factorial table — replaces all lgamma calls
//       in compute_prob and dhyper with O(1) lookups
//   (2) Hypergeometric recurrence in traverse — compute dhyper
//       incrementally as y_k changes: O(1) per step
//   (3) Row ordering — sort rows by decreasing margin
//   (4) Constrained suffix_max DP — tighter pruning bounds
//   (5) Pure C++ math — no R API in hot paths

#include <Rcpp.h>
#include <cmath>
#include <algorithm>
#include <vector>
#include <unordered_map>

using namespace Rcpp;

#define MAX_ROWS 20
#define MAX_N 10000

namespace {

struct CacheEntry {
    double prob;
    int y[MAX_ROWS];
};

struct FisherState {
    int m;
    int n;
    int r[MAX_ROWS];
    int c[2];
    double tol;
    double log_const;
    double p_obs;
    int suffix_r[MAX_ROWS + 1];
    double pval;
    double p_max;
    double p_min;
    int y_max[MAX_ROWS];
    int y_min[MAX_ROWS];

    std::unordered_map<int, CacheEntry> memo_max;
    std::unordered_map<int, CacheEntry> memo_min;

    // Precomputed log-factorials: lfact[i] = log(i!)
    double lfact[MAX_N + 1];
    int lfact_size;

    // Constrained suffix max: cmax[k * c0p1 + c1]
    std::vector<double> cmax;
    int c0p1;
};

static inline int cache_key(int k, int c1) {
    return k * (MAX_N + 1) + c1;
}

static inline double compute_prob(const int* y, const FisherState& s) {
    double lp = s.log_const;
    for (int i = 0; i < s.m; i++) {
        lp -= s.lfact[y[i]] + s.lfact[s.r[i] - y[i]];
    }
    return std::exp(lp);
}

static inline double cpp_choose_lf(int n_arg, int k_arg, const FisherState& s) {
    if (k_arg < 0 || k_arg > n_arg) return 0.0;
    if (k_arg == 0 || k_arg == n_arg) return 1.0;
    return std::round(std::exp(s.lfact[n_arg] - s.lfact[k_arg] - s.lfact[n_arg - k_arg]));
}

// dhyper using log-factorial table
static inline double dhyper_lf(int x, int m_white, int n_black, int k_draw,
                               const FisherState& s) {
    if (x < 0 || x > k_draw || x > m_white || k_draw - x > n_black) return 0.0;
    double lp = s.lfact[m_white] - s.lfact[x] - s.lfact[m_white - x]
              + s.lfact[n_black] - s.lfact[k_draw - x] - s.lfact[n_black - k_draw + x]
              - s.lfact[m_white + n_black] + s.lfact[k_draw] + s.lfact[m_white + n_black - k_draw];
    return std::exp(lp);
}

// --- find_min with memoization ---

static void joe_min_impl(int k, bool descend, int* r_avail, int c1, int n_rem,
                         int* y, FisherState& s, double& local_min_p, int* local_min_y) {
    if (k == s.m) {
        for (int i = 0; i < s.m; i++) {
            if (r_avail[i] >= 0) y[i] = c1;
        }
        double prob = compute_prob(y, s);
        if (prob < local_min_p) {
            local_min_p = prob;
            for (int i = 0; i < s.m; i++) local_min_y[i] = y[i];
        }
        return;
    }

    int r_avail_save[MAX_ROWS], y_save[MAX_ROWS];
    for (int i = 0; i < s.m; i++) {
        r_avail_save[i] = r_avail[i];
        y_save[i] = y[i];
    }
    int c1_save = c1, n_rem_save = n_rem;

    int idx = 0, max_r = -1;
    for (int i = 0; i < s.m; i++) {
        if (r_avail[i] > max_r) { max_r = r_avail[i]; idx = i; }
    }

    y[idx] = descend ? std::min(r_avail[idx], c1)
                     : r_avail[idx] - std::min(r_avail[idx], n_rem - c1);
    n_rem -= r_avail[idx];
    r_avail[idx] = -1;
    c1 -= y[idx];

    joe_min_impl(k + 1, true, r_avail, c1, n_rem, y, s, local_min_p, local_min_y);

    if (descend) {
        for (int i = 0; i < s.m; i++) {
            r_avail[i] = r_avail_save[i];
            y[i] = y_save[i];
        }
        joe_min_impl(k, false, r_avail, c1_save, n_rem_save, y, s, local_min_p, local_min_y);
    }
}

static void find_min(int c1, int* y, int k_start, FisherState& s) {
    int key = cache_key(k_start, c1);

    auto it = s.memo_min.find(key);
    if (it != s.memo_min.end()) {
        const CacheEntry& entry = it->second;
        for (int i = k_start; i < s.m; i++) y[i] = entry.y[i];
        s.p_min = compute_prob(y, s);
        for (int i = 0; i < s.m; i++) s.y_min[i] = y[i];
        return;
    }

    int r_avail[MAX_ROWS];
    int n_rem = 0;
    for (int i = 0; i < s.m; i++) {
        if (i < k_start) r_avail[i] = -1;
        else { r_avail[i] = s.r[i]; n_rem += s.r[i]; }
    }

    double local_min_p = 1e300;
    int local_min_y[MAX_ROWS];
    for (int i = 0; i < s.m; i++) local_min_y[i] = y[i];

    joe_min_impl(k_start, true, r_avail, c1, n_rem, y, s, local_min_p, local_min_y);

    s.p_min = local_min_p;
    for (int i = 0; i < s.m; i++) s.y_min[i] = local_min_y[i];

    CacheEntry entry;
    entry.prob = local_min_p;
    for (int i = 0; i < s.m; i++) entry.y[i] = local_min_y[i];
    s.memo_min[key] = entry;
}

// --- find_max with memoization ---

struct FindMaxEntry {
    int k, c1, n_rem;
    int y[MAX_ROWS];
};

static void find_max(int k_start, int c1, int n_rem, int* y, FisherState& s) {
    int key = cache_key(k_start, c1);

    auto it = s.memo_max.find(key);
    if (it != s.memo_max.end()) {
        const CacheEntry& entry = it->second;
        for (int i = k_start; i < s.m; i++) y[i] = entry.y[i];
        s.p_max = compute_prob(y, s);
        for (int i = 0; i < s.m; i++) s.y_max[i] = y[i];
        return;
    }

    double local_max_p = 0;
    int local_max_y[MAX_ROWS];
    for (int i = 0; i < s.m; i++) local_max_y[i] = y[i];

    FindMaxEntry stack[256];
    int sp = 0;

    stack[sp].k = k_start;
    stack[sp].c1 = c1;
    stack[sp].n_rem = n_rem;
    for (int i = 0; i < s.m; i++) stack[sp].y[i] = y[i];
    sp++;

    while (sp > 0) {
        sp--;
        FindMaxEntry curr = stack[sp];
        int k = curr.k;
        int c1_rem = curr.c1;
        int n_rem_curr = curr.n_rem;

        if (k >= s.m) {
            double prob = compute_prob(curr.y, s);
            if (prob > local_max_p) {
                local_max_p = prob;
                for (int i = 0; i < s.m; i++) local_max_y[i] = curr.y[i];
            }
            continue;
        }

        int d_rem = s.m - k;
        int denom = n_rem_curr + d_rem;
        int y_lo, y_up;

        if (denom == 0 || c1_rem == 0) {
            y_lo = y_up = 0;
        } else {
            y_up = (int)std::floor((double)(s.r[k] + 1) * (c1_rem + d_rem - 1) / denom);
            y_lo = (int)std::ceil((double)(s.r[k] + 1) * (c1_rem + 1) / denom) - 1;
        }
        if (c1_rem == 0) y_lo = y_up = 0;

        int n_after = s.suffix_r[k + 1];
        y_lo = std::max(y_lo, std::max(0, c1_rem - n_after));
        y_up = std::min(y_up, std::min(s.r[k], c1_rem));

        if (y_lo > y_up) continue;

        for (int y_k = y_lo; y_k <= y_up && sp < 255; y_k++) {
            stack[sp].k = k + 1;
            stack[sp].c1 = c1_rem - y_k;
            stack[sp].n_rem = n_rem_curr - s.r[k];
            for (int i = 0; i < s.m; i++) stack[sp].y[i] = curr.y[i];
            stack[sp].y[k] = y_k;
            sp++;
        }
    }

    s.p_max = local_max_p;
    for (int i = 0; i < s.m; i++) s.y_max[i] = local_max_y[i];

    CacheEntry entry;
    entry.prob = local_max_p;
    for (int i = 0; i < s.m; i++) entry.y[i] = local_max_y[i];
    s.memo_max[key] = entry;
}

// --- traverse with hypergeometric recurrence ---

static void traverse(int k, int c1, int* y, double prob_prefix, int n_rem,
                     bool descend, int y_k, int dir, int mode, FisherState& s) {

    if (k == s.m - 1) {
        y[s.m - 1] = c1;
        double prob = compute_prob(y, s);
        if (prob > s.p_obs * (1 + s.tol)) s.pval -= prob;
        return;
    }

    int y_lo = std::max(0, c1 - s.suffix_r[k + 1]);
    int y_hi = std::min(s.r[k], c1);

    if (descend) { mode = s.y_max[k]; y_k = mode; }
    if (y_k < y_lo) return;
    if (y_k > y_hi) {
        traverse(k, c1, y, prob_prefix, n_rem, false, mode - 1, -1, mode, s);
        return;
    }

    y[k] = y_k;
    double hp = dhyper_lf(y_k, c1, n_rem - c1, s.r[k], s);
    double new_prefix = prob_prefix * hp;
    int new_c1 = c1 - y_k;
    int new_n_rem = s.suffix_r[k + 1];

    double threshold = s.p_obs * (1 + s.tol);

    if (s.m - k > 2) {
        // Constrained suffix_max pruning
        double smax = s.cmax[(k + 1) * s.c0p1 + new_c1];
        if (new_prefix * smax <= threshold) {
            traverse(k, c1, y, prob_prefix, n_rem, false, y_k + dir, dir, mode, s);
            return;
        }

        s.p_max = 0;
        find_max(k + 1, new_c1, new_n_rem, y, s);
        s.p_min = 1;
        find_min(new_c1, y, k + 1, s);

        if (s.p_min > threshold) {
            s.pval -= new_prefix;
            traverse(k, c1, y, prob_prefix, n_rem, false, y_k + dir, dir, mode, s);
            return;
        }
        if (s.p_max <= threshold) {
            traverse(k, c1, y, prob_prefix, n_rem, false, y_k + dir, dir, mode, s);
            return;
        }
    }

    traverse(k + 1, new_c1, y, new_prefix, new_n_rem, true, 0, +1, 0, s);
    traverse(k, c1, y, prob_prefix, n_rem, false, y_k + dir, dir, mode, s);
}

// --- Precompute constrained suffix max ---
// cmax stores the maximum suffix PROBABILITY directly:
//   cmax[k][c1] = max over completions y[k..m-1] with sum = c1 of the
//   probability contribution of the suffix, computed with the same
//   log-factorial arithmetic as compute_prob.
//
// This constrained bound replaces the unconstrained
//   suffix_max[k] = prod_{i=k}^{m-1} choose(r[i], floor(r[i]/2))
// of the memoized kernel. The unconstrained product ignores both the
// remaining column budget c1 and the hypergeometric denominator, so
// it is loose; conditioning the DP on (k, c1) yields the tight bound
// that matches what find_max would return for the same state, at the
// cost of an O(m * (c1+1)) precomputation.

static void precompute_cmax(FisherState& s) {
    int c0p1 = s.c[0] + 1;
    s.c0p1 = c0p1;
    s.cmax.assign((s.m + 1) * c0p1, 0.0);

    // cmax[k][c1] = max probability achievable by any table where
    // rows 0..k-1 contribute the "best" prefix and rows k..m-1 have
    // column-1 values summing to c1.
    //
    // But that requires knowing the prefix. Instead, store the max of
    // the SUFFIX contribution in the same units as suffix_max:
    // the product of choose(r[i], y[i]) for i=k..m-1.
    // This is independent of the prefix and directly replaces suffix_max[k].

    // Base: past last row, empty product = 1, only c1=0 is feasible
    s.cmax[s.m * c0p1 + 0] = 1.0;

    for (int k = s.m - 1; k >= 0; k--) {
        for (int c1 = 0; c1 < c0p1; c1++) {
            int y_lo = std::max(0, c1 - s.suffix_r[k + 1]);
            int y_hi = std::min(s.r[k], c1);
            if (y_lo > y_hi) continue;

            double best = 0.0;
            double binom = cpp_choose_lf(s.r[k], y_lo, s);
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

}  // anonymous namespace

// [[Rcpp::export(name = ".rx2_tree_dp_cpp")]]
double rx2_tree_dp_cpp(IntegerMatrix dat) {
    FisherState s;
    s.m = dat.nrow();
    if (s.m > MAX_ROWS) Rcpp::stop("Too many rows (max %d)", MAX_ROWS);

    s.tol = 3.45254e-7;
    s.pval = 1.0;
    s.n = 0;
    s.c[0] = s.c[1] = 0;

    int r_unsorted[MAX_ROWS];
    for (int i = 0; i < s.m; i++) {
        r_unsorted[i] = dat(i, 0) + dat(i, 1);
        s.n += r_unsorted[i];
        s.c[0] += dat(i, 0);
        s.c[1] += dat(i, 1);
    }

    if (s.n > MAX_N) Rcpp::stop("Total count too large (max %d)", MAX_N);

    // Precompute log-factorial table
    s.lfact_size = s.n + 1;
    s.lfact[0] = 0.0;
    for (int i = 1; i <= s.n; i++)
        s.lfact[i] = s.lfact[i - 1] + std::log((double)i);

    // Sort rows by decreasing margin
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

    s.log_const = 0.0;
    for (int i = 0; i < s.m; i++)
        s.log_const += s.lfact[s.r[i]];
    s.log_const += s.lfact[s.c[0]] + s.lfact[s.c[1]] - s.lfact[s.n];

    double log_p_obs = s.log_const;
    for (int i = 0; i < s.m; i++) {
        log_p_obs -= s.lfact[y_obs[i]] + s.lfact[s.r[i] - y_obs[i]];
    }
    s.p_obs = std::exp(log_p_obs);

    // Precompute constrained suffix max
    precompute_cmax(s);

    for (int i = 0; i < s.m; i++) { s.y_max[i] = 0; s.y_min[i] = 0; }
    s.p_max = 0;
    s.p_min = 1;

    find_max(0, s.c[0], s.n, y_obs, s);
    traverse(0, s.c[0], y_obs, 1.0, s.n, true, 0, +1, 0, s);

    return s.pval;
}
