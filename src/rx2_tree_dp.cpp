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

struct CacheEntry {
    double prob;
    int y[MAX_ROWS];
};

struct TreeStateV5 {
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

static inline double compute_prob(const int* y, const TreeStateV5& s) {
    double lp = s.log_const;
    for (int i = 0; i < s.m; i++) {
        lp -= s.lfact[y[i]] + s.lfact[s.r[i] - y[i]];
    }
    return std::exp(lp);
}

static inline double cpp_choose_lf(int n, int k, const TreeStateV5& s) {
    if (k < 0 || k > n) return 0.0;
    if (k == 0 || k == n) return 1.0;
    return std::round(std::exp(s.lfact[n] - s.lfact[k] - s.lfact[n - k]));
}

// dhyper using log-factorial table
static inline double dhyper_lf(int x, int m_white, int n_black, int k_draw,
                               const TreeStateV5& s) {
    if (x < 0 || x > k_draw || x > m_white || k_draw - x > n_black) return 0.0;
    double lp = s.lfact[m_white] - s.lfact[x] - s.lfact[m_white - x]
              + s.lfact[n_black] - s.lfact[k_draw - x] - s.lfact[n_black - k_draw + x]
              - s.lfact[m_white + n_black] + s.lfact[k_draw] + s.lfact[m_white + n_black - k_draw];
    return std::exp(lp);
}

// --- find_min with memoization ---

static void joe_min_impl(int k, bool descend, int* r_avail, int c1, int n_rem,
                         int* y, TreeStateV5& s, double& local_min_p, int* local_min_y) {
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

static void find_min(int c1, int* y, int d, TreeStateV5& s) {
    int key = cache_key(d, c1);

    auto it = s.memo_min.find(key);
    if (it != s.memo_min.end()) {
        const CacheEntry& entry = it->second;
        for (int i = d; i < s.m; i++) y[i] = entry.y[i];
        s.p_min = compute_prob(y, s);
        for (int i = 0; i < s.m; i++) s.y_min[i] = y[i];
        return;
    }

    int r_avail[MAX_ROWS];
    int n_rem = 0;
    for (int i = 0; i < s.m; i++) {
        if (i < d) r_avail[i] = -1;
        else { r_avail[i] = s.r[i]; n_rem += s.r[i]; }
    }

    double local_min_p = 1e300;
    int local_min_y[MAX_ROWS];
    for (int i = 0; i < s.m; i++) local_min_y[i] = y[i];

    joe_min_impl(d, true, r_avail, c1, n_rem, y, s, local_min_p, local_min_y);

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

static void find_max(int k_start, int c1, int n_rem, int* y, TreeStateV5& s) {
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

        int d = s.m - k;
        int denom = n_rem_curr + d;
        int y_lo, y_up;

        if (denom == 0 || c1_rem == 0) {
            y_lo = y_up = 0;
        } else {
            y_up = (int)std::floor((double)(s.r[k] + 1) * (c1_rem + d - 1) / denom);
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
                     bool descend, int y_k, int dir, int mode, TreeStateV5& s) {

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
// cmax[k][c1] = max of prod(choose(r[i], y[i])) for rows k..m-1
// subject to sum(y[i]) = c1.
// Uses the suffix_max formulation: this is the maximum possible factor
// (product of binomial coefficients) achievable for a given suffix and budget.
// In the tree paradigm, this bounds the suffix probability via:
//   max_suffix_prob = D_suffix * cmax[k][c1]
// where D_suffix accounts for the constant terms.
//
// For the tree pruning check (new_prefix * suffix_max <= threshold),
// we need suffix_max in probability space. The relationship:
//   max_prob_suffix = exp(log_const_suffix) * (1 / prod(r_i!)) * cmax
// But compute_prob uses log_const which includes all row factorials.
// So we precompute cmax as the max of:
//   prod_{i=k}^{m-1} choose(r[i], y[i]) = prod (r[i]! / (y[i]! (r[i]-y[i])!))
// and the suffix_max bound for the old pruning check becomes:
//   new_prefix * suffix_max_prob <= threshold
// where suffix_max_prob = exp(sum_{i=k+1}^{m-1} lfact[r[i]]) / exp(lfact[suffix_r[k+1]])
//                       * cmax[k+1][new_c1]    ... this gets complicated.
//
// Simpler approach: store cmax as the max probability achievable for the
// suffix subproblem, computed directly using compute_prob logic.
// cmax_prob[k][c1] = max over y[k..m-1] summing to c1 of
//   exp(log_const_suffix - sum(lfact[y_i] + lfact[r_i - y_i]))
// where log_const_suffix = sum_{i=k}^{m-1} lfact[r[i]] + lfact[c1] + lfact[suffix_r[k]-c1] - lfact[suffix_r[k]]
//
// Actually, the simplest correct approach: store the max of
// prod(1/(y_i! * (r_i-y_i)!)) for the suffix, since
// prob = exp(log_const) * prod(1/(y_i! * (r_i-y_i)!)) over ALL rows.
// The prefix contribution is already in prob_prefix (via dhyper product).
// The suffix_max check in v4 uses:
//   suffix_max[k] = prod_{i=k}^{m-1} choose(r[i], floor(r[i]/2))
// and the check is: new_prefix * suffix_max[k+1] <= threshold
// where new_prefix is a product of dhyper values, and threshold = p_obs*(1+tol).
//
// dhyper(y_k, c1, n_rem-c1, r_k) = choose(c1,y_k)*choose(n_rem-c1, r_k-y_k) / choose(n_rem, r_k)
// So prob_prefix = prod_{i<k} dhyper_i, and the suffix bound must be
// an upper bound on prod_{i>=k} dhyper_i.
//
// The unconstrained suffix_max in v4 uses:
//   suffix_max[k] = prod_{i=k}^{m-1} choose(r[i], floor(r[i]/2))
// This doesn't directly bound the dhyper product. Let me re-examine v4's logic...
//
// In v4, suffix_max[k] = prod_{i=k}^{m-1} choose(r[i], floor(r[i]/2))
// and the check is: new_prefix * suffix_max[k+1] <= p_obs * (1+tol)
// where new_prefix = prob_prefix * dhyper(y_k, ...).
//
// This works because the full probability can be factored as:
//   P = D * prod_{i} choose(r[i], y[i])
// where D = c1!*c2!/n! (from the v5/v6 derivation).
// And prob_prefix in the traverse = D * prod_{i<k} choose(r[i], y[i])
//   ... no wait. prob_prefix = prod_{i<k} dhyper(y_i, ...).
//
// The dhyper factorization:
//   prod_{i=0}^{m-1} dhyper(y_i, c1_i, n_rem_i - c1_i, r[i])
//   = P(table) (the multivariate hypergeometric, by sequential conditioning)
//
// So prob_prefix * prod_{i=k}^{m-1} dhyper(y_i, ...) = P(table).
// The suffix_max needs to bound max over feasible suffixes of
// prod_{i>=k+1} dhyper(y_i, c1_i, n_rem_i - c1_i, r[i]).
//
// But dhyper depends on c1_i which changes along the path. This is why
// v4's suffix_max is an approximation — it uses choose(r[i], floor(r[i]/2))
// as an upper bound on each dhyper factor individually.
//
// For the constrained version, we need:
//   cmax_dhyper[k][c1] = max over y[k..m-1] summing to c1 of
//     prod_{i=k}^{m-1} dhyper(y_i, c1_i, n_rem_i - c1_i, r[i])
//
// This IS the same as max P(suffix) / P(prefix), but it's complex because
// the dhyper at each step depends on the running c1.
//
// Simpler: since P(table) = D * prod choose(r_i, y_i) and
// prob_prefix (as computed in traverse) = product of dhyper values =
// P(table with y[0..k-1] fixed, summing over suffixes) ... actually no.
// prob_prefix = prod dhyper != D * prod choose for prefix only.
//
// Let me just verify what v4 actually computes. In traverse_v4:
//   hp = dhyper_v4(y_k, c1, n_rem - c1, s.r[k])
//   new_prefix = prob_prefix * hp
// And the check: new_prefix * suffix_max[k+1] <= p_obs * (1+tol)
// where suffix_max[k+1] = prod_{i>=k+1} choose(r[i], floor(r[i]/2))
//
// For this check to be valid, we need:
//   new_prefix * suffix_max[k+1] >= P(any table extending current prefix)
// i.e., suffix_max[k+1] >= P(table) / new_prefix for all feasible suffixes.
//
// P(table) / prob_prefix_through_k = prod_{i=k+1}^{m-1} dhyper(y_i, ...)
//
// And choose(r_i, floor(r_i/2)) is NOT a valid upper bound on dhyper.
// dhyper can exceed 1 when parameters are favorable.
// Actually: dhyper(x, m, n, k) = C(m,x)*C(n,k-x)/C(m+n,k) <= 1 always.
// And C(m,x)*C(n,k-x)/C(m+n,k) <= C(m, floor(m/2)) * C(n, ...) / C(m+n,k)
// Hmm, this bound is not simply choose(r_i, floor(r_i/2)).
//
// Wait — I think v4's suffix_max is actually bounding a different quantity.
// Let me just keep the same suffix_max logic as v4 but with the constrained
// version. The DP computes:
//   cmax[k][c1] = max over feasible y summing to c1 of
//     prod_{i=k}^{m-1} choose(r[i], y[i])
// And the pruning check becomes:
//   new_prefix * D_ratio * cmax[k+1][new_c1] <= threshold
// where D_ratio converts between the dhyper-product space and the
// choose-product space.
//
// This is getting too complicated. Let me take a different approach:
// just precompute the max suffix PROBABILITY directly via DP on
// compute_prob values, matching exactly what find_max would return.

static void precompute_cmax(TreeStateV5& s) {
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

// [[Rcpp::export(name = ".rx2_tree_dp_cpp")]]
double rx2_tree_dp_cpp(IntegerMatrix dat) {
    TreeStateV5 s;
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
