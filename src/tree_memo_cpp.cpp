// tree_v4_opt.cpp
// Version 4: Memoization of find_max and find_min, log-factorial table
// Built on v3_opt's fixed-size array approach
//
// Key insight: find_max(k, c1) and find_min(k, c1) results depend only on
// k (starting row), c1 (remaining column 1 total), and r[k:m-1] (fixed).
// By caching these results, we avoid redundant computation.

#include <Rcpp.h>
#include <cmath>
#include <algorithm>
#include <unordered_map>
#include <cstring>

#ifdef PROFILE_V4
#include <chrono>
#endif

using namespace Rcpp;

#ifdef PROFILE_V4
struct ProfileCounters {
    double time_traverse;
    double time_find_max;
    double time_find_min;
    double time_compute_prob;
    long long calls_traverse;
    long long calls_find_max;
    long long calls_find_min;
    long long calls_compute_prob;
    long long cache_hits_max;
    long long cache_hits_min;
    long long cache_misses_max;
    long long cache_misses_min;
    long long nodes_pruned_suffix;
    long long nodes_pruned_max;
    long long nodes_pruned_min;

    ProfileCounters() : time_traverse(0), time_find_max(0),
        time_find_min(0), time_compute_prob(0),
        calls_traverse(0), calls_find_max(0), calls_find_min(0),
        calls_compute_prob(0), cache_hits_max(0), cache_hits_min(0),
        cache_misses_max(0), cache_misses_min(0),
        nodes_pruned_suffix(0), nodes_pruned_max(0),
        nodes_pruned_min(0) {}
};

static ProfileCounters* g_prof = nullptr;
#endif

#define MAX_ROWS 20
#define MAX_N 10001  // max total count for log-factorial table

struct CacheEntry {
    int y[MAX_ROWS];  // Only indices k to m-1 are valid
};

struct FisherStateV4 {
    int m;
    int n;
    int r[MAX_ROWS];
    int c[2];
    double tol;
    double log_const;
    double p_obs;
    double suffix_max[MAX_ROWS + 1];
    int suffix_r[MAX_ROWS + 1];
    double lfact[MAX_N + 1];  // lfact[i] = log(i!)
    double pval;
    double p_max;
    double p_min;
    int y_max[MAX_ROWS];

    // Memoization caches: key = k * 10001 + c1
    std::unordered_map<int, CacheEntry> memo_max;
    std::unordered_map<int, CacheEntry> memo_min;
};

inline int cache_key(int k, int c1) {
    return k * 10001 + c1;
}

inline double compute_prob_v4(const int* y, const FisherStateV4& s) {
#ifdef PROFILE_V4
    if (g_prof) g_prof->calls_compute_prob++;
#endif
    double log_prob = s.log_const;
    for (int i = 0; i < s.m; i++) {
        log_prob -= s.lfact[y[i]];
        log_prob -= s.lfact[s.r[i] - y[i]];
    }
    return std::exp(log_prob);
}

// Forward declaration
void joe_min_v4_impl(int k, bool descend, int* r_avail, int c1, int n_rem,
                     int* y, FisherStateV4& s, double& local_min_p, int* local_min_y);

void joe_min_v4_impl(int k, bool descend, int* r_avail, int c1, int n_rem,
                     int* y, FisherStateV4& s, double& local_min_p, int* local_min_y) {
    if (k == s.m) {
        for (int i = 0; i < s.m; i++) {
            if (r_avail[i] >= 0) y[i] = c1;
        }
        double prob = compute_prob_v4(y, s);
        if (prob < local_min_p) {
            local_min_p = prob;
            std::memcpy(local_min_y, y, s.m * sizeof(int));
        }
        return;
    }

    int r_avail_save[MAX_ROWS], y_save[MAX_ROWS];
    std::memcpy(r_avail_save, r_avail, s.m * sizeof(int));
    std::memcpy(y_save, y, s.m * sizeof(int));
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

    joe_min_v4_impl(k + 1, true, r_avail, c1, n_rem, y, s, local_min_p, local_min_y);

    if (descend) {
        std::memcpy(r_avail, r_avail_save, s.m * sizeof(int));
        std::memcpy(y, y_save, s.m * sizeof(int));
        joe_min_v4_impl(k, false, r_avail, c1_save, n_rem_save, y, s, local_min_p, local_min_y);
    }
}

void find_min_v4(int c1, int* y, int d, FisherStateV4& s) {
#ifdef PROFILE_V4
    auto _t0 = g_prof ? std::chrono::high_resolution_clock::now()
                       : std::chrono::high_resolution_clock::time_point();
    if (g_prof) g_prof->calls_find_min++;
#endif
    int key = cache_key(d, c1);

    auto it = s.memo_min.find(key);
    if (it != s.memo_min.end()) {
#ifdef PROFILE_V4
        if (g_prof) g_prof->cache_hits_min++;
#endif
        const CacheEntry& entry = it->second;
        std::memcpy(y + d, entry.y + d, (s.m - d) * sizeof(int));
        s.p_min = compute_prob_v4(y, s);
#ifdef PROFILE_V4
        if (g_prof) {
            auto _t1 = std::chrono::high_resolution_clock::now();
            g_prof->time_find_min += std::chrono::duration<double>(
                _t1 - _t0).count();
        }
#endif
        return;
    }

#ifdef PROFILE_V4
    if (g_prof) g_prof->cache_misses_min++;
#endif

    int r_avail[MAX_ROWS];
    int n_rem = 0;
    for (int i = 0; i < s.m; i++) {
        if (i < d) r_avail[i] = -1;
        else { r_avail[i] = s.r[i]; n_rem += s.r[i]; }
    }

    double local_min_p = 1e300;
    int local_min_y[MAX_ROWS];
    std::memcpy(local_min_y, y, s.m * sizeof(int));

    joe_min_v4_impl(d, true, r_avail, c1, n_rem, y, s, local_min_p, local_min_y);

    s.p_min = local_min_p;

    CacheEntry entry;
    std::memcpy(entry.y + d, local_min_y + d, (s.m - d) * sizeof(int));
    s.memo_min[key] = entry;

#ifdef PROFILE_V4
    if (g_prof) {
        auto _t1 = std::chrono::high_resolution_clock::now();
        g_prof->time_find_min += std::chrono::duration<double>(
            _t1 - _t0).count();
    }
#endif
}

struct FindMaxEntry {
    int k, c1, n_rem;
    int y[MAX_ROWS];
};

void find_max_v4(int k_start, int c1, int* y, FisherStateV4& s) {
#ifdef PROFILE_V4
    auto _t0 = g_prof ? std::chrono::high_resolution_clock::now()
                       : std::chrono::high_resolution_clock::time_point();
    if (g_prof) g_prof->calls_find_max++;
#endif
    int key = cache_key(k_start, c1);

    auto it = s.memo_max.find(key);
    if (it != s.memo_max.end()) {
#ifdef PROFILE_V4
        if (g_prof) g_prof->cache_hits_max++;
#endif
        const CacheEntry& entry = it->second;
        std::memcpy(y + k_start, entry.y + k_start,
                    (s.m - k_start) * sizeof(int));
        s.p_max = compute_prob_v4(y, s);
        std::memcpy(s.y_max, y, s.m * sizeof(int));
#ifdef PROFILE_V4
        if (g_prof) {
            auto _t1 = std::chrono::high_resolution_clock::now();
            g_prof->time_find_max += std::chrono::duration<double>(
                _t1 - _t0).count();
        }
#endif
        return;
    }

#ifdef PROFILE_V4
    if (g_prof) g_prof->cache_misses_max++;
#endif

    int n_rem = s.suffix_r[k_start];
    double local_max_p = 0;
    int local_max_y[MAX_ROWS];
    std::memcpy(local_max_y, y, s.m * sizeof(int));

    FindMaxEntry stack[256];
    int sp = 0;

    stack[sp].k = k_start;
    stack[sp].c1 = c1;
    stack[sp].n_rem = n_rem;
    std::memcpy(stack[sp].y, y, s.m * sizeof(int));
    sp++;

    while (sp > 0) {
        sp--;
        FindMaxEntry curr = stack[sp];
        int k = curr.k;
        int c1_rem = curr.c1;
        int n_rem_curr = curr.n_rem;

        if (k >= s.m) {
            double prob = compute_prob_v4(curr.y, s);
            if (prob > local_max_p) {
                local_max_p = prob;
                std::memcpy(local_max_y, curr.y, s.m * sizeof(int));
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

        int n_after = s.suffix_r[k + 1];
        y_lo = std::max(y_lo, std::max(0, c1_rem - n_after));
        y_up = std::min(y_up, std::min(s.r[k], c1_rem));

        if (y_lo > y_up) continue;

        for (int y_k = y_lo; y_k <= y_up && sp < 255; y_k++) {
            stack[sp].k = k + 1;
            stack[sp].c1 = c1_rem - y_k;
            stack[sp].n_rem = n_rem_curr - s.r[k];
            std::memcpy(stack[sp].y, curr.y, s.m * sizeof(int));
            stack[sp].y[k] = y_k;
            sp++;
        }
    }

    s.p_max = local_max_p;
    std::memcpy(s.y_max, local_max_y, s.m * sizeof(int));

    CacheEntry entry;
    std::memcpy(entry.y + k_start, local_max_y + k_start,
                (s.m - k_start) * sizeof(int));
    s.memo_max[key] = entry;

#ifdef PROFILE_V4
    if (g_prof) {
        auto _t1 = std::chrono::high_resolution_clock::now();
        g_prof->time_find_max += std::chrono::duration<double>(
            _t1 - _t0).count();
    }
#endif
}

inline double dhyper_v4(int x, int m, int n, int k, const double* lfact) {
    if (x < 0 || x > k || x > m || k - x > n) return 0.0;
    double lp = lfact[m] - lfact[x] - lfact[m - x]
              + lfact[n] - lfact[k - x] - lfact[n - k + x]
              - lfact[m + n] + lfact[k] + lfact[m + n - k];
    return std::exp(lp);
}

void traverse_v4(int k, int c1, int* y, double prob_prefix, int n_rem,
                 bool descend, int y_k, int dir, int mode, FisherStateV4& s) {
#ifdef PROFILE_V4
    if (g_prof) g_prof->calls_traverse++;
#endif

    if (k == s.m - 1) {
        y[s.m - 1] = c1;
        double hp = dhyper_v4(c1, c1, n_rem - c1, s.r[s.m - 1], s.lfact);
        double prob = prob_prefix * hp;
        if (prob > s.p_obs * (1 + s.tol)) s.pval -= prob;
        return;
    }

    int y_lo = std::max(0, c1 - s.suffix_r[k + 1]);
    int y_hi = std::min(s.r[k], c1);

    if (descend) { mode = s.y_max[k]; y_k = mode; }
    if (y_k < y_lo) return;
    if (y_k > y_hi) {
        traverse_v4(k, c1, y, prob_prefix, n_rem, false, mode - 1, -1, mode, s);
        return;
    }

    y[k] = y_k;
    double hp = dhyper_v4(y_k, c1, n_rem - c1, s.r[k], s.lfact);
    double new_prefix = prob_prefix * hp;
    int new_c1 = c1 - y_k;
    int new_n_rem = s.suffix_r[k + 1];

    if (s.m - k > 2) {
        if (new_prefix * s.suffix_max[k + 1] <= s.p_obs * (1 + s.tol)) {
#ifdef PROFILE_V4
            if (g_prof) g_prof->nodes_pruned_suffix++;
#endif
            traverse_v4(k, c1, y, prob_prefix, n_rem, false, y_k + dir, dir, mode, s);
            return;
        }

        s.p_max = 0;
        find_max_v4(k + 1, new_c1, y, s);
        s.p_min = 1;
        find_min_v4(new_c1, y, k + 1, s);

        if (s.p_min > s.p_obs * (1 + s.tol)) {
            s.pval -= new_prefix;
#ifdef PROFILE_V4
            if (g_prof) g_prof->nodes_pruned_min++;
#endif
            traverse_v4(k, c1, y, prob_prefix, n_rem, false, y_k + dir, dir, mode, s);
            return;
        }
        if (s.p_max <= s.p_obs * (1 + s.tol)) {
#ifdef PROFILE_V4
            if (g_prof) g_prof->nodes_pruned_max++;
#endif
            traverse_v4(k, c1, y, prob_prefix, n_rem, false, y_k + dir, dir, mode, s);
            return;
        }
    }

    traverse_v4(k + 1, new_c1, y, new_prefix, new_n_rem, true, 0, +1, 0, s);
    traverse_v4(k, c1, y, prob_prefix, n_rem, false, y_k + dir, dir, mode, s);
}

// [[Rcpp::export(name = ".tree_memo_cpp")]]
double tree_memo_cpp(IntegerMatrix dat) {
    FisherStateV4 s;
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

    bool flipped = s.c[0] > s.c[1];
    if (flipped) std::swap(s.c[0], s.c[1]);

    int order[MAX_ROWS];
    for (int i = 0; i < s.m; i++) order[i] = i;
    std::sort(order, order + s.m, [&](int a, int b) {
        return r_unsorted[a] < r_unsorted[b];
    });

    int y_obs[MAX_ROWS];
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

    s.lfact[0] = 0.0;
    for (int i = 1; i <= s.n; i++)
        s.lfact[i] = s.lfact[i - 1] + std::log((double)i);

    s.log_const = 0.0;
    for (int i = 0; i < s.m; i++)
        s.log_const += s.lfact[s.r[i]];
    s.log_const += s.lfact[s.c[0]] + s.lfact[s.c[1]] - s.lfact[s.n];

    double log_p_obs = s.log_const;
    for (int i = 0; i < s.m; i++) {
        log_p_obs -= s.lfact[y_obs[i]];
        log_p_obs -= s.lfact[s.r[i] - y_obs[i]];
    }
    s.p_obs = std::exp(log_p_obs);

    for (int i = 0; i < s.m; i++) s.y_max[i] = 0;
    s.p_max = 0;
    s.p_min = 1;

    find_max_v4(0, s.c[0], y_obs, s);
    traverse_v4(0, s.c[0], y_obs, 1.0, s.n, true, 0, +1, 0, s);

    return s.pval;
}

// [[Rcpp::export(name = ".tree_memo_profile")]]
Rcpp::List tree_memo_profile(IntegerMatrix dat) {
#ifndef PROFILE_V4
    Rcpp::stop("PROFILE_V4 not enabled. Uncomment in src/Makevars and rebuild.");
    return Rcpp::List();
#else
    ProfileCounters prof;
    g_prof = &prof;

    auto wall_start = std::chrono::high_resolution_clock::now();

    FisherStateV4 s;
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

    bool flipped = s.c[0] > s.c[1];
    if (flipped) std::swap(s.c[0], s.c[1]);

    int order[MAX_ROWS];
    for (int i = 0; i < s.m; i++) order[i] = i;
    std::sort(order, order + s.m, [&](int a, int b) {
        return r_unsorted[a] < r_unsorted[b];
    });

    int y_obs[MAX_ROWS];
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

    s.lfact[0] = 0.0;
    for (int i = 1; i <= s.n; i++)
        s.lfact[i] = s.lfact[i - 1] + std::log((double)i);

    s.log_const = 0.0;
    for (int i = 0; i < s.m; i++)
        s.log_const += s.lfact[s.r[i]];
    s.log_const += s.lfact[s.c[0]] + s.lfact[s.c[1]] - s.lfact[s.n];

    double log_p_obs = s.log_const;
    for (int i = 0; i < s.m; i++) {
        log_p_obs -= s.lfact[y_obs[i]];
        log_p_obs -= s.lfact[s.r[i] - y_obs[i]];
    }
    s.p_obs = std::exp(log_p_obs);

    for (int i = 0; i < s.m; i++) s.y_max[i] = 0;
    s.p_max = 0;
    s.p_min = 1;

    auto traverse_start = std::chrono::high_resolution_clock::now();

    find_max_v4(0, s.c[0], y_obs, s);
    traverse_v4(0, s.c[0], y_obs, 1.0, s.n, true, 0, +1, 0, s);

    auto traverse_end = std::chrono::high_resolution_clock::now();
    prof.time_traverse = std::chrono::duration<double>(
        traverse_end - traverse_start).count();

    auto wall_end = std::chrono::high_resolution_clock::now();
    double wall_total = std::chrono::duration<double>(
        wall_end - wall_start).count();

    g_prof = nullptr;

    return Rcpp::List::create(
        Rcpp::Named("pvalue") = s.pval,
        Rcpp::Named("wall_total") = wall_total,
        Rcpp::Named("time_traverse") = prof.time_traverse,
        Rcpp::Named("time_find_max") = prof.time_find_max,
        Rcpp::Named("time_find_min") = prof.time_find_min,
        Rcpp::Named("calls_traverse") = (double)prof.calls_traverse,
        Rcpp::Named("calls_find_max") = (double)prof.calls_find_max,
        Rcpp::Named("calls_find_min") = (double)prof.calls_find_min,
        Rcpp::Named("calls_compute_prob") = (double)prof.calls_compute_prob,
        Rcpp::Named("cache_hits_max") = (double)prof.cache_hits_max,
        Rcpp::Named("cache_hits_min") = (double)prof.cache_hits_min,
        Rcpp::Named("cache_misses_max") = (double)prof.cache_misses_max,
        Rcpp::Named("cache_misses_min") = (double)prof.cache_misses_min,
        Rcpp::Named("nodes_pruned_suffix") = (double)prof.nodes_pruned_suffix,
        Rcpp::Named("nodes_pruned_max") = (double)prof.nodes_pruned_max,
        Rcpp::Named("nodes_pruned_min") = (double)prof.nodes_pruned_min,
        Rcpp::Named("memo_max_size") = (double)s.memo_max.size(),
        Rcpp::Named("memo_min_size") = (double)s.memo_min.size()
    );
#endif
}
