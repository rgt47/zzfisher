// rx2_net_sweep.cpp
// State-sweep kernel for Fisher's exact test on m x 2 tables:
// Mehta-Patel-Joe-Clarkson-Requena specialized to two columns, with
// the hashing replaced by dense arrays. Implements the design of
// docs/c2_separability_whitepaper_2026-08-04.md, Sections 3B and
// 3B.1:
//
//   - States are (stage k, spent budget a); the sweep proceeds
//     stage by stage, never revisiting a state (Mehta-Patel 1983).
//   - At each state, the paths reaching it are aggregated into
//     classes by prefix log-mass q, merged within a tight epsilon
//     (the PAST device of Clarkson-Fan-Joe 1993); classes sorted
//     by q classify into three contiguous runs by the monotone
//     tests of Proposition 6D, using the exact suffix extrema
//     maxh/minh of the c = 2 separability DP.
//   - A class entirely in the complement region U contributes
//     (class weight) x (total suffix mass) in O(1), the suffix
//     mass being closed-form by Vandermonde. A class entirely in
//     the significance region S is dropped. Only straddling
//     classes propagate, per state, not per path.
//   - At the last two rows, Proposition 6I: leaf masses and their
//     prefix sums are computed once per state; each class's
//     qualifying completions form an interval, nested across
//     classes sorted by q, extracted by expanding two pointers;
//     each class's whole contribution is one prefix-sum
//     difference.
//
// The class-merge epsilon is 1e-12 on the log scale: it merges
// paths whose prefix products are equal up to floating-point noise
// (the source of the measured 1.6-13x multiplicity, chiefly tied
// margins) without merging genuinely distinct masses, so exactness
// against a complete-enumeration oracle is preserved. FEXACT's own
// tolerance enters only through the significance threshold
// log_thresh = log p_obs + log1p(tol), exactly as in the tree
// kernel.
//
// No fixed workspace: storage grows with the live class count.

#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <cmath>

using namespace Rcpp;

namespace {

struct ClassEntry {
  double q;   // prefix log-mass: sum_{i < k} h_i(y_i)
  double w;   // linear multiplicity relative to q (paths merged in)
};

struct SweepState {
  int m, n;
  std::vector<int> r;          // row margins, sorted ascending
  int c1;                      // smaller column margin (budget)
  double log_const;
  double p_obs;
  double log_thresh;
  std::vector<double> lfact;   // lfact[v] = log v!
  std::vector<double> maxh;    // (m+1) x (c1+1) exact suffix maxima
  std::vector<double> minh;
  std::vector<int> suffix_r;   // suffix_r[k] = sum_{i >= k} r[i]
  std::vector<double> slf;     // slf[k] = sum_{i >= k} lfact[r[i]]
  int acap;                    // c1 + 1
  double acc;                  // accumulated mass of U (linear)
};

inline double h_term(const SweepState& s, int k, int y) {
  return -s.lfact[y] - s.lfact[s.r[k] - y];
}

inline double lchoose_lf(const SweepState& s, int nn, int kk) {
  return s.lfact[nn] - s.lfact[kk] - s.lfact[nn - kk];
}

// Exact suffix extrema of sum h_i under a budget, by backward DP
// (same construction as the redesigned tree kernel).
void build_extrema(SweepState& s) {
  int m = s.m, acap = s.acap;
  s.maxh.assign((size_t)(m + 1) * acap, -INFINITY);
  s.minh.assign((size_t)(m + 1) * acap, INFINITY);
  s.maxh[(size_t)m * acap] = 0.0;
  s.minh[(size_t)m * acap] = 0.0;
  for (int k = m - 1; k >= 0; k--) {
    int amax = std::min(acap - 1, s.suffix_r[k]);
    for (int a = 0; a <= amax; a++) {
      double bmax = -INFINITY, bmin = INFINITY;
      int y_lo = std::max(0, a - s.suffix_r[k + 1]);
      int y_hi = std::min(s.r[k], a);
      for (int y = y_lo; y <= y_hi; y++) {
        double nx = s.maxh[(size_t)(k + 1) * acap + (a - y)];
        double nn = s.minh[(size_t)(k + 1) * acap + (a - y)];
        double h = h_term(s, k, y);
        if (std::isfinite(nx) && h + nx > bmax) bmax = h + nx;
        if (std::isfinite(nn) && h + nn < bmin) bmin = h + nn;
      }
      s.maxh[(size_t)k * acap + a] = bmax;
      s.minh[(size_t)k * acap + a] = bmin;
    }
  }
}

// Sort a state's classes by q and merge within eps; multiplicities
// combine exactly: w_rep += w * exp(q - q_rep).
void sort_merge(std::vector<ClassEntry>& cl, double eps) {
  std::sort(cl.begin(), cl.end(),
            [](const ClassEntry& x, const ClassEntry& y) {
              return x.q < y.q;
            });
  size_t out = 0;
  for (size_t i = 0; i < cl.size(); i++) {
    if (out > 0 && cl[i].q - cl[out - 1].q <= eps) {
      cl[out - 1].w += cl[i].w * std::exp(cl[i].q - cl[out - 1].q);
    } else {
      cl[out++] = cl[i];
    }
  }
  cl.resize(out);
}

void run_sweep(SweepState& s, double eps_merge) {
  const int c1 = s.c1;
  std::vector<std::vector<ClassEntry>> cur(c1 + 1), nxt(c1 + 1);
  cur[0].push_back({0.0, 1.0});

  for (int k = 0; k < s.m; k++) {
    const int rem = s.m - k;
    bool any_next = false;

    for (int a = 0; a <= c1; a++) {
      std::vector<ClassEntry>& cl = cur[a];
      if (cl.empty()) continue;
      const int b = c1 - a;
      if (b > s.suffix_r[k]) { cl.clear(); continue; }
      const double mx = s.maxh[(size_t)k * s.acap + b];
      if (!std::isfinite(mx)) { cl.clear(); continue; }
      const double mn = s.minh[(size_t)k * s.acap + b];

      sort_merge(cl, eps_merge);

      // Proposition 6I state-level leaf structures (rem == 2),
      // built lazily on the first straddling class.
      bool leaf_ready = false;
      int y_lo_st = 0, y_hi_st = -1, ym = 0;
      std::vector<double> flog, F;
      int cur_lo = 0, cur_hi = -1;   // empty interval

      for (const ClassEntry& e : cl) {
        const double base = s.log_const + e.q;

        // 6D run 1: entirely in the significance region -> drop.
        if (base + mx <= s.log_thresh) continue;

        // 6D run 3: entirely in the complement region -> bulk.
        if (base + mn > s.log_thresh) {
          s.acc += e.w * std::exp(base
                     + lchoose_lf(s, s.suffix_r[k], b) - s.slf[k]);
          continue;
        }

        // Straddling.
        if (rem >= 3) {
          const int y_lo = std::max(0, b - s.suffix_r[k + 1]);
          const int y_hi = std::min(s.r[k], b);
          for (int y = y_lo; y <= y_hi; y++) {
            nxt[a + y].push_back({e.q + h_term(s, k, y), e.w});
          }
          any_next = true;
          continue;
        }

        // rem == 2: Proposition 6I. (rem == 1 cannot straddle:
        // maxh == minh there, so one of the runs above caught it.)
        if (!leaf_ready) {
          y_lo_st = std::max(0, b - s.r[k + 1]);
          y_hi_st = std::min(s.r[k], b);
          const int w = y_hi_st - y_lo_st + 1;
          flog.resize(w);
          F.assign(w + 1, 0.0);
          for (int y = y_lo_st; y <= y_hi_st; y++) {
            const double fl = h_term(s, k, y)
                            + h_term(s, k + 1, b - y);
            flog[y - y_lo_st] = fl;
            F[y - y_lo_st + 1] = F[y - y_lo_st] + std::exp(fl);
          }
          ym = (int)(std::max_element(flog.begin(), flog.end())
                     - flog.begin());
          cur_lo = ym + 1; cur_hi = ym;      // empty
          leaf_ready = true;
        }
        // Classes are sorted ascending in q, so the per-class
        // threshold c decreases and the qualifying interval only
        // grows: expand the two pointers, never contract.
        const double c = s.log_thresh - base;
        if (cur_hi < cur_lo) {
          if (flog[ym] > c) { cur_lo = ym; cur_hi = ym; }
        }
        if (cur_hi >= cur_lo) {
          while (cur_lo > 0 && flog[cur_lo - 1] > c) cur_lo--;
          while (cur_hi < (int)flog.size() - 1 &&
                 flog[cur_hi + 1] > c) cur_hi++;
          s.acc += e.w * std::exp(base)
                 * (F[cur_hi + 1] - F[cur_lo]);
        }
      }
      cl.clear();
    }

    if (!any_next) break;
    std::swap(cur, nxt);
  }
}

}  // anonymous namespace

// eps_merge: class-merge width on the log scale. The default
// 1e-12 merges only floating-point ties and preserves exactness.
// FEXACT's own grouping width is its tolerance 3.45254e-7; passing
// that reproduces FEXACT's aggregation semantics at the cost of a
// bounded relative error of the same order in the merged masses.
// [[Rcpp::export(name = ".rx2_net_sweep_cpp")]]
double rx2_net_sweep_cpp(IntegerMatrix dat,
                         double eps_merge = 1e-12) {
  SweepState s;
  s.m = dat.nrow();
  if (s.m < 2) return 1.0;

  std::vector<int> r_unsorted(s.m);
  int c0 = 0, c1col = 0;
  s.n = 0;
  for (int i = 0; i < s.m; i++) {
    r_unsorted[i] = dat(i, 0) + dat(i, 1);
    s.n += r_unsorted[i];
    c0 += dat(i, 0);
    c1col += dat(i, 1);
  }
  const bool flipped = c0 > c1col;
  s.c1 = flipped ? c1col : c0;
  const int c2 = s.n - s.c1;
  s.acap = s.c1 + 1;

  std::vector<int> order(s.m);
  for (int i = 0; i < s.m; i++) order[i] = i;
  std::sort(order.begin(), order.end(), [&](int x, int y) {
    return r_unsorted[x] < r_unsorted[y];
  });
  s.r.resize(s.m);
  std::vector<int> y_obs(s.m);
  for (int i = 0; i < s.m; i++) {
    s.r[i] = r_unsorted[order[i]];
    y_obs[i] = flipped ? dat(order[i], 1) : dat(order[i], 0);
  }

  s.suffix_r.assign(s.m + 1, 0);
  for (int i = s.m - 1; i >= 0; i--)
    s.suffix_r[i] = s.suffix_r[i + 1] + s.r[i];

  s.lfact.assign(s.n + 1, 0.0);
  for (int v = 1; v <= s.n; v++)
    s.lfact[v] = s.lfact[v - 1] + std::log((double)v);

  s.slf.assign(s.m + 1, 0.0);
  for (int i = s.m - 1; i >= 0; i--)
    s.slf[i] = s.slf[i + 1] + s.lfact[s.r[i]];

  s.log_const = s.slf[0] + s.lfact[s.c1] + s.lfact[c2]
              - s.lfact[s.n];

  double log_p_obs = s.log_const;
  for (int i = 0; i < s.m; i++)
    log_p_obs += -s.lfact[y_obs[i]] - s.lfact[s.r[i] - y_obs[i]];
  s.p_obs = std::exp(log_p_obs);
  s.log_thresh = log_p_obs + std::log1p(3.45254e-7);

  build_extrema(s);

  s.acc = 0.0;
  run_sweep(s, eps_merge);

  double pval = 1.0 - s.acc;
  if (pval < 0.0) pval = 0.0;
  if (pval > 1.0) pval = 1.0;
  return pval;
}
