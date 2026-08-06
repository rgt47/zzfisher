// rxc_tree_memo.cpp
// First-generation m x c Fisher's exact test kernel.
//
// Notation follows the Concrete Mathematics convention of the
// manuscripts: m rows, c columns, n grand total, row margins r_i,
// column margins c_j, cell coordinates (i, j), path index t.
//
// Algorithm (from Rxc1 paper pp.21-28, network/path reduction):
//   - Linearize the (m-1)(c-1) free cells along the path
//     (1,1) -> (2,1) -> ... -> (m-1,1) -> (1,2) -> ... -> (m-1,c-1)
//     using t = i + (j-1)(m-1) (1-indexed in the paper; 0-indexed
//     here).
//   - Recursively place each free cell over its feasible range
//     [max(0, rresid - sum_{j'>j} cresid[j']), min(rresid, cresid)].
//   - When the path completes, fill the last column (cells (i,c-1) =
//     remaining row residual) and last row (cells (m-1,j) = remaining
//     column residual). This determines the table.
//   - Test statistic T(Y) = 1/Pr(Y); start antipvalue at 1 and subtract
//     Pr(Y) for tables strictly more probable than observed
//     (Pr(Y) > Pr_obs * (1 + tol)). Result is the two-sided p-value.
//
// Tolerance constant matches the rx2 family: 3.45254e-7 (R fisher.test).
//
// Pruning: acceptance-region skip with a capped water-filling
// bound (the Joe-style exact-relaxation analog of rx2_tree_memo's
// suffix extrema, valid for any c) plus a concavity-based sibling
// cutoff (the r x c analog of the m x 2 cascade). See the bound
// comment below. The complement-side bulk subtraction of the
// m x 2 kernels has no r x c analogue (no closed-form suffix
// mass for c >= 3).

#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <cmath>

using namespace Rcpp;

namespace {

const double RXC_TOL = 3.45254e-7;

struct RxcState {
  int m, c, n;                       // rows, columns, grand total
  std::vector<int> rmarg, cmarg;     // observed row/column margins
  std::vector<int> rresid, cresid;   // residual row/column margins
  std::vector<int> y;                // current table, m*c ints, row-major
  double log_const;                  // log(prod r_i! prod c_j! / n!)
  double p_obs;
  double log_thresh;                 // log(p_obs * (1 + tol))
  double pval;                       // antipvalue, starts at 1
  std::vector<double> lf;            // lf[v] = log(v!)
  // Per-depth caps scratch. One buffer per path position: a single
  // shared buffer would be clobbered by the recursive calls made
  // between two siblings, leaving the parent's later children
  // bounded with a descendant's caps (observed defect: p-values
  // wrong by up to 0.8 from wrongly skipped live subtrees).
  std::vector<std::vector<int>> caps_at;
  long n_skipped;                    // subtrees skipped by the bound
};

inline int yidx(const RxcState& s, int i, int j) { return i * s.c + j; }


// --- Acceptance-region skip via capped water-filling -------------
//
// The m x 2 kernel prunes with exact subtree extrema, which exist
// there because two columns make log Pr separable into concave
// one-dimensional terms. That separability is unavailable for
// c >= 3: the cells of a partial table are tied by both row and
// column constraints simultaneously.
//
// What transfers is the relaxation, and for the current column it
// can respect the row residuals. Maximize
//
//   sum_t -log(y_t!)  s.t.  sum_t y_t = budget, 0 <= y_t <= cap_t,
//
// where the caps are the residuals of the rows still open in the
// current column. The objective is separable and concave with
// identical terms, so the maximum is capped water-filling:
// allocate as evenly as possible, saturating small caps (the
// greedy-exactness lemma of the c = 2 white paper, Section 2,
// applied with box constraints). Dropping the coupling of a row's
// residual across columns keeps this a valid relaxation, and it is
// never looser than the uncapped balanced split (caps of +infinity
// recover it). Later columns use the uncapped split: their budgets
// do not depend on the value placed at the current cell, so that
// part of the bound is a per-node constant, and a fully capped
// per-child bound measured 2.7x slower end to end than this
// split — the extra tightness did not pay for its cost. When the
// bound on every completion's log-probability does not exceed
// log(Pr_obs * (1 + tol)), the subtree lies in the significance
// region S and is skipped exactly.
//
// Convexity gives a sibling cutoff on top. The child bound as a
// function of the value x placed at the current cell is
// -log(x!) plus value functions of concave programs whose
// right-hand sides (the column budget and row-i cap) move linearly
// in x, hence concave in x. The children that survive the bound
// therefore form a contiguous interval, and the scan over x stops
// at the first sub-threshold child after the interval — the r x c
// analog of the m x 2 cascade.
//
// The complement side still has no analogue: bulk subtraction
// needs a subtree's total mass, which for c = 2 is Vandermonde's
// convolution but for c >= 3 has no closed form.
// See docs/c2_separability_whitepaper_2026-08-04.md, Section 5.

// Max of sum -log(y_t!) over sum y_t = budget, 0 <= y_t <= caps[t],
// caps pre-sorted ascending. -infinity if the caps cannot absorb
// the budget.
double water_fill_lp(const RxcState& s, int budget,
                     const std::vector<int>& caps, int q) {
  if (q <= 0) return (budget == 0) ? 0.0 : -INFINITY;
  double bound = 0.0;
  int rem_budget = budget;
  for (int t = 0; t < q; ++t) {
    const int cells_left = q - t;
    const int level = rem_budget / cells_left;
    if (caps[t] <= level) {
      // Small cap saturates; redistribute the rest more thickly.
      bound -= s.lf[caps[t]];
      rem_budget -= caps[t];
    } else {
      // No remaining cap binds: balanced split of the remainder.
      const int r2 = rem_budget % cells_left;
      bound -= (double)(cells_left - r2) * s.lf[level]
             + (double)r2 * s.lf[level + 1];
      return bound;
    }
  }
  return (rem_budget == 0) ? bound : -INFINITY;
}

// Uncapped balanced split: water-filling with no binding caps.
// Used for later columns, whose budgets do not depend on the value
// placed at the current cell, so this part of the bound is a
// per-node constant.
inline double balanced_split_lp(const RxcState& s, int q, int total) {
  if (q <= 0) return (total == 0) ? 0.0 : -INFINITY;
  const int base = total / q, rem = total % q;
  return -((double)(q - rem) * s.lf[base] +
           (double)rem * s.lf[base + 1]);
}

// After path traversal completes, fill last column and last row from
// margin residuals. Returns false if any computed cell is negative.
//
// Implementation note: this function does NOT mutate s.cresid or
// s.rresid (only s.y). An earlier version decremented s.cresid[c-1]
// in place; on the infeasibility-return code path that mutation was
// not reverted and corrupted subsequent iterations. The current form
// computes the last (m-1, c-1) cell from a local sum and reads the
// last-row cells directly from s.cresid[j], which the caller has not
// yet mutated.
bool finish_table(RxcState& s) {
  const int mm1 = s.m - 1;
  const int cm1 = s.c - 1;
  // Last column for rows 0..m-2: cell value = remaining row residual.
  int sum_last_col = 0;
  for (int i = 0; i < mm1; ++i) {
    int v = s.rresid[i];
    if (v < 0) return false;
    s.y[yidx(s, i, cm1)] = v;
    sum_last_col += v;
  }
  // (m-1, c-1) cell from last-column total minus what was just placed.
  int v_corner = s.cmarg[cm1] - sum_last_col;
  if (v_corner < 0) return false;
  s.y[yidx(s, mm1, cm1)] = v_corner;
  // Last row for cols 0..c-2: cell value = remaining column residual
  // (s.cresid[j] is the residual after path traversal in cols < c-1).
  for (int j = 0; j < cm1; ++j) {
    int v = s.cresid[j];
    if (v < 0) return false;
    s.y[yidx(s, mm1, j)] = v;
  }
  return true;
}

// Children are bounded before recursing, so a subtree is skipped
// without entering it, and by concavity of the child bound in x
// the surviving children form an interval: the scan breaks at the
// first sub-threshold child after a surviving one.
void traverse(int t, int tmax, RxcState& s, double placed_lp) {
  if (t > tmax) {
    if (!finish_table(s)) return;
    // placed_lp carries -log(x!) of every free cell, so only the
    // determined last column and last row need summing.
    double det = 0.0;
    const int mm1 = s.m - 1, cm1 = s.c - 1;
    for (int ii = 0; ii < s.m; ++ii)
      det += s.lf[s.y[yidx(s, ii, cm1)]];
    for (int jj = 0; jj < cm1; ++jj)
      det += s.lf[s.y[yidx(s, mm1, jj)]];
    const double p = std::exp(s.log_const + placed_lp - det);
    if (p - s.p_obs > s.p_obs * RXC_TOL) s.pval -= p;
    return;
  }
  // Path step t (1-indexed) -> cell (i, j) zero-indexed.
  const int tt = t - 1;
  const int i = tt % (s.m - 1);
  const int j = tt / (s.m - 1);

  // Lower bound: row i residual must be absorbable by remaining-cols
  // capacity (cols j+1..c-1), so x >= rresid - sum of those cresid.
  int sum_rest_cols = 0;
  for (int jj = j + 1; jj < s.c; ++jj) sum_rest_cols += s.cresid[jj];
  const int x_lo = std::max(0, s.rresid[i] - sum_rest_cols);
  const int x_up = std::min(s.rresid[i], s.cresid[j]);

  // Per-node constants of the child bound: the later columns'
  // uncapped balanced splits (their budgets do not involve x) and
  // the current column's caps (rows i+1..m-1; row i is excluded,
  // so the caps do not involve x either), sorted once.
  double later = 0.0;
  for (int jj = j + 1; jj < s.c; ++jj)
    later += balanced_split_lp(s, s.m, s.cresid[jj]);
  std::vector<int>& caps = s.caps_at[t];
  caps.clear();
  for (int ii = i + 1; ii < s.m; ++ii)
    caps.push_back(s.rresid[ii]);
  std::sort(caps.begin(), caps.end());
  const int q_cur = (int)caps.size();
  const double node_const = s.log_const + placed_lp + later;

  bool seen_live = false;
  for (int x = x_lo; x <= x_up; ++x) {
    const double bnd = node_const - s.lf[x] +
      water_fill_lp(s, s.cresid[j] - x, caps, q_cur);
    if (bnd <= s.log_thresh) {
      ++s.n_skipped;
      if (seen_live) break;   // concavity: past the live interval
      continue;
    }
    seen_live = true;
    s.y[yidx(s, i, j)] = x;
    s.rresid[i] -= x;
    s.cresid[j] -= x;
    traverse(t + 1, tmax, s, placed_lp - s.lf[x]);
    s.rresid[i] += x;
    s.cresid[j] += x;
  }
  s.y[yidx(s, i, j)] = 0;
}

}  // anonymous namespace

// [[Rcpp::export(name = ".rxc_tree_memo_cpp")]]
double rxc_tree_memo_cpp(IntegerMatrix dat) {
  RxcState s;
  s.m = dat.nrow();
  s.c = dat.ncol();
  if (s.m < 2) Rcpp::stop("dat must have at least 2 rows");
  if (s.c < 2) Rcpp::stop("dat must have at least 2 columns");

  s.rmarg.assign(s.m, 0);
  s.cmarg.assign(s.c, 0);
  s.y.assign(s.m * s.c, 0);
  for (int i = 0; i < s.m; ++i) {
    for (int j = 0; j < s.c; ++j) {
      const int v = dat(i, j);
      if (v < 0) Rcpp::stop("dat entries must be non-negative");
      s.rmarg[i] += v;
      s.cmarg[j] += v;
    }
  }
  s.n = 0;
  for (int v : s.rmarg) s.n += v;

  double sum_log_r = 0.0, sum_log_c = 0.0;
  for (int v : s.rmarg) sum_log_r += std::lgamma(v + 1);
  for (int v : s.cmarg) sum_log_c += std::lgamma(v + 1);
  s.log_const = sum_log_r + sum_log_c - std::lgamma(s.n + 1);

  double sum_log_obs = 0.0;
  for (int i = 0; i < s.m; ++i)
    for (int j = 0; j < s.c; ++j)
      sum_log_obs += std::lgamma(dat(i, j) + 1);
  s.p_obs = std::exp(s.log_const - sum_log_obs);

  s.rresid = s.rmarg;
  s.cresid = s.cmarg;
  s.pval = 1.0;
  s.n_skipped = 0;
  s.log_thresh = std::log(s.p_obs) + std::log1p(RXC_TOL);

  // Log-factorial table, sized to the grand total (no cell exceeds
  // it), so the balanced-split bound is a table lookup.
  s.lf.assign(s.n + 2, 0.0);
  for (int v = 1; v <= s.n + 1; ++v)
    s.lf[v] = s.lf[v - 1] + std::log((double)v);

  const int tmax = (s.m - 1) * (s.c - 1);
  s.caps_at.assign(tmax + 1, std::vector<int>());
  for (auto& v : s.caps_at) v.reserve(s.m);
  traverse(1, tmax, s, 0.0);

  // Numeric guard: pval should be in (0, 1]. Clamp to handle round-off
  // around the boundary cases.
  if (s.pval < 0.0) s.pval = 0.0;
  if (s.pval > 1.0) s.pval = 1.0;
  return s.pval;
}
