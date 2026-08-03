// rxck_tree_memo.cpp
// First-generation m x c x k Fisher's exact test kernel implementing
// the no-three-way-interaction model (Bartlett, [AB][AC][BC]).
//
// Notation follows the Concrete Mathematics convention of the
// manuscripts: m rows, c columns, k layers, cell coordinates
// (i, j, l), path index t. The letter k is reserved for the layer
// count and never used as a coordinate. The three pairwise margins
// carry the Bartlett factor labels: ab (row x column), ac
// (row x layer), bc (column x layer).
//
// All three pairwise margins are conditioned on:
//   y_ij. = sum_l y[i,j,l]   (row x column, AB)
//   y_i.l = sum_j y[i,j,l]   (row x layer, AC)
//   y_.jl = sum_i y[i,j,l]   (column x layer, BC)
//
// Under the null of no three-way interaction, the conditional
// distribution given these pairwise margins is
//   P(Y | margins) proportional to 1 / prod_{i,j,l} y_ijl!
// so the test statistic T(Y) = sum lgamma(y_ijl + 1) orders tables
// from least extreme (small T) to most extreme (large T).
//
// Path linearization (Rxc1 paper p.23, with handwritten correction):
//   t = i + (j-1)(m-1) + (l-1)(m-1)(c-1)     (1-indexed)
//   total free cells = (m-1)(c-1)(k-1)
//
// After placing all free cells, the dependent cells in the last layer,
// last column, and last row are computed from the pairwise margins;
// each computed value is checked for non-negativity AND for consistency
// with the alternative margin equation (since the over-determined cells
// must agree under any feasible solution).
//
// p-value = sum exp(-T(Y)) over { Y : T(Y) >= T(Y_obs) }
//           divided by the same sum over the full reference set.
// Both sums are accumulated by enumeration; the normalization constant
// of the conditional distribution is not closed form for general 3-way
// tables, so this ratio approach avoids needing it explicitly.
//
// First-gen scope notes:
//   - Pruning: feasibility-only via per-pairwise-margin residuals.
//     No multivariate-hypergeometric mode pruning yet.
//   - Validated against a pure-R brute-force enumerator for the 2x2x2
//     case in inst/tinytest/test_htest.R (one free cell, easy oracle).
//   - Larger m*c*k becomes intractable quickly (factorial growth in
//     reference set); this kernel is intended for small tables only.

#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <algorithm>

using namespace Rcpp;

static const double RXCK_TOL = 3.45254e-7;

namespace rxck {

struct State {
  int m, c, k;                 // rows, columns, layers
  // Pairwise margins (computed from observed table):
  std::vector<int> ab;         // row x col, indexed [i*c + j]
  std::vector<int> ac;         // row x layer
  std::vector<int> bc;         // col x layer
  // Running residuals during path traversal:
  std::vector<int> ab_resid;
  std::vector<int> ac_resid;
  std::vector<int> bc_resid;
  // Current table; flat m*c*k, indexed [(i*c + j)*k + l].
  std::vector<int> y;

  double T_obs;        // sum lgamma(y_obs + 1)
  double mass_total;   // accumulated sum exp(-T)
  double mass_extreme; // accumulated sum exp(-T) over { Y : T >= T_obs }
};

inline int yidx(const State& s, int i, int j, int l) {
  return (i * s.c + j) * s.k + l;
}
inline int abidx(const State& s, int i, int j) { return i * s.c + j; }
inline int acidx(const State& s, int i, int l) { return i * s.k + l; }
inline int bcidx(const State& s, int j, int l) { return j * s.k + l; }

double table_T(const State& s) {
  double sum = 0.0;
  for (int v : s.y) sum += std::lgamma(v + 1);
  return sum;
}

// Compute dependent cells from pairwise margins. Returns false if any
// resulting cell is negative or if a consistency check (over-determined
// cell values agreeing across margins) fails.
bool finish_table(State& s) {
  const int mm1 = s.m - 1;
  const int cm1 = s.c - 1;
  const int km1 = s.k - 1;

  // 1. Last layer for free (i, j) cells: y[i, j, k-1] from AB margin.
  for (int i = 0; i < mm1; ++i) {
    for (int j = 0; j < cm1; ++j) {
      int sum = 0;
      for (int l = 0; l < km1; ++l) sum += s.y[yidx(s, i, j, l)];
      int v = s.ab[abidx(s, i, j)] - sum;
      if (v < 0) return false;
      s.y[yidx(s, i, j, km1)] = v;
    }
  }
  // 2. Last column for free (i, l) cells (l < k-1): from AC margin.
  for (int i = 0; i < mm1; ++i) {
    for (int l = 0; l < km1; ++l) {
      int sum = 0;
      for (int j = 0; j < cm1; ++j) sum += s.y[yidx(s, i, j, l)];
      int v = s.ac[acidx(s, i, l)] - sum;
      if (v < 0) return false;
      s.y[yidx(s, i, cm1, l)] = v;
    }
  }
  // 3. Last row for free (j, l) cells (j < c-1, l < k-1): from BC margin.
  for (int j = 0; j < cm1; ++j) {
    for (int l = 0; l < km1; ++l) {
      int sum = 0;
      for (int i = 0; i < mm1; ++i) sum += s.y[yidx(s, i, j, l)];
      int v = s.bc[bcidx(s, j, l)] - sum;
      if (v < 0) return false;
      s.y[yidx(s, mm1, j, l)] = v;
    }
  }
  // 4. y[i, c-1, k-1] for i < m-1: from AB margin (sum over layers
  //    in the last column for row i must equal ab[i, c-1]).
  //    Cross-check with AC margin (sum over columns in last layer
  //    for row i must equal ac[i, k-1]).
  for (int i = 0; i < mm1; ++i) {
    int sum_l = 0;
    for (int l = 0; l < km1; ++l) sum_l += s.y[yidx(s, i, cm1, l)];
    int v = s.ab[abidx(s, i, cm1)] - sum_l;
    if (v < 0) return false;
    s.y[yidx(s, i, cm1, km1)] = v;
    int ac_check = 0;
    for (int j = 0; j < s.c; ++j) ac_check += s.y[yidx(s, i, j, km1)];
    if (ac_check != s.ac[acidx(s, i, km1)]) return false;
  }
  // 5. y[m-1, j, k-1] for j < c-1: from BC margin.
  //    Cross-check with AB margin.
  for (int j = 0; j < cm1; ++j) {
    int sum_i = 0;
    for (int i = 0; i < mm1; ++i) sum_i += s.y[yidx(s, i, j, km1)];
    int v = s.bc[bcidx(s, j, km1)] - sum_i;
    if (v < 0) return false;
    s.y[yidx(s, mm1, j, km1)] = v;
    int ab_check = 0;
    for (int l = 0; l < s.k; ++l) ab_check += s.y[yidx(s, mm1, j, l)];
    if (ab_check != s.ab[abidx(s, mm1, j)]) return false;
  }
  // 6. y[m-1, c-1, l] for l < k-1: from BC margin.
  //    Cross-check with AC margin.
  for (int l = 0; l < km1; ++l) {
    int sum_i = 0;
    for (int i = 0; i < mm1; ++i) sum_i += s.y[yidx(s, i, cm1, l)];
    int v = s.bc[bcidx(s, cm1, l)] - sum_i;
    if (v < 0) return false;
    s.y[yidx(s, mm1, cm1, l)] = v;
    int ac_check = 0;
    for (int j = 0; j < s.c; ++j) ac_check += s.y[yidx(s, mm1, j, l)];
    if (ac_check != s.ac[acidx(s, mm1, l)]) return false;
  }
  // 7. Final cell y[m-1, c-1, k-1]: any margin works; use AB.
  //    Cross-check with AC and BC.
  {
    int sum_l = 0;
    for (int l = 0; l < km1; ++l) sum_l += s.y[yidx(s, mm1, cm1, l)];
    int v = s.ab[abidx(s, mm1, cm1)] - sum_l;
    if (v < 0) return false;
    s.y[yidx(s, mm1, cm1, km1)] = v;
    int ac_check = 0;
    for (int j = 0; j < s.c; ++j) ac_check += s.y[yidx(s, mm1, j, km1)];
    if (ac_check != s.ac[acidx(s, mm1, km1)]) return false;
    int bc_check = 0;
    for (int i = 0; i < s.m; ++i) bc_check += s.y[yidx(s, i, cm1, km1)];
    if (bc_check != s.bc[bcidx(s, cm1, km1)]) return false;
  }
  return true;
}

void traverse(int t, int tmax, State& s) {
  if (t > tmax) {
    if (!finish_table(s)) return;
    double T = table_T(s);
    double w = std::exp(-T);
    s.mass_total += w;
    if (T - s.T_obs >= -RXCK_TOL) s.mass_extreme += w;
    return;
  }
  // Decode t (1-indexed) -> (i, j, l) zero-indexed via
  //   t_0 = i + j*(m-1) + l*(m-1)*(c-1)
  const int t0 = t - 1;
  const int mm1 = s.m - 1;
  const int cm1 = s.c - 1;
  const int per_layer = mm1 * cm1;
  const int l = t0 / per_layer;
  const int rem = t0 % per_layer;
  const int j = rem / mm1;
  const int i = rem % mm1;

  // Upper bound = min of the three pairwise residuals at this cell.
  // (These were decremented by previous placements at cells that share
  // each pairwise margin.)
  const int x_up = std::min({
    s.ab_resid[abidx(s, i, j)],
    s.ac_resid[acidx(s, i, l)],
    s.bc_resid[bcidx(s, j, l)]
  });
  if (x_up < 0) return;

  for (int x = 0; x <= x_up; ++x) {
    s.y[yidx(s, i, j, l)] = x;
    s.ab_resid[abidx(s, i, j)] -= x;
    s.ac_resid[acidx(s, i, l)] -= x;
    s.bc_resid[bcidx(s, j, l)] -= x;
    traverse(t + 1, tmax, s);
    s.ab_resid[abidx(s, i, j)] += x;
    s.ac_resid[acidx(s, i, l)] += x;
    s.bc_resid[bcidx(s, j, l)] += x;
  }
  s.y[yidx(s, i, j, l)] = 0;
}

} // namespace rxck

// [[Rcpp::export(name = ".rxck_tree_memo_cpp")]]
double rxck_tree_memo_cpp(IntegerVector dat) {
  RObject dim_attr = dat.attr("dim");
  if (dim_attr.isNULL())
    Rcpp::stop("dat must be a 3D integer array (no dim attribute)");
  IntegerVector d(dim_attr);
  if (d.size() != 3)
    Rcpp::stop("dat must be a 3D integer array");

  rxck::State s;
  s.m = d[0];
  s.c = d[1];
  s.k = d[2];
  if (s.m < 2 || s.c < 2 || s.k < 2)
    Rcpp::stop("each dimension must be at least 2");

  const int total = s.m * s.c * s.k;
  s.y.assign(total, 0);
  s.ab.assign(s.m * s.c, 0);
  s.ac.assign(s.m * s.k, 0);
  s.bc.assign(s.c * s.k, 0);

  // Source layout is column-major (R's array storage):
  //   dat[i + m*j + m*c*l] -> cell (i, j, l).
  // Internal layout uses [(i*c + j)*k + l].
  for (int l = 0; l < s.k; ++l) {
    for (int j = 0; j < s.c; ++j) {
      for (int i = 0; i < s.m; ++i) {
        int v = dat[i + s.m * j + s.m * s.c * l];
        if (v < 0) Rcpp::stop("dat entries must be non-negative");
        s.y[rxck::yidx(s, i, j, l)] = v;
        s.ab[rxck::abidx(s, i, j)] += v;
        s.ac[rxck::acidx(s, i, l)] += v;
        s.bc[rxck::bcidx(s, j, l)] += v;
      }
    }
  }
  s.T_obs = rxck::table_T(s);

  // Reset y so the recursive traversal builds fresh tables.
  std::fill(s.y.begin(), s.y.end(), 0);
  s.ab_resid = s.ab;
  s.ac_resid = s.ac;
  s.bc_resid = s.bc;
  s.mass_total = 0.0;
  s.mass_extreme = 0.0;

  const int tmax = (s.m - 1) * (s.c - 1) * (s.k - 1);
  rxck::traverse(1, tmax, s);

  if (s.mass_total <= 0.0)
    Rcpp::stop("rxck enumeration produced no feasible tables; "
               "check that pairwise margins are mutually consistent");

  double pval = s.mass_extreme / s.mass_total;
  if (pval < 0.0) pval = 0.0;
  if (pval > 1.0) pval = 1.0;
  return pval;
}
