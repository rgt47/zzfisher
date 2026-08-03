// rxck_tree_memo.cpp
// First-generation r x c x k Fisher's exact test kernel implementing
// the no-three-way-interaction model (Bartlett, [AB][AC][BC]).
//
// All three pairwise margins are conditioned on:
//   y_ij. = sum_k y[i,j,k]   (row x column)
//   y_i.k = sum_j y[i,j,k]   (row x layer)
//   y_.jk = sum_i y[i,j,k]   (column x layer)
//
// Under the null of no three-way interaction, the conditional
// distribution given these pairwise margins is
//   P(Y | margins) proportional to 1 / prod_{i,j,k} y_ijk!
// so the test statistic T(Y) = sum lgamma(y_ijk + 1) orders tables
// from least extreme (small T) to most extreme (large T).
//
// Path linearization (Rxc1 paper p.23, with handwritten correction):
//   l = r + (c-1)(R-1) + k(R-1)(C-1)       (1-indexed)
//   total free cells = (R-1)(C-1)(K-1)
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
//   - Larger R*C*K becomes intractable quickly (factorial growth in
//     reference set); this kernel is intended for small tables only.

#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <algorithm>

using namespace Rcpp;

static const double RXCK_TOL = 3.45254e-7;

namespace rxck {

struct State {
  int R, C, K;
  // Pairwise margins (computed from observed table):
  std::vector<int> RC; // row x col, indexed [i*C + j]
  std::vector<int> RK; // row x layer
  std::vector<int> CK; // col x layer
  // Running residuals during path traversal:
  std::vector<int> RC_resid;
  std::vector<int> RK_resid;
  std::vector<int> CK_resid;
  // Current table; flat R*C*K, indexed [(i*C + j)*K + k].
  std::vector<int> y;

  double T_obs;        // sum lgamma(y_obs + 1)
  double mass_total;   // accumulated sum exp(-T)
  double mass_extreme; // accumulated sum exp(-T) over { Y : T >= T_obs }
};

inline int yidx(const State& s, int i, int j, int k) {
  return (i * s.C + j) * s.K + k;
}
inline int rcidx(const State& s, int i, int j) { return i * s.C + j; }
inline int rkidx(const State& s, int i, int k) { return i * s.K + k; }
inline int ckidx(const State& s, int j, int k) { return j * s.K + k; }

double table_T(const State& s) {
  double sum = 0.0;
  for (int v : s.y) sum += std::lgamma(v + 1);
  return sum;
}

// Compute dependent cells from pairwise margins. Returns false if any
// resulting cell is negative or if a consistency check (over-determined
// cell values agreeing across margins) fails.
bool finish_table(State& s) {
  const int R1 = s.R - 1;
  const int C1 = s.C - 1;
  const int K1 = s.K - 1;

  // 1. Last layer for free (i, j) cells: y[i, j, K-1] from RC margin.
  for (int i = 0; i < R1; ++i) {
    for (int j = 0; j < C1; ++j) {
      int sum = 0;
      for (int k = 0; k < K1; ++k) sum += s.y[yidx(s, i, j, k)];
      int v = s.RC[rcidx(s, i, j)] - sum;
      if (v < 0) return false;
      s.y[yidx(s, i, j, K1)] = v;
    }
  }
  // 2. Last column for free (i, k) cells (k < K-1): from RK margin.
  for (int i = 0; i < R1; ++i) {
    for (int k = 0; k < K1; ++k) {
      int sum = 0;
      for (int j = 0; j < C1; ++j) sum += s.y[yidx(s, i, j, k)];
      int v = s.RK[rkidx(s, i, k)] - sum;
      if (v < 0) return false;
      s.y[yidx(s, i, C1, k)] = v;
    }
  }
  // 3. Last row for free (j, k) cells (j < C-1, k < K-1): from CK margin.
  for (int j = 0; j < C1; ++j) {
    for (int k = 0; k < K1; ++k) {
      int sum = 0;
      for (int i = 0; i < R1; ++i) sum += s.y[yidx(s, i, j, k)];
      int v = s.CK[ckidx(s, j, k)] - sum;
      if (v < 0) return false;
      s.y[yidx(s, R1, j, k)] = v;
    }
  }
  // 4. y[i, C-1, K-1] for i < R-1: from RC margin (sum over layers
  //    in the last column for row i must equal RC[i, C-1]).
  //    Cross-check with RK margin (sum over columns in last layer
  //    for row i must equal RK[i, K-1]).
  for (int i = 0; i < R1; ++i) {
    int sum_k = 0;
    for (int k = 0; k < K1; ++k) sum_k += s.y[yidx(s, i, C1, k)];
    int v = s.RC[rcidx(s, i, C1)] - sum_k;
    if (v < 0) return false;
    s.y[yidx(s, i, C1, K1)] = v;
    int rk_check = 0;
    for (int j = 0; j < s.C; ++j) rk_check += s.y[yidx(s, i, j, K1)];
    if (rk_check != s.RK[rkidx(s, i, K1)]) return false;
  }
  // 5. y[R-1, j, K-1] for j < C-1: from CK margin.
  //    Cross-check with RC margin.
  for (int j = 0; j < C1; ++j) {
    int sum_i = 0;
    for (int i = 0; i < R1; ++i) sum_i += s.y[yidx(s, i, j, K1)];
    int v = s.CK[ckidx(s, j, K1)] - sum_i;
    if (v < 0) return false;
    s.y[yidx(s, R1, j, K1)] = v;
    int rc_check = 0;
    for (int k = 0; k < s.K; ++k) rc_check += s.y[yidx(s, R1, j, k)];
    if (rc_check != s.RC[rcidx(s, R1, j)]) return false;
  }
  // 6. y[R-1, C-1, k] for k < K-1: from CK margin.
  //    Cross-check with RK margin.
  for (int k = 0; k < K1; ++k) {
    int sum_i = 0;
    for (int i = 0; i < R1; ++i) sum_i += s.y[yidx(s, i, C1, k)];
    int v = s.CK[ckidx(s, C1, k)] - sum_i;
    if (v < 0) return false;
    s.y[yidx(s, R1, C1, k)] = v;
    int rk_check = 0;
    for (int j = 0; j < s.C; ++j) rk_check += s.y[yidx(s, R1, j, k)];
    if (rk_check != s.RK[rkidx(s, R1, k)]) return false;
  }
  // 7. Final cell y[R-1, C-1, K-1]: any margin works; use RC.
  //    Cross-check with RK and CK.
  {
    int sum_k = 0;
    for (int k = 0; k < K1; ++k) sum_k += s.y[yidx(s, R1, C1, k)];
    int v = s.RC[rcidx(s, R1, C1)] - sum_k;
    if (v < 0) return false;
    s.y[yidx(s, R1, C1, K1)] = v;
    int rk_check = 0;
    for (int j = 0; j < s.C; ++j) rk_check += s.y[yidx(s, R1, j, K1)];
    if (rk_check != s.RK[rkidx(s, R1, K1)]) return false;
    int ck_check = 0;
    for (int i = 0; i < s.R; ++i) ck_check += s.y[yidx(s, i, C1, K1)];
    if (ck_check != s.CK[ckidx(s, C1, K1)]) return false;
  }
  return true;
}

void traverse(int l, int lmax, State& s) {
  if (l > lmax) {
    if (!finish_table(s)) return;
    double T = table_T(s);
    double w = std::exp(-T);
    s.mass_total += w;
    if (T - s.T_obs >= -RXCK_TOL) s.mass_extreme += w;
    return;
  }
  // Decode l (1-indexed) -> (r, c, k) zero-indexed via
  //   l_0 = r + c*(R-1) + k*(R-1)*(C-1)
  const int l0 = l - 1;
  const int R1 = s.R - 1;
  const int C1 = s.C - 1;
  const int per_layer = R1 * C1;
  const int k = l0 / per_layer;
  const int rem = l0 % per_layer;
  const int c = rem / R1;
  const int r = rem % R1;

  // Upper bound = min of the three pairwise residuals at this cell.
  // (These were decremented by previous placements at cells that share
  // each pairwise margin.)
  const int x_up = std::min({
    s.RC_resid[rcidx(s, r, c)],
    s.RK_resid[rkidx(s, r, k)],
    s.CK_resid[ckidx(s, c, k)]
  });
  if (x_up < 0) return;

  for (int x = 0; x <= x_up; ++x) {
    s.y[yidx(s, r, c, k)] = x;
    s.RC_resid[rcidx(s, r, c)] -= x;
    s.RK_resid[rkidx(s, r, k)] -= x;
    s.CK_resid[ckidx(s, c, k)] -= x;
    traverse(l + 1, lmax, s);
    s.RC_resid[rcidx(s, r, c)] += x;
    s.RK_resid[rkidx(s, r, k)] += x;
    s.CK_resid[ckidx(s, c, k)] += x;
  }
  s.y[yidx(s, r, c, k)] = 0;
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
  s.R = d[0];
  s.C = d[1];
  s.K = d[2];
  if (s.R < 2 || s.C < 2 || s.K < 2)
    Rcpp::stop("each dimension must be at least 2");

  const int total = s.R * s.C * s.K;
  s.y.assign(total, 0);
  s.RC.assign(s.R * s.C, 0);
  s.RK.assign(s.R * s.K, 0);
  s.CK.assign(s.C * s.K, 0);

  // Source layout is column-major (R's array storage):
  //   dat[i + R*j + R*C*k] -> cell (i, j, k).
  // Internal layout uses [(i*C + j)*K + k].
  for (int k = 0; k < s.K; ++k) {
    for (int j = 0; j < s.C; ++j) {
      for (int i = 0; i < s.R; ++i) {
        int v = dat[i + s.R * j + s.R * s.C * k];
        if (v < 0) Rcpp::stop("dat entries must be non-negative");
        s.y[rxck::yidx(s, i, j, k)] = v;
        s.RC[rxck::rcidx(s, i, j)] += v;
        s.RK[rxck::rkidx(s, i, k)] += v;
        s.CK[rxck::ckidx(s, j, k)] += v;
      }
    }
  }
  s.T_obs = rxck::table_T(s);

  // Reset y so the recursive traversal builds fresh tables.
  std::fill(s.y.begin(), s.y.end(), 0);
  s.RC_resid = s.RC;
  s.RK_resid = s.RK;
  s.CK_resid = s.CK;
  s.mass_total = 0.0;
  s.mass_extreme = 0.0;

  const int lmax = (s.R - 1) * (s.C - 1) * (s.K - 1);
  rxck::traverse(1, lmax, s);

  if (s.mass_total <= 0.0)
    Rcpp::stop("rxck enumeration produced no feasible tables; "
               "check that pairwise margins are mutually consistent");

  double pval = s.mass_extreme / s.mass_total;
  if (pval < 0.0) pval = 0.0;
  if (pval > 1.0) pval = 1.0;
  return pval;
}
