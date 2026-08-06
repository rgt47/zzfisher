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
// Scope notes:
//   - Because the p-value is a ratio of two sums over the full
//     reference set, significance-region skipping cannot reduce the
//     traversal (every feasible table feeds the denominator); the
//     bound machinery of the m x 2 and r x c kernels does not apply
//     here. What does transfer from the r x c pass is
//     determine-as-you-go: a layer's dependent last-column and
//     last-row cells are fully determined once its free cells are
//     placed, so they are computed, feasibility-checked, and charged
//     to the AB margins at the layer boundary instead of at the
//     leaf. Infeasible or inconsistent prefixes are pruned entire
//     layers early, and the test statistic accumulates as a running
//     prefix over a log-factorial table (no per-leaf lgamma sweep
//     over all m*c*k cells).
//   - Validated against a pure-R brute-force enumerator for the 2x2x2
//     case in inst/tinytest/test_htest.R (one free cell, easy oracle),
//     and by fingerprint against the leaf-check predecessor over
//     random small arrays.
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

  double T_obs;        // sum log(y_obs!)
  double mass_total;   // accumulated sum exp(-T)
  double mass_extreme; // accumulated sum exp(-T) over { Y : T >= T_obs }
  std::vector<double> lf;  // lf[v] = log(v!)
  // Per-layer close scratch (a layer is closed at most once per
  // active path, so per-layer buffers cannot alias; heap churn per
  // close measured a 25 percent regression before this).
  std::vector<std::vector<int>> close_ic, close_mj;
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

// Close free layer ll: its dependent last-column cells (i, c-1, ll),
// last-row cells (m-1, j, ll), and corner (m-1, c-1, ll) are fully
// determined by the layer's AC and BC residuals. Feasibility (the
// dependent values fit their AB margins) and the layer's
// over-determination consistency (the corner agrees from both the
// BC and AC sides) are checked BEFORE any mutation; on success the
// dependent values are charged to the AB residuals — the only
// margins later stages read — and their log-factorials summed into
// t_add. Returns false (nothing mutated) if the layer is
// infeasible, pruning the prefix k-2 layers before the old
// leaf-time finish_table() could.
bool close_layer(State& s, int ll, double& t_add, int& v_corner) {
  const int mm1 = s.m - 1;
  const int cm1 = s.c - 1;
  std::vector<int>& v_ic = s.close_ic[ll];
  std::vector<int>& v_mj = s.close_mj[ll];
  int sum_ic = 0, sum_mj = 0;
  for (int i = 0; i < mm1; ++i) {
    const int v = s.ac_resid[acidx(s, i, ll)];
    if (s.ab_resid[abidx(s, i, cm1)] < v) return false;
    v_ic[i] = v;
    sum_ic += v;
  }
  for (int j = 0; j < cm1; ++j) {
    const int v = s.bc_resid[bcidx(s, j, ll)];
    if (s.ab_resid[abidx(s, mm1, j)] < v) return false;
    v_mj[j] = v;
    sum_mj += v;
  }
  v_corner = s.bc_resid[bcidx(s, cm1, ll)] - sum_ic;
  if (v_corner < 0) return false;
  // Corner over-determination: BC and AC sides must agree.
  if (s.ac_resid[acidx(s, mm1, ll)] - sum_mj != v_corner) return false;
  if (s.ab_resid[abidx(s, mm1, cm1)] < v_corner) return false;

  t_add = 0.0;
  for (int i = 0; i < mm1; ++i) {
    s.ab_resid[abidx(s, i, cm1)] -= v_ic[i];
    s.y[yidx(s, i, cm1, ll)] = v_ic[i];
    t_add += s.lf[v_ic[i]];
  }
  for (int j = 0; j < cm1; ++j) {
    s.ab_resid[abidx(s, mm1, j)] -= v_mj[j];
    s.y[yidx(s, mm1, j, ll)] = v_mj[j];
    t_add += s.lf[v_mj[j]];
  }
  s.ab_resid[abidx(s, mm1, cm1)] -= v_corner;
  s.y[yidx(s, mm1, cm1, ll)] = v_corner;
  t_add += s.lf[v_corner];
  return true;
}

void open_layer(State& s, int ll, int v_corner) {
  const int mm1 = s.m - 1;
  const int cm1 = s.c - 1;
  const std::vector<int>& v_ic = s.close_ic[ll];
  const std::vector<int>& v_mj = s.close_mj[ll];
  for (int i = 0; i < mm1; ++i)
    s.ab_resid[abidx(s, i, cm1)] += v_ic[i];
  for (int j = 0; j < cm1; ++j)
    s.ab_resid[abidx(s, mm1, j)] += v_mj[j];
  s.ab_resid[abidx(s, mm1, cm1)] += v_corner;
}

// After every free layer is placed and closed, the last layer is
// determined outright: y[i, j, k-1] = ab_resid[i, j] for all
// (i, j). Non-negativity is automatic (placements and closes never
// drive an AB residual negative); what can still fail is the
// over-determination against the last layer's AC and BC margins,
// which the closes never touched. On success the completed table's
// statistic is t_pref plus the last layer's log-factorials.
void leaf(State& s, double t_pref) {
  const int km1 = s.k - 1;
  for (int i = 0; i < s.m; ++i) {
    int row_sum = 0;
    for (int j = 0; j < s.c; ++j)
      row_sum += s.ab_resid[abidx(s, i, j)];
    if (row_sum != s.ac_resid[acidx(s, i, km1)]) return;
  }
  for (int j = 0; j < s.c; ++j) {
    int col_sum = 0;
    for (int i = 0; i < s.m; ++i)
      col_sum += s.ab_resid[abidx(s, i, j)];
    if (col_sum != s.bc_resid[bcidx(s, j, km1)]) return;
  }
  double T = t_pref;
  for (int i = 0; i < s.m; ++i)
    for (int j = 0; j < s.c; ++j)
      T += s.lf[s.ab_resid[abidx(s, i, j)]];
  const double w = std::exp(-T);
  s.mass_total += w;
  if (T - s.T_obs >= -RXCK_TOL) s.mass_extreme += w;
}

void traverse_cell(int t, int tmax, State& s, double t_pref,
                   int i, int j, int l);

// t_pref carries log(x!) of every placed free cell plus the
// contributions of every closed layer. A layer is closed when the
// path reaches the first cell of the next layer (or the end of the
// path), so infeasible prefixes die at layer boundaries.
void traverse(int t, int tmax, State& s, double t_pref) {
  const int mm1 = s.m - 1;
  const int cm1 = s.c - 1;
  const int per_layer = mm1 * cm1;

  if (t > tmax) {
    double t_add = 0.0;
    int v_corner = 0;
    if (!close_layer(s, s.k - 2, t_add, v_corner)) return;
    leaf(s, t_pref + t_add);
    open_layer(s, s.k - 2, v_corner);
    return;
  }
  // Decode t (1-indexed) -> (i, j, l) zero-indexed via
  //   t_0 = i + j*(m-1) + l*(m-1)*(c-1)
  const int t0 = t - 1;
  const int l = t0 / per_layer;
  const int rem = t0 % per_layer;
  const int j = rem / mm1;
  const int i = rem % mm1;

  // Layer boundary: close the completed previous layer first.
  if (i == 0 && j == 0 && l > 0) {
    double t_add = 0.0;
    int v_corner = 0;
    if (!close_layer(s, l - 1, t_add, v_corner)) return;
    traverse_cell(t, tmax, s, t_pref + t_add, i, j, l);
    open_layer(s, l - 1, v_corner);
    return;
  }
  traverse_cell(t, tmax, s, t_pref, i, j, l);
}

void traverse_cell(int t, int tmax, State& s, double t_pref,
                   int i, int j, int l) {
  // Upper bound = min of the three pairwise residuals at this cell.
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
    traverse(t + 1, tmax, s, t_pref + s.lf[x]);
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

  int n_total = 0;
  for (int v : s.y) n_total += v;
  s.lf.assign(n_total + 2, 0.0);
  for (int v = 1; v <= n_total + 1; ++v)
    s.lf[v] = s.lf[v - 1] + std::log((double)v);
  s.close_ic.assign(s.k, std::vector<int>(s.m - 1, 0));
  s.close_mj.assign(s.k, std::vector<int>(s.c - 1, 0));

  // Reset y so the recursive traversal builds fresh tables.
  std::fill(s.y.begin(), s.y.end(), 0);
  s.ab_resid = s.ab;
  s.ac_resid = s.ac;
  s.bc_resid = s.bc;
  s.mass_total = 0.0;
  s.mass_extreme = 0.0;

  const int tmax = (s.m - 1) * (s.c - 1) * (s.k - 1);
  rxck::traverse(1, tmax, s, 0.0);

  if (s.mass_total <= 0.0)
    Rcpp::stop("rxck enumeration produced no feasible tables; "
               "check that pairwise margins are mutually consistent");

  double pval = s.mass_extreme / s.mass_total;
  if (pval < 0.0) pval = 0.0;
  if (pval > 1.0) pval = 1.0;
  return pval;
}
