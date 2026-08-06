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
// Pruning: feasibility-only. Joe/Requena majorization-based max/min
// pruning (analog of rx2_tree_memo's find_max/find_min) is deferred to
// a second generation; it requires a multi-column hypergeometric mode
// formula (the 9/13/25 research challenge note).

#include <Rcpp.h>
#include <vector>
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
  long n_skipped;                    // subtrees skipped by the bound
};

inline int yidx(const RxcState& s, int i, int j) { return i * s.c + j; }

inline double table_p(const RxcState& s) {
  double sum = 0.0;
  for (int v : s.y) sum += std::lgamma(v + 1);
  return std::exp(s.log_const - sum);
}

// --- Acceptance-region skip via a balanced-split relaxation ------
//
// The m x 2 kernel prunes with exact subtree extrema, which exist
// there because two columns make log Pr separable into concave
// one-dimensional terms. That separability is unavailable for
// c >= 3: the cells of a partial table are tied by both row and
// column constraints simultaneously.
//
// What does transfer is a relaxation. For q free cells summing to
// s, dropping all constraints except the total gives
//
//   max { -sum_t log(y_t!) : sum y_t = s, y_t >= 0 }
//     = -(q - rho) * log(floor(s/q)!) - rho * log(ceil(s/q)!),
//     rho = s mod q,
//
// attained at the balanced split (proof: -log(y!) is concave, so
// moving a unit from a cell exceeding another by two or more
// strictly increases the objective; hence maxima have all cells
// within one of each other). Relaxing constraints only enlarges
// the feasible set, so applying this per remaining column yields a
// valid upper bound on the log-probability of every completion.
// When that bound does not exceed log(Pr_obs * (1 + tol)), the
// entire subtree lies in the acceptance region and contributes
// nothing to the antipvalue, so it may be skipped exactly.
//
// Unlike the m x 2 case the bound is not the true extremum, so no
// interval structure over siblings is claimed and no scan cutoff
// is applied. The complement side has no analogue at all: bulk
// subtraction there needs a subtree's total mass, which for c = 2
// is Vandermonde's convolution but for c >= 3 has no closed form.
// See docs/c2_separability_whitepaper_2026-08-04.md, Section 5.
inline double balanced_split_lp(const RxcState& s, int q, int total) {
  if (q <= 0) return (total == 0) ? 0.0 : -INFINITY;
  const int base = total / q, rem = total % q;
  return -((double)(q - rem) * s.lf[base] +
           (double)rem * s.lf[base + 1]);
}

// Upper bound on log Pr over all completions of the current state,
// given the log-mass already committed to placed and determined
// cells. Column j has (m - 1 - i) free cells left in this pass
// plus its determined last-row cell, which the relaxation absorbs
// by treating the column residual as spread over m - i cells;
// later columns are unconstrained apart from their margins.
double completion_upper_lp(const RxcState& s, int i, int j,
                           double placed_lp) {
  double bound = placed_lp;
  bound += balanced_split_lp(s, s.m - i, s.cresid[j]);
  for (int jj = j + 1; jj < s.c; ++jj)
    bound += balanced_split_lp(s, s.m, s.cresid[jj]);
  return s.log_const + bound;
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

void traverse(int t, int tmax, RxcState& s, double placed_lp) {
  if (t > tmax) {
    if (!finish_table(s)) return;
    double p = table_p(s);
    if (p - s.p_obs > s.p_obs * RXC_TOL) s.pval -= p;
    return;
  }
  // Path step t (1-indexed) -> cell (i, j) zero-indexed.
  const int tt = t - 1;
  const int i = tt % (s.m - 1);
  const int j = tt / (s.m - 1);

  // Acceptance-region skip: if no completion can exceed the
  // threshold, the whole subtree contributes nothing.
  if (completion_upper_lp(s, i, j, placed_lp) <= s.log_thresh) {
    ++s.n_skipped;
    return;
  }

  // Lower bound: row i residual must be absorbable by remaining-cols
  // capacity (cols j+1..c-1), so x >= rresid - sum of those cresid.
  int sum_rest_cols = 0;
  for (int jj = j + 1; jj < s.c; ++jj) sum_rest_cols += s.cresid[jj];
  const int x_lo = std::max(0, s.rresid[i] - sum_rest_cols);
  const int x_up = std::min(s.rresid[i], s.cresid[j]);

  for (int x = x_lo; x <= x_up; ++x) {
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
  traverse(1, tmax, s, 0.0);

  // Numeric guard: pval should be in (0, 1]. Clamp to handle round-off
  // around the boundary cases.
  if (s.pval < 0.0) s.pval = 0.0;
  if (s.pval > 1.0) s.pval = 1.0;
  return s.pval;
}
