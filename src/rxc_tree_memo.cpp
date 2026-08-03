// rxc_tree_memo.cpp
// First-generation r x c Fisher's exact test kernel.
//
// Algorithm (from Rxc1 paper pp.21-28, network/path reduction):
//   - Linearize the (R-1)(C-1) free cells along the path
//     (1,1) -> (2,1) -> ... -> (R-1,1) -> (1,2) -> ... -> (R-1,C-1)
//     using k = r + (c-1)(R-1) (1-indexed in the paper; 0-indexed here).
//   - Recursively place each free cell over its feasible range
//     [max(0, R_resid - sum_{c'>c} C_resid[c']), min(R_resid, C_resid)].
//   - When the path completes, fill the last column (cells (r,C-1) =
//     remaining row residual) and last row (cells (R-1,c) = remaining
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

static const double RXC_TOL = 3.45254e-7;

namespace {

struct RxC {
  int R, C, N;
  std::vector<int> Rm, Cm;     // observed row/column margins
  std::vector<int> Rr, Cr;     // residual row/column margins
  std::vector<int> y;          // current table, R*C ints, row-major
  double cns;                  // log(prod R_i! prod C_j! / N!)
  double pobs;
  double pval;                 // antipvalue, starts at 1
};

inline int idx(const RxC& s, int r, int c) { return r * s.C + c; }

inline double table_p(const RxC& s) {
  double sum = 0.0;
  for (int v : s.y) sum += std::lgamma(v + 1);
  return std::exp(s.cns - sum);
}

// After path traversal completes, fill last column and last row from
// margin residuals. Returns false if any computed cell is negative.
//
// Implementation note: this function does NOT mutate s.Cr or s.Rr
// (only s.y). An earlier version decremented s.Cr[C-1] in place; on
// the infeasibility-return code path that mutation was not reverted
// and corrupted subsequent iterations. The current form computes the
// last (R-1, C-1) cell from a local sum and reads the last-row cells
// directly from s.Cr[c], which the caller has not yet mutated.
bool finish_table(RxC& s) {
  const int Rm1 = s.R - 1;
  const int Cm1 = s.C - 1;
  // Last column for rows 0..R-2: cell value = remaining row residual.
  int sum_last_col = 0;
  for (int r = 0; r < Rm1; ++r) {
    int v = s.Rr[r];
    if (v < 0) return false;
    s.y[idx(s, r, Cm1)] = v;
    sum_last_col += v;
  }
  // (R-1, C-1) cell from last-column total minus what was just placed.
  int v_corner = s.Cm[Cm1] - sum_last_col;
  if (v_corner < 0) return false;
  s.y[idx(s, Rm1, Cm1)] = v_corner;
  // Last row for cols 0..C-2: cell value = remaining column residual
  // (s.Cr[c] is the residual after path traversal in cols < C-1).
  for (int c = 0; c < Cm1; ++c) {
    int v = s.Cr[c];
    if (v < 0) return false;
    s.y[idx(s, Rm1, c)] = v;
  }
  return true;
}

void traverse(int k, int kmax, RxC& s) {
  if (k > kmax) {
    if (!finish_table(s)) return;
    double p = table_p(s);
    if (p - s.pobs > s.pobs * RXC_TOL) s.pval -= p;
    return;
  }
  // Path step k (1-indexed) -> cell (r, c) zero-indexed.
  const int kk = k - 1;
  const int r = kk % (s.R - 1);
  const int c = kk / (s.R - 1);

  // Lower bound: row r residual must be absorbable by remaining-cols
  // capacity (cols c+1..C-1), so x >= R_resid - sum of those C_resid.
  int sum_rest_cols = 0;
  for (int cc = c + 1; cc < s.C; ++cc) sum_rest_cols += s.Cr[cc];
  const int x_lo = std::max(0, s.Rr[r] - sum_rest_cols);
  const int x_up = std::min(s.Rr[r], s.Cr[c]);

  for (int x = x_lo; x <= x_up; ++x) {
    s.y[idx(s, r, c)] = x;
    s.Rr[r] -= x;
    s.Cr[c] -= x;
    traverse(k + 1, kmax, s);
    s.Rr[r] += x;
    s.Cr[c] += x;
  }
  s.y[idx(s, r, c)] = 0;
}

} // namespace

// [[Rcpp::export(name = ".rxc_tree_memo_cpp")]]
double rxc_tree_memo_cpp(IntegerMatrix dat) {
  RxC s;
  s.R = dat.nrow();
  s.C = dat.ncol();
  if (s.R < 2) Rcpp::stop("dat must have at least 2 rows");
  if (s.C < 2) Rcpp::stop("dat must have at least 2 columns");

  s.Rm.assign(s.R, 0);
  s.Cm.assign(s.C, 0);
  s.y.assign(s.R * s.C, 0);
  for (int i = 0; i < s.R; ++i) {
    for (int j = 0; j < s.C; ++j) {
      const int v = dat(i, j);
      if (v < 0) Rcpp::stop("dat entries must be non-negative");
      s.Rm[i] += v;
      s.Cm[j] += v;
    }
  }
  s.N = 0;
  for (int v : s.Rm) s.N += v;

  double sumLogR = 0.0, sumLogC = 0.0;
  for (int v : s.Rm) sumLogR += std::lgamma(v + 1);
  for (int v : s.Cm) sumLogC += std::lgamma(v + 1);
  s.cns = sumLogR + sumLogC - std::lgamma(s.N + 1);

  double sumLogObs = 0.0;
  for (int i = 0; i < s.R; ++i)
    for (int j = 0; j < s.C; ++j)
      sumLogObs += std::lgamma(dat(i, j) + 1);
  s.pobs = std::exp(s.cns - sumLogObs);

  s.Rr = s.Rm;
  s.Cr = s.Cm;
  s.pval = 1.0;

  const int kmax = (s.R - 1) * (s.C - 1);
  traverse(1, kmax, s);

  // Numeric guard: pval should be in (0, 1]. Clamp to handle round-off
  // around the boundary cases.
  if (s.pval < 0.0) s.pval = 0.0;
  if (s.pval > 1.0) s.pval = 1.0;
  return s.pval;
}
