// rxck_net_ci.cpp
// Exact test of CONDITIONAL INDEPENDENCE in an m x c x k array:
// rows independent of columns given the layer (the log-linear
// model [AC][BC]), the stratified Fisher's exact test.
//
// Structure. Conditioning on each stratum's row and column
// margins makes the conditional law a PRODUCT of independent
// per-stratum r x c multivariate hypergeometrics, so
// log P(Y) = sum_l log P_l(Y_l) and the reference set is the
// product of the per-stratum fibers. The algorithm is therefore
// the sweep architecture with strata as stages and no state
// vector at all:
//   1. Per stratum, enumerate its fiber once, collecting the
//      distribution of log P_l as a class list (value, mass)
//      hashed within 1e-12 (floating-point ties only). Each
//      stratum's masses sum to one.
//   2. Suffix extrema of the total log-probability are exact by
//      simple sums of the per-stratum minima and maxima.
//   3. Sweep the strata, convolving class lists with three-case
//      classification against the exact suffix extrema: a class
//      certified inside the complement region U = {P > P_obs
//      (1+tol)} contributes its own mass in O(1), because the
//      remaining strata's conditional masses sum to one (the
//      product-form Pagano-Halvorsen identity); a class
//      certified in the significance region is dropped; only
//      straddling classes convolve forward.
//   4. p = 1 - mass(U), the minimum-likelihood two-sided value
//      with the FEXACT tolerance convention.
// This is the stage-wise full-distribution discipline of
// Mehta-Patel-Gray (1985) fused with the classification devices
// of the two-way kernels. With k = 1 the test reduces to the
// ordinary r x c Fisher exact test, which the test suite checks
// against fisher.test().
//
// Scaling note: each stratum's fiber is enumerated completely
// (its distribution is needed), so strata must individually be
// small; the cross-stratum work is convolution of class lists,
// polynomial in their sizes, guarded by MAX_CLASSES.

#include <Rcpp.h>
#include <vector>
#include <unordered_map>
#include <cmath>

using namespace Rcpp;

namespace {

struct CiCls {
  double q;   // log total probability value of the class members
  double w;   // ABSOLUTE probability mass of the class
};

using CiMap = std::unordered_map<long long, CiCls>;

inline long long ci_qkey(double q) { return llround(q * 1e12); }

// Members of a class share q to within 1e-12, so masses add
// directly.
inline void ci_add(CiMap& m, double q, double w) {
  auto r = m.emplace(ci_qkey(q), CiCls{q, w});
  if (!r.second) r.first->second.w += w;
}

struct CiState {
  std::vector<double> lf;
};

// Enumerate one stratum's r x c fiber, adding (log P_l, mass)
// classes. Free cells (i, j), i < m-1, j < c-1, column-major;
// last column and row determined by residual margins.
struct StratumEnum {
  const CiState* st;
  int m, c;
  std::vector<int> rres, cres, cmarg;
  double log_const;
  CiMap* out;

  void run(int i, int j, double acc) {
    if (j >= c - 1) {
      double det = 0.0;
      int corner = cmarg[c - 1];
      for (int g = 0; g < m - 1; ++g) {
        if (rres[g] < 0) return;
        det += st->lf[rres[g]];
        corner -= rres[g];
      }
      if (corner < 0) return;
      det += st->lf[corner];
      for (int h = 0; h < c - 1; ++h) {
        if (cres[h] < 0) return;
        det += st->lf[cres[h]];
      }
      const double q = log_const + acc - det;
      ci_add(*out, q, std::exp(q));
      return;
    }
    int in = i + 1, jn = j;
    if (in > m - 2) { in = 0; jn = j + 1; }
    int rest = 0;
    for (int h = j + 1; h < c; ++h) rest += cres[h];
    const int lo = std::max(0, rres[i] - rest);
    const int hi = std::min(rres[i], cres[j]);
    for (int x = lo; x <= hi; ++x) {
      rres[i] -= x;
      cres[j] -= x;
      run(in, jn, acc - st->lf[x]);
      rres[i] += x;
      cres[j] += x;
    }
  }
};

}  // anonymous namespace

// [[Rcpp::export(name = ".rxck_net_ci_cpp")]]
double rxck_net_ci_cpp(IntegerVector dat) {
  RObject dim_attr = dat.attr("dim");
  if (dim_attr.isNULL())
    Rcpp::stop("dat must be a 3D integer array");
  IntegerVector d(dim_attr);
  if (d.size() != 3)
    Rcpp::stop("dat must be a 3D integer array");
  const int m = d[0], c = d[1], k = d[2];
  if (m < 2 || c < 2 || k < 1)
    Rcpp::stop("dat must be at least 2 x 2 x 1");

  int n_total = 0;
  for (int t = 0; t < dat.size(); ++t) {
    if (dat[t] < 0) Rcpp::stop("dat entries must be non-negative");
    n_total += dat[t];
  }
  CiState st;
  st.lf.assign(n_total + 2, 0.0);
  for (int v = 1; v <= n_total + 1; ++v)
    st.lf[v] = st.lf[v - 1] + std::log((double)v);

  // Per-stratum class lists, observed log-probabilities, extrema.
  std::vector<std::vector<CiCls>> strat(k);
  std::vector<double> smin(k), smax(k);
  double log_p_obs = 0.0;
  for (int l = 0; l < k; ++l) {
    StratumEnum en;
    en.st = &st;
    en.m = m;
    en.c = c;
    en.rres.assign(m, 0);
    en.cres.assign(c, 0);
    int nl = 0;
    double obs_lf = 0.0;
    for (int i = 0; i < m; ++i)
      for (int j = 0; j < c; ++j) {
        const int v = dat[i + m * j + m * c * l];
        en.rres[i] += v;
        en.cres[j] += v;
        nl += v;
        obs_lf += st.lf[v];
      }
    en.cmarg = en.cres;
    en.log_const = -st.lf[nl];
    for (int i = 0; i < m; ++i) en.log_const += st.lf[en.rres[i]];
    for (int j = 0; j < c; ++j) en.log_const += st.lf[en.cres[j]];
    log_p_obs += en.log_const - obs_lf;

    CiMap cm;
    en.out = &cm;
    if (m == 1 || nl == 0) {
      ci_add(cm, 0.0, 1.0);   // degenerate stratum: single table
    } else {
      en.run(0, 0, 0.0);
    }
    std::vector<CiCls>& cl = strat[l];
    cl.reserve(cm.size());
    double mn = INFINITY, mx = -INFINITY;
    for (const auto& kv : cm) {
      cl.push_back(kv.second);
      if (kv.second.q < mn) mn = kv.second.q;
      if (kv.second.q > mx) mx = kv.second.q;
    }
    smin[l] = mn;
    smax[l] = mx;
  }
  const double log_thresh = log_p_obs + std::log1p(3.45254e-7);

  // Exact suffix extrema across strata.
  std::vector<double> sufmin(k + 1, 0.0), sufmax(k + 1, 0.0);
  for (int l = k - 1; l >= 0; --l) {
    sufmin[l] = sufmin[l + 1] + smin[l];
    sufmax[l] = sufmax[l + 1] + smax[l];
  }

  // Sweep the strata: convolve with three-case classification.
  const size_t MAX_CLASSES = 5000000;
  CiMap cur, nxt;
  ci_add(cur, 0.0, 1.0);
  double mass_U = 0.0;
  for (int l = 0; l < k && !cur.empty(); ++l) {
    nxt.clear();
    for (const auto& kv : cur) {
      const CiCls& e = kv.second;
      if (e.q + sufmax[l] <= log_thresh) continue;   // all in S
      if (e.q + sufmin[l] > log_thresh) {            // all in U
        mass_U += e.w;   // remaining strata's masses sum to one
        continue;
      }
      for (const CiCls& f : strat[l])
        ci_add(nxt, e.q + f.q, e.w * f.w);
    }
    if (nxt.size() > MAX_CLASSES)
      Rcpp::stop("class limit exceeded in stratified CI sweep");
    cur.swap(nxt);
  }
  // Classes surviving all strata are complete tables.
  for (const auto& kv : cur)
    if (kv.second.q > log_thresh) mass_U += kv.second.w;

  double pval = 1.0 - mass_U;
  if (pval < 0.0) pval = 0.0;
  if (pval > 1.0) pval = 1.0;
  return pval;
}
