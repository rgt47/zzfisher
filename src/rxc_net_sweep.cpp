// rxc_net_sweep.cpp
// Fused state-aggregation kernel for Fisher's exact test on
// m x c tables: the Mehta-Patel network with exact backward
// extrema, the class-merge device executed by sorting, and
// closed-form bulk disposal, on the complement-region ledger.
// This is the general-c form of the m x 2 sweep kernel
// (rx2_net_sweep.cpp) and the headline construction of the
// flagship r x c paper.
//
// Design, with attributions:
//   - Stages are rows; a state is the vector of column residuals
//     after the placed rows. Two prefixes with equal residuals
//     are one state carrying summed mass (Mehta-Patel 1983).
//     The table is TRANSPOSED first if it has more columns than
//     rows, so the state dimension is min(m, c) - 1; this is
//     what makes the next item affordable.
//   - Completion extrema are computed EXACTLY, by backward
//     induction over the enumerated state graph (the device of
//     Mehta-Patel 1980, where the state was a scalar; here it is
//     a short vector). No searched bounds, no relaxations: every
//     classification uses the true maximum and minimum.
//   - Paths reaching a state with equal accumulated log-mass
//     merge into one class (Clarkson-Fan-Joe); classes are kept
//     sorted and merged within 1e-12 on the log scale, so only
//     floating-point ties merge and exactness is preserved.
//   - A class certified interior to the complement region
//     U = {P > P_obs(1+tol)} is disposed of in one step: the
//     total path weight of a state's completions is the
//     multinomial coefficient S!/prod(v_j!) (the many-column
//     Vandermonde; Mehta-Patel 1983, eq. 2.4). A class certified
//     inside the significance region is dropped. Only straddling
//     classes propagate.
//   - The p-value is accumulated as 1 - mass(U), so everything
//     discarded is discarded without tolerance bookkeeping.
//
// Storage grows with the live state and class count; there is no
// workspace parameter. Not yet done (candidates for a second
// pass): permutation-collapsing of tied-margin states, the
// child-creation classification of the m x 2 kernel, and a
// last-two-rows closed form.

#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <functional>
#include <cmath>

using namespace Rcpp;

namespace {

struct Cls {
  double q;   // accumulated log path weight: sum lr_g - sum l(x_gj)
  double w;   // linear multiplicity relative to q
};

struct FusedState {
  int m, c, n;
  std::vector<int> r;        // row margins, in processing order
  std::vector<int> cm;       // column margins
  std::vector<int> suffix_r; // suffix_r[i] = sum_{g>=i} r[g]
  double log_D;              // log(prod c_j! / n!)
  double log_thresh;         // log(p_obs (1 + tol))
  std::vector<double> lf;    // lf[v] = log v!
  double pval;
};

// Mixed-radix key machinery for residual vectors: radix c_j + 1,
// weight W[j] for coordinate j. A child key is the parent key
// minus sum_j x_j W[j], computed incrementally during
// enumeration, so no child vector is ever materialized on the
// hot path.
struct KeyCoder {
  std::vector<long long> W;
  std::vector<int> radix;
  void init(const std::vector<int>& cm) {
    const int c = (int)cm.size();
    W.assign(c, 1);
    radix.assign(c, 0);
    for (int j = 0; j < c; ++j) radix[j] = cm[j] + 1;
    for (int j = c - 2; j >= 0; --j) W[j] = W[j + 1] * radix[j + 1];
  }
  long long encode(const std::vector<int>& v) const {
    long long k = 0;
    for (size_t j = 0; j < v.size(); ++j) k += v[j] * W[j];
    return k;
  }
  void decode(long long k, std::vector<int>& v) const {
    for (size_t j = 0; j < W.size(); ++j) {
      v[j] = (int)(k / W[j]);
      k -= (long long)v[j] * W[j];
    }
  }
};

// Enumerate compositions x of `total` over c parts with
// 0 <= x_j <= v[j], calling f(key_delta, neg_lfac) where
// key_delta = sum_j x_j W[j] and neg_lfac = -sum_j l(x_j).
// Template recursion: no std::function, no allocation.
template <typename F>
void comps_rec(const FusedState& s, const int* v,
               const long long* W, const int* sufcap,
               int j, int rest, double acc, long long kd, F& f) {
  if (j == s.c - 1) {
    if (rest <= v[j])
      f(kd + rest * W[j], acc - s.lf[rest]);
    return;
  }
  const int lo = std::max(0, rest - sufcap[j + 1]);
  const int hi = std::min(v[j], rest);
  for (int t = lo; t <= hi; ++t)
    comps_rec(s, v, W, sufcap, j + 1, rest - t,
              acc - s.lf[t], kd + (long long)t * W[j], f);
}

template <typename F>
void compositions(const FusedState& s, const std::vector<int>& v,
                  const KeyCoder& kc, int total, F&& f) {
  int sufcap[16];
  sufcap[s.c] = 0;
  for (int j = s.c - 1; j >= 0; --j)
    sufcap[j] = sufcap[j + 1] + v[j];
  comps_rec(s, v.data(), kc.W.data(), sufcap, 0, total, 0.0,
            0LL, f);
}

// One stage of the state graph.
struct Stage {
  std::unordered_map<long long, int> idx;
  std::vector<std::vector<int>> vecs;
  std::vector<long long> keys;
  std::vector<double> vmax, vmin;      // exact completion extrema
  std::vector<std::vector<Cls>> cls;   // forward class lists
  int add(long long k, const KeyCoder& kc, int c) {
    auto it = idx.find(k);
    if (it != idx.end()) return it->second;
    int id = (int)vecs.size();
    idx.emplace(k, id);
    std::vector<int> v(c);
    kc.decode(k, v);
    vecs.push_back(std::move(v));
    keys.push_back(k);
    return id;
  }
};

// Class container per state: classes hashed by quantized log
// mass (grid 1e-12, so only floating-point ties merge and
// exactness is preserved). Profiling showed sorting the class
// lists — whose only purpose was this merge — took 39 percent of
// the run time at general c, where distinct path products
// abound; hashing is the device FEXACT always used here, and the
// m x 2 kernels could avoid it only because two columns produce
// few classes.
using ClsMap = std::unordered_map<long long, Cls>;

inline long long qkey(double q) {
  return llround(q * 1e12);
}

inline void cls_add(ClsMap& m, double q, double w) {
  auto r = m.emplace(qkey(q), Cls{q, w});
  if (!r.second) {
    Cls& e = r.first->second;
    e.w += w * std::exp(q - e.q);
  }
}

}  // anonymous namespace

// [[Rcpp::export(name = ".rxc_net_sweep_cpp")]]
double rxc_net_sweep_cpp(IntegerMatrix dat_in) {
  // Transpose so the state dimension (columns) is the smaller side.
  IntegerMatrix dat = dat_in;
  if (dat.ncol() > dat.nrow()) {
    IntegerMatrix t(dat_in.ncol(), dat_in.nrow());
    for (int i = 0; i < dat_in.nrow(); ++i)
      for (int j = 0; j < dat_in.ncol(); ++j)
        t(j, i) = dat_in(i, j);
    dat = t;
  }

  FusedState s;
  s.m = dat.nrow();
  s.c = dat.ncol();
  if (s.m < 2 || s.c < 2)
    Rcpp::stop("dat must be at least 2 x 2");
  if (s.c >= 16)
    Rcpp::stop("state dimension >= 16 unsupported");

  std::vector<int> r_raw(s.m, 0);
  s.cm.assign(s.c, 0);
  s.n = 0;
  for (int i = 0; i < s.m; ++i)
    for (int j = 0; j < s.c; ++j) {
      const int v = dat(i, j);
      if (v < 0) Rcpp::stop("dat entries must be non-negative");
      r_raw[i] += v;
      s.cm[j] += v;
      s.n += v;
    }

  // Process rows in descending margin order (residuals shrink
  // early, which keeps later stages' state counts down), and
  // sort columns descending so tied column margins are
  // contiguous: states differing only by a permutation within a
  // tied group have isomorphic completion structure and are
  // collapsed to a canonical representative (Mehta-Patel 1983,
  // implementation device 2).
  std::vector<int> ord(s.m);
  for (int i = 0; i < s.m; ++i) ord[i] = i;
  std::sort(ord.begin(), ord.end(),
            [&](int a, int b) { return r_raw[a] > r_raw[b]; });
  s.r.resize(s.m);
  for (int i = 0; i < s.m; ++i) s.r[i] = r_raw[ord[i]];
  std::sort(s.cm.begin(), s.cm.end(), std::greater<int>());

  s.suffix_r.assign(s.m + 1, 0);
  for (int i = s.m - 1; i >= 0; --i)
    s.suffix_r[i] = s.suffix_r[i + 1] + s.r[i];

  s.lf.assign(s.n + 2, 0.0);
  for (int v = 1; v <= s.n + 1; ++v)
    s.lf[v] = s.lf[v - 1] + std::log((double)v);

  s.log_D = -s.lf[s.n];
  for (int j = 0; j < s.c; ++j) s.log_D += s.lf[s.cm[j]];

  double log_p_obs = s.log_D;
  for (int i = 0; i < s.m; ++i) log_p_obs += s.lf[r_raw[i]];
  for (int i = 0; i < s.m; ++i)
    for (int j = 0; j < s.c; ++j)
      log_p_obs -= s.lf[dat(i, j)];
  s.log_thresh = log_p_obs + std::log1p(3.45254e-7);

  KeyCoder kc;
  kc.init(s.cm);
  long long KS = 1;
  for (int j = 0; j < s.c; ++j) {
    KS *= (s.cm[j] + 1);
    if (KS > (1LL << 26))
      Rcpp::stop("state key space too large for the dense sweep");
  }

  // Tied-column groups (columns sorted descending, so ties are
  // contiguous runs). A state is canonical when its residuals
  // are non-increasing within every tied group; canon() sorts a
  // state into that form. For tables with no ties canon is the
  // identity and costs one flag test.
  std::vector<int> grp_start(s.c);
  bool has_ties = false;
  for (int j = 0; j < s.c; ++j) {
    if (j > 0 && s.cm[j] == s.cm[j - 1]) {
      grp_start[j] = grp_start[j - 1];
      has_ties = true;
    } else {
      grp_start[j] = j;
    }
  }
  std::vector<int> cbuf(s.c);
  auto canon = [&](long long k) -> long long {
    if (!has_ties) return k;
    kc.decode(k, cbuf);
    for (int j = 0; j < s.c; ) {
      int e = j + 1;
      while (e < s.c && grp_start[e] == grp_start[j]) ++e;
      if (e - j > 1)
        std::sort(cbuf.begin() + j, cbuf.begin() + e,
                  std::greater<int>());
      j = e;
    }
    return kc.encode(cbuf);
  };
  auto is_canon = [&](const std::vector<int>& v) -> bool {
    if (!has_ties) return true;
    for (int j = 1; j < s.c; ++j)
      if (grp_start[j] == grp_start[j - 1] && v[j] > v[j - 1])
        return false;
    return true;
  };

  // Every bounded vector summing to the stage total is a
  // reachable state (transportation polytopes with matching
  // totals are never empty), so canonical states are enumerated
  // directly per stage. Arc handling is hybrid, chosen by
  // measurement: for tables WITH tied margins, arcs are
  // enumerated once into flat CSR storage with the child
  // canonicalized at build time, so the backward pass and every
  // straddling class scan flat arrays (canonicalizing per
  // class-arc instead costs the sort each time); for tie-free
  // tables, arcs are re-enumerated on the fly, because a CSR
  // store measured 2.3x slower there (tens of millions of stored
  // arcs cost more in memory traffic than the recursion costs in
  // arithmetic) and canon is the identity anyway.
  std::vector<std::vector<long long>> cst(s.m);
  std::vector<std::vector<int>> coff(s.m);
  std::vector<std::vector<long long>> arc_ck(s.m);
  std::vector<std::vector<double>> arc_st(s.m);
  std::vector<std::unordered_map<long long, int>> sidx(s.m);
  std::vector<int> vbuf(s.c);
  if (has_ties) {
    const long long MAX_ARCS = 1LL << 26;
    for (int i = 0; i < s.m; ++i) {
      const double lr = s.lf[s.r[i]];
      coff[i].push_back(0);
      compositions(s, s.cm, kc, s.suffix_r[i],
        [&](long long sk, double) {
          kc.decode(sk, vbuf);
          if (!is_canon(vbuf)) return;
          sidx[i].emplace(sk, (int)cst[i].size());
          cst[i].push_back(sk);
          compositions(s, vbuf, kc, s.r[i],
            [&](long long kd, double neg_lfac) {
              arc_ck[i].push_back(canon(sk - kd));
              arc_st[i].push_back(lr + neg_lfac);
            });
          coff[i].push_back((int)arc_ck[i].size());
          if ((long long)arc_ck[i].size() > MAX_ARCS)
            Rcpp::stop("arc storage limit exceeded");
        });
    }
  }

  // ---- Exact completion extrema, backward (Mehta-Patel 1980's
  // induction, vector states, dense keys, canonical states
  // only). ----
  std::vector<std::vector<double>> vmax(s.m + 1), vmin(s.m + 1);
  vmax[s.m].assign(1, 0.0);   // terminal state has key 0
  vmin[s.m].assign(1, 0.0);
  for (int i = s.m - 1; i >= 0; --i) {
    vmax[i].assign(KS, -INFINITY);
    vmin[i].assign(KS, INFINITY);
    const double lr = s.lf[s.r[i]];
    const std::vector<double>& nmx = vmax[i + 1];
    const std::vector<double>& nmn = vmin[i + 1];
    if (has_ties) {
      for (size_t u = 0; u < cst[i].size(); ++u) {
        double mx = -INFINITY, mn = INFINITY;
        for (int t = coff[i][u]; t < coff[i][u + 1]; ++t) {
          const long long ck = arc_ck[i][t];
          const double step = arc_st[i][t];
          const double cmx = nmx[ck];
          if (std::isfinite(cmx)) {
            if (step + cmx > mx) mx = step + cmx;
            const double cmn = nmn[ck];
            if (step + cmn < mn) mn = step + cmn;
          }
        }
        vmax[i][cst[i][u]] = mx;
        vmin[i][cst[i][u]] = mn;
      }
    } else {
      compositions(s, s.cm, kc, s.suffix_r[i],
        [&](long long sk, double) {
          kc.decode(sk, vbuf);
          double mx = -INFINITY, mn = INFINITY;
          compositions(s, vbuf, kc, s.r[i],
            [&](long long kd, double neg_lfac) {
              const long long ck = sk - kd;
              const double step = lr + neg_lfac;
              const double cmx = nmx[ck];
              if (std::isfinite(cmx)) {
                if (step + cmx > mx) mx = step + cmx;
                const double cmn = nmn[ck];
                if (step + cmn < mn) mn = step + cmn;
              }
            });
          vmax[i][sk] = mx;
          vmin[i][sk] = mn;
        });
    }
  }

  // ---- Forward sweep with hashed class maps per canonical
  // state. ----
  std::unordered_map<long long, ClsMap> cur, nxt;
  cls_add(cur[kc.encode(s.cm)], 0.0, 1.0);
  s.pval = 1.0;

  for (int i = 0; i < s.m && !cur.empty(); ++i) {
    const double lr = s.lf[s.r[i]];
    const std::vector<double>& nmx = vmax[i + 1];
    const std::vector<double>& nmn = vmin[i + 1];
    for (auto& kvp : cur) {
      const long long sk = kvp.first;
      ClsMap& cl = kvp.second;
      if (cl.empty()) continue;
      kc.decode(sk, vbuf);
      const double vmx = vmax[i][sk];
      const double vmn = vmin[i][sk];
      double lbulk = s.lf[s.suffix_r[i]];
      for (int j = 0; j < s.c; ++j) lbulk -= s.lf[vbuf[j]];

      for (const auto& ep : cl) {
        const Cls& e = ep.second;
        const double base = s.log_D + e.q;
        if (!std::isfinite(vmx) || base + vmx <= s.log_thresh)
          continue;
        if (base + vmn > s.log_thresh) {
          s.pval -= e.w * std::exp(base + lbulk);
          continue;
        }
        // Straddling: expand row i, classifying each child at
        // creation (dropped and bulk children never enter the
        // class lists). With ties, scan the stored arcs; without,
        // re-enumerate.
        auto handle_child = [&](long long ck, double step) {
          const double q2 = e.q + step;
          const double b2 = s.log_D + q2;
          const double cmx = nmx[ck];
          if (!std::isfinite(cmx) || b2 + cmx <= s.log_thresh)
            return;
          if (b2 + nmn[ck] > s.log_thresh) {
            // Bulk mass: multinomial Vandermonde on the child.
            double lb = s.lf[s.suffix_r[i + 1]];
            long long tt = ck;
            for (int j = 0; j < s.c; ++j) {
              const int vj = (int)(tt / kc.W[j]);
              tt -= (long long)vj * kc.W[j];
              lb -= s.lf[vj];
            }
            s.pval -= e.w * std::exp(b2 + lb);
            return;
          }
          cls_add(nxt[ck], q2, e.w);
        };
        if (has_ties) {
          const int u = sidx[i].at(sk);
          for (int t = coff[i][u]; t < coff[i][u + 1]; ++t)
            handle_child(arc_ck[i][t], arc_st[i][t]);
        } else {
          compositions(s, vbuf, kc, s.r[i],
            [&](long long kd, double neg_lfac) {
              handle_child(sk - kd, lr + neg_lfac);
            });
        }
      }
    }
    cur.swap(nxt);
    nxt.clear();
  }

  if (s.pval < 0.0) s.pval = 0.0;
  if (s.pval > 1.0) s.pval = 1.0;
  return s.pval;
}
