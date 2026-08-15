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
// workspace parameter. Permutation-collapsing of tied-margin
// states and the child-creation classification of the m x 2
// kernel are both implemented below. A last-two-rows closed form
// was built and measured neutral to slightly negative, because
// classification at creation already makes the final stages
// free; it is recorded in the flagship paper and not retained.
//
// READER'S MAP, keyed to the r x c paper (analysis/report/02-rxc,
// "The State-Aggregation Routine" and its subsections). The
// paper's notation: an m x c table Y = (y_ij) with row margins
// r_i, column margins c_j, total n; log P(Y) = log_D
// + sum_i [log r_i! - sum_j log y_ij!] with
// log_D = log[prod_j c_j! / n!]; the complement region
// U = { Y : P(Y) > P_obs (1 + tau) }; and the ledger
// p = 1 - mass(U).
//
//   paper object / subsection          code location
//   ------------------------------     ------------------------
//   state = residual column margins    mixed-radix key sk; the
//     ("Representation")               decoded vector is vbuf
//   transposition, row ordering,       setup block of
//     permutation collapsing           rxc_net_sweep_cpp();
//     ("Representation")               canon(), is_canon()
//   transitions (one row placement)    compositions() /
//     ("Representation")               comps_rec()
//   exact completion extrema           backward pass filling
//     ("Certificates from backward     vmax, vmin
//     induction")
//   stage-level mass identity,         lmult, filled in the
//     Mehta-Patel eq. (2.4)            backward pass; read by
//     ("Certificates...")              the bulk move
//   class (q, w)                       struct Cls; container
//     ("Classes and their merging")    FlatCls
//   drop / bulk / propagate            the per-arc run
//     ("The sweep")                    classification in the
//                                      forward sweep
//   the ledger p = 1 - mass(U)         s.pval, decremented by
//     (Section "The Convex             every bulk disposal
//     Complement Region")

#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <cmath>

using namespace Rcpp;

namespace {

// The largest supported state dimension: sufcap in
// compositions() is sized to it, and the entry guard enforces it.
constexpr int kMaxC = 16;

// A class (paper: "Classes and their merging"): all path prefixes
// reaching one state with equal accumulated log mass, carried as
// a single object. The full table probability of any completion
// through this class is log_D + q + (log mass of the remaining
// rows), which is why every classification below is a comparison
// of q against a threshold.
struct Cls {
  double q;   // accumulated log path weight of the placed rows:
              //   sum_i [ log r_i! - sum_j log x_ij! ]
  double w;   // linear multiplicity relative to q: the number of
              // prefixes merged in, each weighted exp(q' - q)
};

// Global quantities of the instance, in the paper's notation.
struct FusedState {
  int m, c, n;               // rows, columns, grand total (after
                             // the transposition of the setup)
  std::vector<int> r;        // row margins r_i, processing order
  std::vector<int> cm;       // column margins c_j, descending
  std::vector<int> suffix_r; // suffix_r[i] = sum_{g >= i} r_g,
                             // the total of the unplaced rows
  double log_D;              // log[prod_j c_j! / n!], the
                             // row-free factor of log P(Y)
  double log_thresh;         // log[P_obs (1 + tau)], the
                             // boundary of U
  std::vector<double> lf;    // lf[v] = log v!
  double pval;               // the ledger: starts at 1, ends at
                             // 1 - mass(U)
};

// Mixed-radix key machinery for residual vectors (paper:
// "Representation"). A state, the vector of residual column
// margins after the placed rows, is addressed by the integer
// sum_j v_j W[j] with radix c_j + 1 and weight W[j] for
// coordinate j; the root state is (c_1, ..., c_c) and the
// terminal state is the zero vector, key 0. A child key is the
// parent key minus sum_j x_j W[j], computed incrementally during
// enumeration, so no child vector is ever materialized on the
// hot path.
struct KeyCoder {
  std::vector<long long> W;
  void init(const std::vector<int>& cm) {
    const int c = (int)cm.size();
    W.assign(c, 1);
    for (int j = c - 2; j >= 0; --j)
      W[j] = W[j + 1] * (cm[j + 1] + 1);
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
// 0 <= x_j <= v[j]: one composition is one placement of the
// current row within the residual margins, that is, one
// transition (arc) of the state graph, and the bounds on each
// x_j are the feasibility window of the paper's "Paths and
// partial tables". Calls f(key_delta, neg_lfac) where
// key_delta = sum_j x_j W[j] and neg_lfac = -sum_j log x_j!.
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
  int sufcap[kMaxC];
  sufcap[s.c] = 0;
  for (int j = s.c - 1; j >= 0; --j)
    sufcap[j] = sufcap[j + 1] + v[j];
  comps_rec(s, v.data(), kc.W.data(), sufcap, 0, total, 0.0,
            0LL, f);
}

inline long long qkey(double q) {
  return llround(q * 1e12);
}

// Class container per state: an open-addressing table on the
// quantized log mass (grid 1e-12, so only floating-point ties
// merge and exactness is preserved). Hashing is the device FEXACT
// always used here; an earlier draft used std::unordered_map, and
// sample profiling on the classical suite put 43 percent of run
// time in its node-based emplace and another 22 percent in the
// exp() of the merge arithmetic (5.3 million merges on the
// largest problem, 27 paths per surviving class). The flat table
// removes the per-node allocation and is reused across stages by
// epoch stamping, so clearing is O(1); and colliding keys agree
// to 1e-12 by construction, so exp(q - e.q) is replaced by its
// first-order expansion, with relative error below 1e-24, ten
// orders beneath the kernel's measured 1e-14 agreement with
// FEXACT. A sorted-run merge was prototyped in place of the
// table and measured 4.7x slower on the largest problem: hashing
// deduplicates at insertion, so the 27x tie multiplicity never
// reaches memory, while run merging materializes every appended
// class before folding; deferred deduplication is the wrong
// trade at high multiplicity.
struct FlatCls {
  struct Slot { long long key; int idx; unsigned stamp; };
  std::vector<Slot> tab;
  std::vector<Cls> items;
  unsigned epoch = 0;
  size_t mask = 0;

  void reset() {
    items.clear();
    if (tab.empty()) {
      tab.assign(32, Slot{0, 0, 0});
      mask = 31;
      epoch = 1;
    } else if (++epoch == 0) {  // stamp wrap: wipe once, restart
      std::fill(tab.begin(), tab.end(), Slot{0, 0, 0});
      epoch = 1;
    }
  }
  static inline size_t hashk(long long k) {
    return (size_t)((unsigned long long)k *
                    0x9E3779B97F4A7C15ULL >> 32);
  }
  void grow() {
    tab.assign(tab.size() * 2, Slot{0, 0, 0});
    mask = tab.size() - 1;
    for (int t = 0; t < (int)items.size(); ++t) {
      const long long k = qkey(items[t].q);
      size_t h = hashk(k) & mask;
      while (tab[h].stamp == epoch) h = (h + 1) & mask;
      tab[h] = Slot{k, t, epoch};
    }
  }
  inline void add(double q, double w) {
    const long long k = qkey(q);
    size_t h = hashk(k) & mask;
    while (true) {
      Slot& sl = tab[h];
      if (sl.stamp != epoch) {
        sl = Slot{k, (int)items.size(), epoch};
        items.push_back(Cls{q, w});
        if (items.size() * 2 >= tab.size()) grow();
        return;
      }
      if (sl.key == k) {
        Cls& e = items[sl.idx];
        // |q - e.q| < 1e-12: equal keys. First-order expansion
        // of exp(q - e.q); the quadratic term is below 1e-24.
        e.w += w * (1.0 + (q - e.q));
        return;
      }
      h = (h + 1) & mask;
    }
  }
};

}  // anonymous namespace

// Entry point. Computes the two-sided p-value of Fisher's exact
// test under the minimum-likelihood ordering: the total
// conditional probability, given both margins, of the tables no
// more probable than the observed one, with FEXACT's tolerance
// convention tau = 3.45254e-7 (paper: "The Path Representation",
// and docs/tolerance_analysis_whitepaper.tex). Computed as
// 1 - mass(U) per the complement-region ledger.
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
  if (s.c >= kMaxC)
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

  // Permutation collapsing of tied-margin states (paper:
  // "Representation"; Mehta-Patel 1983, implementation device
  // 2). Tied-column groups are contiguous runs because the
  // columns were sorted descending. A state is canonical when
  // its residuals are non-increasing within every tied group;
  // states differing only by a permutation within a tied group
  // have isomorphic completion structure, so all are folded onto
  // the canonical representative by canon(). For tables with no
  // ties canon is the identity and costs one flag test.
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
  // directly per stage. Arc handling is hybrid, thresholded by
  // measurement on the size of the state key space. Small tied
  // tables re-enumerate arcs on the fly with the child
  // canonicalized per arc: the CSR store they previously built
  // was measured at 16 to 64 percent of their total run time
  // against a backward-pass consumer of 2 to 5 percent, and
  // removing it gained up to 2x on the smallest tables. Large
  // tied tables keep the CSR store: at mp5 scale the per-arc
  // canonicalization of the backward pass costs more than the
  // build (measured 0.84x without the store). Tie-free tables
  // always enumerate on the fly, canon being the identity.
  const bool use_csr = has_ties && KS > (1LL << 13);
  std::vector<std::vector<long long>> csr_state(s.m);
  std::vector<std::vector<int>> csr_off(s.m);
  std::vector<std::vector<long long>> csr_child(s.m);
  std::vector<std::vector<double>> csr_step(s.m);
  std::vector<std::unordered_map<long long, int>> csr_idx(s.m);
  std::vector<int> vbuf(s.c);
  if (use_csr) {
    const long long MAX_ARCS = 1LL << 26;
    for (int i = 0; i < s.m; ++i) {
      const double lr = s.lf[s.r[i]];
      csr_off[i].push_back(0);
      compositions(s, s.cm, kc, s.suffix_r[i],
        [&](long long sk, double) {
          kc.decode(sk, vbuf);
          if (!is_canon(vbuf)) return;
          csr_idx[i].emplace(sk, (int)csr_state[i].size());
          csr_state[i].push_back(sk);
          compositions(s, vbuf, kc, s.r[i],
            [&](long long kd, double neg_lfac) {
              csr_child[i].push_back(canon(sk - kd));
              csr_step[i].push_back(lr + neg_lfac);
            });
          csr_off[i].push_back((int)csr_child[i].size());
          if ((long long)csr_child[i].size() > MAX_ARCS)
            Rcpp::stop("arc storage limit exceeded");
        });
    }
  }

  // ---- Exact completion extrema, backward (paper:
  // "Certificates from backward induction"; Mehta-Patel 1980's
  // device, vector states, dense keys, canonical states only).
  // ----
  // After this pass, vmax[i][sk] and vmin[i][sk] are the TRUE
  // maximum and minimum, over all completions of state sk at
  // stage i, of the remaining rows' log mass
  // sum_{g >= i} [log r_g! - sum_j log y_gj!]. A class with
  // accumulated mass q therefore has completions spanning
  // exactly [log_D + q + vmin, log_D + q + vmax] on the log
  // scale, and every classification in the sweep compares those
  // endpoints against log_thresh. No relaxation is involved,
  // which is why the water-filling and concentration bounds of
  // the paper's Section 3 play no role in this routine.
  //
  // lmult[i][sk] is the log multinomial mass of state sk at stage
  // i, lf[suffix_r[i]] - sum_j lf[v_j]: the total path weight of
  // the state's completions (Mehta-Patel 1983, eq. 2.4; the
  // paper's stage-level mass identity). Filled here because the
  // backward pass visits every state anyway; the forward sweep's
  // bulk move then reads it in one lookup instead of decoding
  // the child key with c divisions per bulk arc.
  std::vector<std::vector<double>> vmax(s.m + 1), vmin(s.m + 1);
  std::vector<std::vector<double>> lmult(s.m + 1);
  vmax[s.m].assign(1, 0.0);   // terminal state has key 0
  vmin[s.m].assign(1, 0.0);
  lmult[s.m].assign(1, 0.0);
  for (int i = s.m - 1; i >= 0; --i) {
    vmax[i].assign(KS, -INFINITY);
    vmin[i].assign(KS, INFINITY);
    lmult[i].assign(KS, 0.0);
    const double lr = s.lf[s.r[i]];
    const double lsuf = s.lf[s.suffix_r[i]];
    const std::vector<double>& nxt_max = vmax[i + 1];
    const std::vector<double>& nxt_min = vmin[i + 1];
    if (use_csr) {
      for (size_t u = 0; u < csr_state[i].size(); ++u) {
        double mx = -INFINITY, mn = INFINITY;
        for (int t = csr_off[i][u]; t < csr_off[i][u + 1]; ++t) {
          const long long ck = csr_child[i][t];
          const double step = csr_step[i][t];
          const double cmx = nxt_max[ck];
          if (std::isfinite(cmx)) {
            if (step + cmx > mx) mx = step + cmx;
            const double cmn = nxt_min[ck];
            if (step + cmn < mn) mn = step + cmn;
          }
        }
        const long long sk = csr_state[i][u];
        vmax[i][sk] = mx;
        vmin[i][sk] = mn;
        kc.decode(sk, vbuf);
        double lm = lsuf;
        for (int j = 0; j < s.c; ++j) lm -= s.lf[vbuf[j]];
        lmult[i][sk] = lm;
      }
    } else {
      // The has_ties split is hoisted out of the enumeration so
      // the tie-free inner loop carries no canon call at all and
      // compiles as it did before canonicalization moved
      // on the fly.
      auto scan = [&](auto&& child_key) {
        compositions(s, s.cm, kc, s.suffix_r[i],
          [&](long long sk, double) {
            kc.decode(sk, vbuf);
            if (!is_canon(vbuf)) return;
            double lm = lsuf;
            for (int j = 0; j < s.c; ++j) lm -= s.lf[vbuf[j]];
            lmult[i][sk] = lm;
            double mx = -INFINITY, mn = INFINITY;
            compositions(s, vbuf, kc, s.r[i],
              [&](long long kd, double neg_lfac) {
                const long long ck = child_key(sk - kd);
                const double step = lr + neg_lfac;
                const double cmx = nxt_max[ck];
                if (std::isfinite(cmx)) {
                  if (step + cmx > mx) mx = step + cmx;
                  const double cmn = nxt_min[ck];
                  if (step + cmn < mn) mn = step + cmn;
                }
              });
            vmax[i][sk] = mx;
            vmin[i][sk] = mn;
          });
      };
      if (has_ties) scan([&](long long k) { return canon(k); });
      else scan([](long long k) { return k; });
    }
  }

  // ---- Forward sweep with pooled flat class tables per
  // canonical state (paper: "The sweep"). ----
  //
  // The three-way classification of the paper, per class-arc
  // pair: a child whose completions all satisfy
  // P <= P_obs(1 + tau) lies wholly in the significance region
  // and is DROPPED, contributing nothing to mass(U); a child
  // whose completions all exceed the threshold lies wholly in U
  // and is BULK-disposed, its entire mass moved to the ledger in
  // one step via the stage-level identity; only children that
  // straddle PROPAGATE to the next stage. The ledger is
  // s.pval = 1 - mass(U), so each bulk disposal is a
  // subtraction from pval.
  //
  // Loop order is arcs-outside-classes. Arcs are a property of the
  // state, not of the class, so enumerating them inside the class
  // loop walked a state's arc list once per class. Hoisting the
  // walk out also hoists the child extrema lookups, the child
  // multinomial lookup, and the destination-table probe, each of
  // which is then paid once per arc rather than once per class-arc
  // pair.
  //
  // The inversion pays a second time through monotonicity. Both
  // classification tests compare log_D + q + step + (extremum)
  // against log_thresh, so at a fixed arc each is a threshold on q
  // alone; with the state's classes held sorted by q the survivors
  // of an arc are a contiguous run, located by two binary searches.
  // Dropped classes then cost nothing, and the bulk run collapses
  // to a single exp() against a precomputed suffix sum. Classes are
  // accumulated by hashing (FlatCls) and sorted once per state per
  // stage, in place, so the O(K log K) is paid once rather than
  // per arc: this is why the sort does not reintroduce the cost
  // that motivated hashed class merging. The tables live in pools
  // indexed by a small state map and are recycled across stages by
  // epoch stamping, so steady state allocates nothing.
  //
  // There is no state-level classification. A class enters a
  // state's table only if it straddled at creation, and the
  // creation test reads the same vmax/vmin entries a state-level
  // re-test would read, so the re-test can never fire. (A class
  // whose merged representative q drifted within the 1e-12 bin
  // could flip a comparison landing inside that bin, but both
  // outcomes are then correct to within the 3.45e-7 tolerance,
  // five orders coarser.) The root is the one state never
  // classified; were the whole reference set disposable, its
  // children would all drop or bulk in the first stage, reaching
  // the same answer one stage later.
  std::unordered_map<long long, int> cur_idx, nxt_idx;
  std::vector<FlatCls> cur_pool, nxt_pool;
  int nxt_used = 0;
  cur_pool.emplace_back();
  cur_pool[0].reset();
  cur_pool[0].add(0.0, 1.0);
  cur_idx.emplace(kc.encode(s.cm), 0);
  s.pval = 1.0;

  std::vector<double> tsum;       // class suffix sums, per state
  std::vector<long long> arc_child;  // arc child keys (non-CSR)
  std::vector<double> arc_step;      // arc steps (non-CSR)
  // upper_bound over classes by q: first class with q > value.
  auto q_gt = [](double v, const Cls& e) { return v < e.q; };

  for (int i = 0; i < s.m && !cur_idx.empty(); ++i) {
    const double lr = s.lf[s.r[i]];
    const std::vector<double>& nxt_max = vmax[i + 1];
    const std::vector<double>& nxt_min = vmin[i + 1];
    const std::vector<double>& nxt_lmult = lmult[i + 1];
    for (auto& kvp : cur_idx) {
      const long long sk = kvp.first;
      std::vector<Cls>& cv = cur_pool[kvp.second].items;
      const int K = (int)cv.size();
      std::sort(cv.begin(), cv.end(),
                [](const Cls& a, const Cls& b) { return a.q < b.q; });

      // Suffix sums anchored at the largest q, so every term is
      // exp() of a non-positive argument and an arc's bulk slice
      // [hi, K) is read off as one entry: always a sum, never a
      // difference of two large sums.
      const double qref = cv[K - 1].q;
      tsum.resize(K + 1);
      tsum[K] = 0.0;
      for (int t = K - 1; t >= 0; --t)
        tsum[t] = tsum[t + 1] + cv[t].w * std::exp(cv[t].q - qref);

      // Arcs of this state, enumerated once. Above the CSR
      // threshold, read the prebuilt run in place; below it,
      // materialize into reusable buffers with the child
      // canonicalized per arc (identity for tie-free tables).
      const long long* ck_p;
      const double* st_p;
      int na;
      if (use_csr) {
        const int u = csr_idx[i].at(sk);
        ck_p = csr_child[i].data() + csr_off[i][u];
        st_p = csr_step[i].data() + csr_off[i][u];
        na = csr_off[i][u + 1] - csr_off[i][u];
      } else {
        kc.decode(sk, vbuf);
        arc_child.clear();
        arc_step.clear();
        if (has_ties) {
          compositions(s, vbuf, kc, s.r[i],
            [&](long long kd, double neg_lfac) {
              arc_child.push_back(canon(sk - kd));
              arc_step.push_back(lr + neg_lfac);
            });
        } else {
          compositions(s, vbuf, kc, s.r[i],
            [&](long long kd, double neg_lfac) {
              arc_child.push_back(sk - kd);
              arc_step.push_back(lr + neg_lfac);
            });
        }
        ck_p = arc_child.data();
        st_p = arc_step.data();
        na = (int)arc_child.size();
      }

      for (int a = 0; a < na; ++a) {
        const long long ck = ck_p[a];
        const double cmx = nxt_max[ck];
        if (!std::isfinite(cmx)) continue;
        const double step = st_p[a];
        const double base = s.log_thresh - s.log_D - step;
        const int lo = (int)(std::upper_bound(
            cv.begin(), cv.end(), base - cmx, q_gt)
            - cv.begin());
        const int hi = (int)(std::upper_bound(
            cv.begin() + lo, cv.end(), base - nxt_min[ck], q_gt)
            - cv.begin());
        if (hi < K)
          // Bulk slice: classes cv[hi..K) lie wholly inside U
          // through this arc, so their entire mass moves to the
          // ledger at once; nxt_lmult[ck] is the child's total
          // completion weight (eq. 2.4) and tsum[hi] the classes'
          // summed weight, one exp() for the whole run.
          s.pval -= std::exp(s.log_D + step + nxt_lmult[ck] + qref) *
                    tsum[hi];
        if (lo < hi) {
          auto r = nxt_idx.emplace(ck, nxt_used);
          if (r.second) {
            if (nxt_used == (int)nxt_pool.size())
              nxt_pool.emplace_back();
            nxt_pool[nxt_used].reset();
            ++nxt_used;
          }
          FlatCls& dst = nxt_pool[r.first->second];
          for (int t = lo; t < hi; ++t)
            dst.add(cv[t].q + step, cv[t].w);
        }
      }
    }
    cur_idx.swap(nxt_idx);
    cur_pool.swap(nxt_pool);
    nxt_idx.clear();
    nxt_used = 0;
  }

  if (s.pval < 0.0) s.pval = 0.0;
  if (s.pval > 1.0) s.pval = 1.0;
  return s.pval;
}
