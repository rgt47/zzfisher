# zzfisher 0.8.0

* New `ss_fxpower_fast()`: a faster variant of `ss_fxpower()`, using
  the same doubling (exponential) search + bisection strategy already
  used by `fisher_power_fast()`'s `n1`-solving path, instead of
  `ss_fxpower()`'s linear scan incrementing `n2` one at a time. The
  underlying power evaluator (`fxpower()`) is completely unchanged --
  this is purely a faster search strategy layered on top of it, and
  the fix is independent of (and now complete alongside) the
  `fisher_power_fast()` fix for the pure-R side: the same
  linear-scan-to-bisection idea now applies to both the `enum` and
  `mode` power algorithms' sample-size solvers.
* Verified bitwise identical to `ss_fxpower()` across a grid of
  targets, including both allocation directions and `k != 1`.
  Benchmarked: 25x faster at a target n ~ 500 growing to **162x**
  faster at a target n ~ 8900 (9.7s down to 0.06s) -- a realistic
  small-effect-size design in the low-event-rate trials this package
  family targets, not a contrived worst case.

# zzfisher 0.7.0

* `fxpower()` gains a two-sided option (`alternative = c("less",
  "two.sided")`, default unchanged at `"less"`): the compiled
  mode-centered kernel (`src/fxpower_rcpp.cpp`) previously had no
  two-sided rejection rule at all, by construction rather than by
  omission. The new `cv_two_sided()` locates the minimum-likelihood
  rejection region's two boundaries via binary search on the null
  hypergeometric's monotonic tails (the same idea used in
  `fisher_power_fast()`'s two-sided fix, ported to C++ using
  `std::upper_bound`), reusing the existing `fsum_unequal()`
  adjacency-recurrence summation and the mode-finding machinery
  already in place for the one-sided path. Verified against
  `fisher_power_fast(alternative = "two.sided")` to ~1e-12 (exact) and
  within the reported eps-trimming error bound (checked, not merely
  assumed) at every configuration tested, up to n = 1000.
* Fixed a bug found and caught by that verification before release:
  the first implementation of the two-sided eps-trimming walked each
  rejection block from its own outer boundary rather than from the
  alternative distribution's own mode (`xfm`) when the mode falls
  inside that block. Since `binom_joint()` is itself unimodal in `x`,
  starting from the wrong end could make the summed values increase
  rather than decrease along the walk, tripping the `f < flim`
  early-stop immediately and silently discarding a whole block's
  contribution -- collapsing power to 0 at large n where trimming
  inevitably engages (n = 5000 triggered it outright; n <= 2000
  happened not to). Fixed by mirroring the one-sided path's existing
  `xfm`-aware block-splitting logic for each side. Regression-tested
  in `inst/tinytest/test_two_sided.R`.

# zzfisher 0.6.1

* `fisher_power_fast()`'s `alternative = "two.sided"` path no longer
  scans the whole hypergeometric support for every candidate table
  (an O(range) `vapply()`/`sum()` per candidate x, i.e. O(range^2) per
  m). Since the hypergeometric density is unimodal, the rejection set
  for any threshold is always a prefix of the increasing left tail
  plus a suffix of the decreasing right tail; `findInterval()`
  (binary search) locates both boundaries for every candidate at once
  against precomputed prefix/suffix sums, for O(range log range) per
  m. Verified against `fisher_power()`'s original two-sided
  implementation across the migration's test grid, including
  boundary-peak edge cases (p1/p2 forcing the mode to either end of
  the support): agreement is exact (0 discrepancy) at every case
  checked, not merely within tolerance. Benchmarked: 4.3x faster at
  n = 200 growing to 31.5x at n = 2000, and now tractable well past
  the n <= 200 practical limit of the original scan (n = 2000
  completes in ~1.1s where the original could not finish inside a 15s
  budget past n = 200).

# zzfisher 0.6.0

* New `fisher_power_fast()`: a faster, numerically identical variant
  of `fisher_power()`, added as a competing implementation in
  `inst/benchmarks/benchmark_power_comparison.R`. Two independent
  fixes, both exact (no eps-trimming): (1) removes a `dhyper()` call
  that was computed unconditionally every iteration of the one-sided
  paths despite only the two-sided rejection rule needing it; (2)
  replaces the linear scan used to solve for `n1` with an exponential
  (doubling) search followed by bisection, which never evaluates power
  near `n1_max` unless the search genuinely needs to (a naive
  bisection against `n1_max` would pay for one expensive evaluation at
  the cap on its very first step). Verified bitwise identical to
  `fisher_power()` at every tested configuration; on one benchmark
  case (p1 = 0.05, p2 = 0.10, power = 0.9, alternative = "less",
  landing at n1 = 503) `fisher_power_fast()` took 0.76s against
  `fisher_power()`'s 24.2s, a ~32x speedup from reducing the number of
  power evaluations during n1-solving alone.

# zzfisher 0.5.0

* New `fisher_power()`: exact power and sample size for the r x 2
  Fisher's exact test, following Conlon & Thomas (1993, Algorithm
  AS 280). A direct re-derivation from the power function's
  mathematical definition (exhaustive enumeration of the null
  rejection region against the product-binomial alternative), not a
  port of the original Fortran/NAG source. Mirrors
  `stats::power.t.test()`'s calling convention: leave `n1` or `power`
  `NULL` and it solves for the missing one, returning a classed
  `power.htest` object. Cross-checked against
  `Exact::power.exact.test(method = "fisher")` to 6 decimal places.
* New `fxpower()`/`fxpower_r()`/`ss_fxpower()`/`power_table()`/
  `ss_table()`: a second, independently developed power-function
  implementation migrated from the research workspace
  `rgt47/fisherpowerunequaln`. `fxpower()` is a compiled C++
  mode-finding algorithm generalizing Algorithm AS 280 to unequal
  allocation, with an `eps` trimming parameter trading a quantified
  error bound for one to two orders of magnitude of speed; `eps = 0`
  gives the exact one-sided power. Verified bitwise identical to
  `fisher_power(..., alternative = "less")` across the migration's
  test grid.

# zzfisher 0.4.1

* Picks up sweep-kernel work completed in the research workspace
  after 0.4.0 was cut. `rx2_net_sweep.cpp` tabulates the row term
  `h_k(y)` once per row instead of recomputing it inside the budget
  loop. `rxc_net_sweep.cpp` gains permutation collapsing of
  tied-margin states and child-creation classification, and carries
  a reader's map keyed to the sections of the r x c paper. No change
  to any returned p-value: verified bitwise identical to the
  research build across 57 tables and all ten dispatchers.
* Restores `importFrom(stats, dhyper)`, dropped from NAMESPACE when
  `wrappers.R` was replaced by the dispatchers in 0.2.0. The pure-R
  `.rx2_tree_dp` kernel calls `dhyper()`, so the package had been
  relying on stats being attached rather than importing it.
* Drops the obsolete `CXX_STD = CXX11` from `src/Makevars`; R
  already ignores it in favor of the default standard.
* Adds drift detection against the research workspace. The 21 files
  in `tools/sync-manifest.txt` are copied verbatim from
  rgt47/fisherexacttestrx2 and must not be edited here;
  `tools/check-sync.sh` verifies them in CI on every push and pull
  request to `main`, and `tools/sync-provenance.txt` records the
  upstream commit this build came from.

# zzfisher 0.4.0

* Fused state-aggregation kernel for r x c tables
  (`net_sweep_cpp`'s rxc arm now routes to it): Mehta-Patel
  network states (column-residual vectors, table transposed so
  the state dimension is the smaller side), EXACT completion
  extrema by backward induction with dense mixed-radix keys (the
  1980 device, vector states), hashed class merge, three-case
  classification with child-creation filtering, closed-form bulk
  disposal via the multinomial Vandermonde, and permutation
  collapsing of tied-margin states with hybrid arc handling (CSR
  storage with ties, on-the-fly enumeration without). Verified
  to 7.6e-14 over the r x c oracle gate; completes all six
  Mehta-Patel (1983) Table 1 problems, two of which are
  infeasible for the tree kernel.
* New test and kernel `net_ci_cpp`: stratified
  conditional-independence exact test for r x c x k arrays (rows
  independent of columns given the layer, log-linear [AC][BC];
  the stratified Fisher's exact test). Per-stratum
  log-probability distributions convolved across strata with
  exact suffix extrema and O(1) bulk moves; k = 1 reduces to the
  ordinary r x c Fisher test, checked against fisher.test in the
  test suite alongside brute-force 2 x 2 x k oracles.

* Measurement-validity note: all pre-0.4.0 benchmark claims were
  taken on pkgload::load_all() builds, which compile with debug
  flags (-O0 -UNDEBUG) and slowed these kernels up to 16x against
  R's optimized libraries. Every figure below is from an
  R CMD INSTALL (-O2) build.
* Sweep constant-factor pass: children of a straddling class are
  classified at creation against per-stage threshold arrays, so
  dropped and bulk-subtracted children never enter the stage
  buffers and no stage-entry classification pass exists. On the
  36-cell wrapper-stripped grid the sweep beats raw FEXACT in all
  36 cells: median cell ratio 9.0, from 18x on small tables down
  to 1.6-2.0x at the heaviest (m = 8, rho = 10).
* r x c kernel: the balanced-split skip is upgraded to capped
  water-filling for the current column (the Joe-style
  exact-relaxation bound, respecting row residuals) with a
  concavity-based sibling cutoff (the r x c analog of the m x 2
  cascade). Verified to 7.5e-13 against fisher.test over 170
  random tables; 3.7x faster overall on a fixed seeded set, with
  every cell gaining (2.9x to 10.2x).
* r x c x k kernel: determine-as-you-go restructure (a layer's
  dependent cells are computed, checked, and charged to the AB
  margins at the layer boundary; running prefix statistic over a
  log-factorial table replaces the per-leaf O(mck) lgamma sweep).
  Bit-identical to the predecessor over 90 random arrays and
  performance-neutral (0.96-1.01): the ratio-form p-value admits
  no valid skip, and placement bounding already keeps prefixes
  near-feasible. Documented as a measured negative result.
* New kernel `net_sweep_cpp`: a state-sweep (network) algorithm for
  m x 2 tables, per Sections 3B and 3B.1 of the project white paper
  `c2_separability_whitepaper_2026-08-04`. The sweep proceeds stage
  by stage over (stage, spent budget) states; paths reaching a
  state are aggregated into prefix log-mass classes (the
  Clarkson-Fan-Joe past device, with dense arrays in place of
  hashing); classes are classified against the exact suffix
  extrema of the c = 2 separability DP, so a class entirely in the
  complement region is subtracted in O(1) via a closed-form
  Vandermonde suffix mass, a class entirely in the significance
  region is dropped, and only straddling classes propagate. At the
  last two rows, leaf masses and prefix sums are shared per state
  and each class contributes a single prefix-sum difference, with
  nested qualifying intervals extracted by expanding two pointers.
  This is the Mehta-Patel-Joe-Clarkson-Requena algorithm
  specialized to two columns. Unlike FEXACT it has no fixed
  workspace: storage grows with the live class count.
* Verified against `fisher.test()` over 402 random tables (m = 2
  to 10, including tied margins up to the m = 8, rho = 10 regime,
  degenerate tables, and the large-m regime that forces FEXACT
  workspace growth): max absolute discrepancy 8.0e-14.

# zzfisher 0.3.0

* Redesigned the m x 2 tree kernel (`tree_memo_c`) around the
  separability of the objective when the table has two columns.
  The maximum and minimum subtree bounds are now the exact extrema,
  computed once per table by a backward dynamic program over
  (level, residual budget), replacing the branch-and-bound
  `find_max` search, the `joe_min_impl` corner recursion, and both
  memoization hash tables. Every bound query is an array read, and
  all threshold comparisons are performed in log space.
* Added a proven two-sided sibling cutoff. Concavity of the value
  function makes the children requiring attention a contiguous
  interval, so the scan stops at the first skippable child in each
  direction; the same argument truncates the leaf-level walk. This
  is the proved m x 2 form of the cascade that `tree_s4_c` still
  assumes without proof.
* Added an acceptance-region skip to the r x c kernel based on a
  balanced-split relaxation, which is valid for any number of
  columns. Its practical value is not yet quantified.
* No user-visible API change. Verified against a complete-
  enumeration oracle (max discrepancy 8.2e-15 over 96 m x 2 tables,
  2.3e-14 over 34 r x c tables) and against the previous kernel
  (max drift 3.5e-15, summation order only); 204 tests pass.
* Mathematics, with proofs, in the project white paper
  `c2_separability_whitepaper_2026-08-04`.

# zzfisher 0.2.1

* Fixed an undefined-behavior bug (ODR violation): two translation
  units defined incompatible global-scope `CacheEntry` structs used
  as `std::unordered_map` value types, which could corrupt the
  memoization cache of `tree_memo_cpp()` and crash R on some tables.
  All internal kernel structs and helpers now live in anonymous
  namespaces.
* Aligned internal identifiers with the manuscripts' Concrete
  Mathematics notation: hypergeometric helper parameters no longer
  shadow m/n/k (now `m_white`/`n_black`/`k_draw`); `find_min`'s
  start-row argument is `k_start` and `find_max`'s rows-remaining
  local is `d_rem`; historical `_v4`/`_v5` symbol suffixes dropped;
  the rxc kernel uses `m`, `c`, `n`, `rmarg`/`cmarg`,
  `rresid`/`cresid`, coordinates `i`, `j`, and path index `t`; the
  rxck kernel uses coordinates `i`, `j`, `l` (with `k` reserved for
  the layer count), Bartlett margin labels `ab`/`ac`/`bc`, and path
  index `t`. No exported API or numerical behavior changes; all
  kernels verified byte-identical on a 50-table fingerprint and the
  204-assertion test suite.

# zzfisher 0.2.0

* Synchronized all kernels with the research compendium (single
  source of truth). Kernel sources now carry the shape prefix
  (`rx2_*`, `rxc_*`, `rxck_*`); compiled internals are named
  `.rx2_*`, `.rxc_*`, `.rxck_*`.
* User-facing functions are now shape-dispatching entry points:
  a 2D matrix with 2 columns routes to the r x 2 kernel, a 2D
  matrix with more columns to the r x c kernel, and a 3D array to
  the r x c x k kernel (first-generation shared kernels for the
  latter two).
* Renamed `tree_memo_s4_c()` to `tree_s4_c()` to match the
  compendium naming convention.
* Added the manuscript benchmark driver as
  `inst/benchmarks/benchmark_manuscript_tables.R`.
* Extended the tinytest suite to cover r x c and r x c x k shapes.

# zzfisher 0.1.0

* Initial public release.
