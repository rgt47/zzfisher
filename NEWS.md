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
