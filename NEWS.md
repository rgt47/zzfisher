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
