# benchmark_kernel_comparison.R
#
# Times the three production kernels (tree_memo_c, tree_dp_cpp,
# net_dp_cpp) on the same equal-margin (m, rho) grid as Table 2 of
# the manuscript, supporting the Section 6 statement about their
# relative performance. Same seed and table-generation scheme as
# benchmark_manuscript_tables.R; each cell is run in a forked child
# with a wall-clock budget and recorded as NA on timeout.
#
# Output: ./derived_data/manuscript_kernel_comparison.csv (redirect
# with the ZZFISHER_BENCH_OUT environment variable) with median
# seconds per kernel per (m, rho) cell.
#
# Usage:
#   Rscript -e "source(system.file('benchmarks',
#     'benchmark_kernel_comparison.R', package = 'zzfisher'))"

suppressMessages(library(zzfisher))
suppressMessages({
  library(microbenchmark)
  library(parallel)
})

set.seed(20260411)

out_dir <- Sys.getenv("ZZFISHER_BENCH_OUT", "derived_data")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

bench <- function(call_fun, reps) {
  tm <- microbenchmark(call_fun(), times = reps, unit = "s")
  stats::median(tm$time) / 1e9
}

with_budget <- function(f, budget_s) {
  job <- mcparallel(f())
  res <- mccollect(job, wait = FALSE, timeout = budget_s)
  if (is.null(res)) {
    tools::pskill(job$pid)
    tools::pskill(job$pid, tools::SIGKILL)
    suppressWarnings(try(mccollect(job), silent = TRUE))
    return(NULL)
  }
  val <- res[[1]]
  if (inherits(val, "try-error")) return(NULL)
  val
}

rand_table <- function(m, margin_fun) {
  R <- margin_fun(m)
  y <- vapply(R, function(ri) sample.int(ri + 1L, 1L) - 1L, integer(1))
  cbind(y, R - y)
}

kernels <- list(
  tree_memo_c = function(dat) tree_memo_c(dat)$p.value,
  tree_dp_cpp = function(dat) tree_dp_cpp(dat)$p.value,
  net_dp_cpp  = function(dat) net_dp_cpp(dat)$p.value
)

ms <- c(3, 4, 5, 6, 7, 8, 9, 10, 12, 15, 20)
rs <- c(2, 3, 4, 5, 8, 10, 15, 20)
budget_s <- 3
grid <- expand.grid(m = ms, rho = rs)

cell <- function(m, rho) {
  dat <- rand_table(m, function(k) rep(rho, k))
  if (any(colSums(dat) == 0))
    return(setNames(rep(NA_real_, length(kernels)), names(kernels)))
  vapply(names(kernels), function(kn) {
    f <- kernels[[kn]]
    out <- with_budget(function() {
      t_probe <- bench(function() f(dat), 1L)
      reps <- if (t_probe < 0.02) 30L else if (t_probe < 0.2) 10L else 3L
      bench(function() f(dat), reps)
    }, budget_s)
    if (is.null(out)) NA_real_ else out
  }, numeric(1))
}

vals <- mapply(cell, grid$m, grid$rho)
comparison <- cbind(grid, t(vals))

write.csv(comparison,
          file.path(out_dir, "manuscript_kernel_comparison.csv"),
          row.names = FALSE)

ok <- stats::complete.cases(comparison)
cat(sprintf("cells complete for all kernels: %d of %d\n",
            sum(ok), nrow(comparison)))
cat(sprintf("tree_dp_cpp / tree_memo_c ratio: median %.2f, max %.2f\n",
            stats::median(comparison$tree_dp_cpp[ok] /
                          comparison$tree_memo_c[ok]),
            max(comparison$tree_dp_cpp[ok] /
                comparison$tree_memo_c[ok])))
cat(sprintf("net_dp_cpp  / tree_memo_c ratio: median %.2f, max %.2f\n",
            stats::median(comparison$net_dp_cpp[ok] /
                          comparison$tree_memo_c[ok]),
            max(comparison$net_dp_cpp[ok] /
                comparison$tree_memo_c[ok])))
