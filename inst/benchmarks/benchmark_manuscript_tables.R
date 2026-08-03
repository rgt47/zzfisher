# benchmark_manuscript_tables.R
#
# Regenerates the three timing tables reported in
# analysis/report/01-rx2/report.Rmd (Tables 1-3):
#
#   Table 1  win rate of tree_memo_c against fisher.test on 300 random
#            small tables, by number of rows m
#   Table 2  time ratio (fisher.test / tree_memo_c) over a grid of m and
#            common row margin rho
#   Table 3  median run time for two large diverse-margin configurations
#
# The manuscript reports median point estimates from a single run of this
# driver. Absolute timings depend on hardware, BLAS, and R version; the
# qualitative regime structure is stable. This script replicates each
# timing and reports the median together with the interquartile range so
# that the run-to-run variability discussed in Section 6 is visible.
#
# Robustness. The compiled kernels are atomic C calls and cannot be
# interrupted by an R-level time limit. Each potentially expensive call is
# therefore run in a forked child process (parallel::mcparallel) that is
# killed if it exceeds a per-cell wall-clock budget; killed cells are
# recorded as NA, reproducing the "---" / blank entries of Table 2 for
# configurations that are intractable for the tree kernel.
#
# Usage:
#   Rscript -e "source(system.file('benchmarks',
#     'benchmark_manuscript_tables.R', package = 'zzfisher'))"
#
# Output: writes three CSV files to ./derived_data/ and prints a
# dispersion summary to stdout. Set the ZZFISHER_BENCH_OUT environment
# variable to redirect the output directory.

suppressMessages(library(zzfisher))
suppressMessages({
  library(microbenchmark)
  library(parallel)
})

set.seed(20260411)

out_dir <- Sys.getenv("ZZFISHER_BENCH_OUT", "derived_data")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- timing helpers -------------------------------------------------

# Median and IQR (seconds) of `reps` evaluations of call_fun().
bench <- function(call_fun, reps) {
  tm <- microbenchmark(call_fun(), times = reps, unit = "s")
  c(median = stats::median(tm$time) / 1e9, iqr = stats::IQR(tm$time) / 1e9)
}

# Run f() in a forked child, killing it if it exceeds budget_s seconds.
# Returns f()'s value, or NULL on timeout / error.
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

# Adaptive replication count from a single probe time.
pick_reps <- function(t_probe, fast = 30L, mid = 10L, slow = 3L) {
  if (t_probe < 0.02) fast else if (t_probe < 0.2) mid else slow
}

# A random r x 2 table with given number of rows and a margin sampler.
rand_table <- function(m, margin_fun) {
  R <- margin_fun(m)
  y <- vapply(R, function(ri) sample.int(ri + 1L, 1L) - 1L, integer(1))
  cbind(y, R - y)
}

# ---- Table 1: win rate on 300 random small tables -------------------
cat("Table 1: 300 random small tables ...\n")
table1 <- local({
  ms <- 3:7
  res <- lapply(ms, function(m) {
    wins <- 0L
    ratios <- numeric(0)
    for (i in seq_len(60L)) {
      dat <- rand_table(m, function(k) sample(2:30, k, replace = TRUE))
      if (any(colSums(dat) == 0)) next
      t_fish <- tryCatch(bench(function() fisher.test(dat)$p.value, 25L)["median"],
                         error = function(e) NA_real_)
      t_tree <- bench(function() tree_memo_c(dat)$p.value, 25L)["median"]
      if (is.na(t_fish)) next
      ratios <- c(ratios, t_fish / t_tree)
      if (t_tree < t_fish) wins <- wins + 1L
    }
    data.frame(m = m, n = length(ratios), faster = wins,
               ratio_median = stats::median(ratios),
               ratio_min = min(ratios), ratio_max = max(ratios))
  })
  do.call(rbind, res)
})
write.csv(table1, file.path(out_dir, "manuscript_table1_winrate.csv"),
          row.names = FALSE)

# ---- Table 2: crossover grid (fisher.test / tree_memo_c) ------------
cat("Table 2: crossover grid (forked, budgeted) ...\n")
table2 <- local({
  ms <- c(3, 4, 5, 6, 7, 8, 9, 10, 12, 15, 20)
  rs <- c(2, 3, 4, 5, 8, 10, 15, 20)
  budget_s <- 4
  grid <- expand.grid(m = ms, rho = rs)
  cell <- function(m, rho) {
    dat <- rand_table(m, function(k) rep(rho, k))
    if (any(colSums(dat) == 0)) return(c(ratio = NA, tree = NA, tree_iqr = NA))
    out <- with_budget(function() {
      t_probe <- bench(function() tree_memo_c(dat)$p.value, 1L)["median"]
      reps <- pick_reps(t_probe)
      tr <- bench(function() tree_memo_c(dat)$p.value, reps)
      fi <- tryCatch(bench(function() fisher.test(dat)$p.value, 30L)["median"],
                     error = function(e) NA_real_)
      c(ratio = unname(fi / tr["median"]), tree = unname(tr["median"]),
        tree_iqr = unname(tr["iqr"]))
    }, budget_s)
    if (is.null(out)) c(ratio = NA, tree = NA, tree_iqr = NA) else out
  }
  vals <- mapply(cell, grid$m, grid$rho)
  cbind(grid, t(vals))
})
write.csv(table2, file.path(out_dir, "manuscript_table2_crossover.csv"),
          row.names = FALSE)

# ---- Table 3: large diverse-margin configurations -------------------
cat("Table 3: large diverse-margin configs (this is the slow one) ...\n")
table3 <- local({
  configs <- list(
    list(label = "m=12, rho in 5-20, spread", m = 12, lo = 5, hi = 20),
    list(label = "m=15, rho in 5-15, spread", m = 15, lo = 5, hi = 15)
  )
  res <- lapply(configs, function(cfg) {
    dat <- rand_table(cfg$m,
                      function(k) sample(cfg$lo:cfg$hi, k, replace = TRUE))
    tr <- with_budget(function() {
      t_probe <- bench(function() tree_memo_c(dat)$p.value, 1L)["median"]
      reps <- if (t_probe < 1) 5L else 2L
      bench(function() tree_memo_c(dat)$p.value, reps)
    }, budget_s = 400)
    fi <- tryCatch(bench(function() fisher.test(dat)$p.value, 30L),
                   error = function(e) c(median = NA, iqr = NA))
    data.frame(config = cfg$label,
               tree = if (is.null(tr)) NA else unname(tr["median"]),
               tree_iqr = if (is.null(tr)) NA else unname(tr["iqr"]),
               fisher = unname(fi["median"]),
               fisher_iqr = unname(fi["iqr"]))
  })
  do.call(rbind, res)
})
write.csv(table3, file.path(out_dir, "manuscript_table3_large.csv"),
          row.names = FALSE)

# ---- dispersion summary --------------------------------------------
cat("\n================ DISPERSION SUMMARY ================\n")
cat("\nTable 1 (win rate by m):\n")
print(table1, row.names = FALSE)
cat("\nTable 2 (ratio fisher/tree; tree median seconds; tree IQR):\n")
t2 <- table2
t2$ratio <- round(t2$ratio, 2)
t2$tree <- signif(t2$tree, 3)
t2$tree_iqr <- signif(t2$tree_iqr, 2)
print(t2, row.names = FALSE)
cat("\nTable 3 (large configs, seconds):\n")
print(table3, row.names = FALSE)
cat("\nWrote manuscript_table{1,2,3}_*.csv to", out_dir, "\n")
