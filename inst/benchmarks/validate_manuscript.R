# validate_manuscript.R
#
# Random-table validation study reported in the manuscript's
# Validation subsection (Section 6). For each m in 2..7, generates
# random m x 2 tables with margins sampled uniformly from 2..30 and
# compares the three production kernels (tree_memo_c, tree_dp_cpp,
# net_dp_cpp) against stats::fisher.test().
#
# Tables on which fisher.test() fails at its default workspace
# (FEXACT 'LDSTP too small') are validated against
# fisher.test(workspace = 2e8) instead, and counted separately.
#
# Output: ./derived_data/manuscript_validation.csv (redirect with
# the ZZFISHER_BENCH_OUT environment variable) with
# one row per m: number of tables, default-workspace failures,
# and the maximum absolute p-value discrepancy per kernel across
# both oracles.
#
# Usage:
#   Rscript -e "source(system.file('benchmarks',
#     'validate_manuscript.R', package = 'zzfisher'))"

suppressMessages(library(zzfisher))

set.seed(20260802)

out_dir <- Sys.getenv("ZZFISHER_BENCH_OUT", "derived_data")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

rand_table <- function(m) {
  R <- sample(2:30, m, replace = TRUE)
  y <- vapply(R, function(ri) sample.int(ri + 1L, 1L) - 1L, integer(1))
  cbind(y, R - y)
}

kernels <- list(
  tree_memo_c = function(dat) tree_memo_c(dat)$p.value,
  tree_dp_cpp = function(dat) tree_dp_cpp(dat)$p.value,
  net_dp_cpp  = function(dat) net_dp_cpp(dat)$p.value
)

rows <- lapply(2:7, function(m) {
  n_tab <- 200L
  n_fail <- 0L
  worst <- setNames(numeric(length(kernels)), names(kernels))
  n_used <- 0L
  for (i in seq_len(n_tab)) {
    dat <- rand_table(m)
    if (any(colSums(dat) == 0)) next
    oracle <- tryCatch(fisher.test(dat)$p.value, error = function(e) NULL)
    if (is.null(oracle)) {
      n_fail <- n_fail + 1L
      oracle <- tryCatch(fisher.test(dat, workspace = 2e8)$p.value,
                         error = function(e) NULL)
      if (is.null(oracle)) next
    }
    n_used <- n_used + 1L
    for (kn in names(kernels)) {
      worst[kn] <- max(worst[kn], abs(kernels[[kn]](dat) - oracle))
    }
  }
  data.frame(m = m, n = n_used, fisher_default_fail = n_fail,
             max_diff_tree_memo_c = worst["tree_memo_c"],
             max_diff_tree_dp_cpp = worst["tree_dp_cpp"],
             max_diff_net_dp_cpp = worst["net_dp_cpp"],
             row.names = NULL)
})
validation <- do.call(rbind, rows)

stopifnot(nrow(validation) == 6L, all(validation$n > 0))

write.csv(validation, file.path(out_dir, "manuscript_validation.csv"),
          row.names = FALSE)
cat("Validation summary:\n")
print(validation, row.names = FALSE)
cat(sprintf("\nOverall max discrepancy: %.3e\n",
            max(validation[, grep("max_diff", names(validation))])))
