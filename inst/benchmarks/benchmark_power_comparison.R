# benchmark_power_comparison.R
#
# Times and cross-validates the competing exact power-function
# implementations that ship in zzfisher:
#
#   fisher_power_1s  fisher_power(alternative = "less")   -- pure R,
#                    direct enumeration via phyper()/dbinom() vectors
#   fisher_power_2s  fisher_power(alternative = "two.sided") -- pure
#                    R, elementwise unimodal-tail scan (the most
#                    expensive path: see docs on binary-search
#                    opportunities not yet implemented)
#   fxpower_exact    fxpower(..., eps = 0)   -- compiled C++,
#                    mode-centered zigzag traversal, no trimming
#   fxpower_trim     fxpower(..., eps = 1e-6) -- same kernel with
#                    eps-trimming (a quantified early-termination
#                    bound, this package's fastest path)
#   fxpower_r_ref    fxpower_r() -- pure R reference for fxpower,
#                    capped at n <= 500 per the manuscript's own
#                    feasibility note
#   exact_pkg        Exact::power.exact.test(method = "fisher"),
#                    an independent third-party implementation, run
#                    if the Exact package is installed
#
# fxpower()/Exact::power.exact.test() are one-sided only, so every
# comparison uses fisher_power(alternative = "less") as the matching
# direction (verified in inst/tinytest/test_fxpower.R to be bitwise
# identical to fxpower(eps = 0) on a separate case grid). fisher_power
# two-sided has no one-sided counterpart here and is benchmarked on
# its own, restricted to smaller n where its O(range^2)-per-m
# elementwise scan stays tractable.
#
# Output: ./derived_data/power_comparison.csv (redirect with the
# ZZFISHER_BENCH_OUT environment variable), one row per
# (config, algorithm) with elapsed time and the computed power, plus
# a printed cross-validation summary (max discrepancy between
# algorithms that should agree) and speed-ratio summary.
#
# Usage:
#   Rscript -e "source(system.file('benchmarks',
#     'benchmark_power_comparison.R', package = 'zzfisher'))"

suppressMessages(library(zzfisher))
suppressMessages({
  library(microbenchmark)
  library(parallel)
})
has_exact <- suppressMessages(requireNamespace("Exact", quietly = TRUE))

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

# --- algorithms ------------------------------------------------------

algorithms <- list(
  fisher_power_1s = list(
    fn = function(cfg) {
      fisher_power(cfg$n1, cfg$n2, p1 = cfg$p1, p2 = cfg$p2,
                   alpha = cfg$alpha, alternative = "less")$power
    },
    n_cap = Inf, group = "one-sided"
  ),
  fisher_power_2s = list(
    fn = function(cfg) {
      fisher_power(cfg$n1, cfg$n2, p1 = cfg$p1, p2 = cfg$p2,
                   alpha = cfg$alpha, alternative = "two.sided")$power
    },
    n_cap = 200, group = "two-sided"
  ),
  fxpower_exact = list(
    fn = function(cfg) {
      fxpower(cfg$n1, cfg$n2, cfg$p1, cfg$p2,
              alpha = cfg$alpha, eps = 0)$power
    },
    n_cap = Inf, group = "one-sided"
  ),
  fxpower_trim = list(
    fn = function(cfg) {
      fxpower(cfg$n1, cfg$n2, cfg$p1, cfg$p2,
              alpha = cfg$alpha, eps = 1e-6)$power
    },
    n_cap = Inf, group = "one-sided"
  ),
  fxpower_r_ref = list(
    fn = function(cfg) {
      fxpower_r(cfg$n1, cfg$n2, cfg$p1, cfg$p2, alpha = cfg$alpha)
    },
    n_cap = 500, group = "one-sided"
  )
)

if (has_exact) {
  algorithms$exact_pkg <- list(
    fn = function(cfg) {
      suppressMessages(
        Exact::power.exact.test(
          p1 = cfg$p1, p2 = cfg$p2, n1 = cfg$n1, n2 = cfg$n2,
          alpha = cfg$alpha, alternative = "less", method = "fisher"
        )$power
      )
    },
    n_cap = Inf, group = "one-sided"
  )
}

# --- configurations ---------------------------------------------------

configs <- list(
  list(n1 = 50,   n2 = 50,   p1 = 0.10, p2 = 0.30, alpha = 0.05),
  list(n1 = 100,  n2 = 100,  p1 = 0.05, p2 = 0.10, alpha = 0.05),
  list(n1 = 200,  n2 = 100,  p1 = 0.05, p2 = 0.10, alpha = 0.05),
  list(n1 = 500,  n2 = 500,  p1 = 0.05, p2 = 0.10, alpha = 0.05),
  list(n1 = 1000, n2 = 500,  p1 = 0.05, p2 = 0.10, alpha = 0.05),
  list(n1 = 2000, n2 = 2000, p1 = 0.05, p2 = 0.10, alpha = 0.05),
  list(n1 = 5000, n2 = 5000, p1 = 0.05, p2 = 0.10, alpha = 0.05)
)
if (nzchar(Sys.getenv("ZZFISHER_BENCH_SMOKE_TEST"))) {
  configs <- list(
    list(n1 = 50, n2 = 50, p1 = 0.10, p2 = 0.30, alpha = 0.05),
    list(n1 = 200, n2 = 100, p1 = 0.05, p2 = 0.10, alpha = 0.05)
  )
}

budget_s <- 15

run_one <- function(cfg, alg_name) {
  alg <- algorithms[[alg_name]]
  if (max(cfg$n1, cfg$n2) > alg$n_cap) {
    return(data.frame(
      n1 = cfg$n1, n2 = cfg$n2, p1 = cfg$p1, p2 = cfg$p2,
      algorithm = alg_name, group = alg$group,
      time_s = NA_real_, power = NA_real_,
      skipped = "n exceeds algorithm's feasibility cap"
    ))
  }

  out <- with_budget(function() {
    power_val <- alg$fn(cfg)
    t_probe <- bench(function() alg$fn(cfg), 1L)
    reps <- if (t_probe < 0.02) 20L else if (t_probe < 0.5) 5L else 1L
    time_s <- bench(function() alg$fn(cfg), reps)
    list(power = power_val, time_s = time_s)
  }, budget_s)

  if (is.null(out)) {
    data.frame(
      n1 = cfg$n1, n2 = cfg$n2, p1 = cfg$p1, p2 = cfg$p2,
      algorithm = alg_name, group = alg$group,
      time_s = NA_real_, power = NA_real_,
      skipped = sprintf("exceeded %ds budget", budget_s)
    )
  } else {
    data.frame(
      n1 = cfg$n1, n2 = cfg$n2, p1 = cfg$p1, p2 = cfg$p2,
      algorithm = alg_name, group = alg$group,
      time_s = out$time_s, power = out$power,
      skipped = NA_character_
    )
  }
}

rows <- list()
for (cfg in configs) {
  for (alg_name in names(algorithms)) {
    rows[[length(rows) + 1L]] <- run_one(cfg, alg_name)
  }
}
comparison <- do.call(rbind, rows)

write.csv(comparison, file.path(out_dir, "power_comparison.csv"),
          row.names = FALSE)

# --- cross-validation: algorithms in the same group should agree -----

cat("\n--- cross-validation (one-sided group) ---\n")
one_sided <- comparison[comparison$group == "one-sided" &
                         is.na(comparison$skipped), ]
for (cfg in configs) {
  sub <- one_sided[one_sided$n1 == cfg$n1 & one_sided$n2 == cfg$n2 &
                    one_sided$p1 == cfg$p1 & one_sided$p2 == cfg$p2, ]
  if (nrow(sub) < 2) next
  ref <- sub$power[sub$algorithm == "fxpower_exact"]
  if (length(ref) == 0) next
  diffs <- abs(sub$power - ref)
  names(diffs) <- sub$algorithm
  cat(sprintf(
    "n1=%d n2=%d p1=%.2f p2=%.2f  max |diff from fxpower_exact| = %.2e\n",
    cfg$n1, cfg$n2, cfg$p1, cfg$p2, max(diffs)
  ))
}

# --- speed ratios ------------------------------------------------------

cat("\n--- speed ratios (relative to fxpower_trim) ---\n")
wide_time <- reshape(
  comparison[, c("n1", "n2", "p1", "p2", "algorithm", "time_s")],
  idvar = c("n1", "n2", "p1", "p2"), timevar = "algorithm",
  direction = "wide"
)
names(wide_time) <- sub("^time_s\\.", "", names(wide_time))
if (all(c("fxpower_trim", "fxpower_exact") %in% names(wide_time))) {
  ok <- stats::complete.cases(wide_time[, c("fxpower_trim", "fxpower_exact")])
  cat(sprintf(
    "fxpower_exact / fxpower_trim:      median %.1fx, max %.1fx\n",
    stats::median(wide_time$fxpower_exact[ok] / wide_time$fxpower_trim[ok]),
    max(wide_time$fxpower_exact[ok] / wide_time$fxpower_trim[ok])
  ))
}
if (all(c("fxpower_trim", "fisher_power_1s") %in% names(wide_time))) {
  ok <- stats::complete.cases(wide_time[, c("fxpower_trim", "fisher_power_1s")])
  cat(sprintf(
    "fisher_power_1s / fxpower_trim:    median %.1fx, max %.1fx\n",
    stats::median(wide_time$fisher_power_1s[ok] / wide_time$fxpower_trim[ok]),
    max(wide_time$fisher_power_1s[ok] / wide_time$fxpower_trim[ok])
  ))
}

cat(sprintf("\nfull results written to %s\n",
            file.path(out_dir, "power_comparison.csv")))
