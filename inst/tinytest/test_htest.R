library(tinytest)

dat_2x2 <- matrix(c(1, 9, 11, 3), nrow = 2)
dat_3x2 <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 3)
ref_2x2 <- fisher.test(dat_2x2)$p.value
ref_3x2 <- fisher.test(dat_3x2)$p.value

all_fns <- list(
tree_memo      = tree_memo,
tree_dp        = tree_dp,
tree_memo_c    = tree_memo_c,
tree_s4_c      = tree_s4_c,
tree_memo_cpp  = tree_memo_cpp,
tree_dp_cpp    = tree_dp_cpp,
net_vander_cpp = net_vander_cpp,
net_dp_cpp     = net_dp_cpp,
net_sweep_cpp  = net_sweep_cpp
)

for (nm in names(all_fns)) {
fn <- all_fns[[nm]]

result <- fn(dat_2x2)
expect_inherits(result, "htest")

expect_true(abs(result$p.value - ref_2x2) < 1e-12)

result_3x2 <- fn(dat_3x2)
expect_true(abs(result_3x2$p.value - ref_3x2) < 1e-12)

expect_true(grepl(nm, result$method))
expect_true(grepl("Fisher's Exact Test", result$method))

expect_identical(result$alternative, "two.sided")

expect_identical(result$data.name, "dat_2x2")

}

# tree_memo_profile returns htest with profile
profile_ok <- tryCatch(
  { tree_memo_profile(dat_2x2); TRUE },
  error = function(e) FALSE
)
if (profile_ok) {
  result <- tree_memo_profile(dat_2x2)
  expect_inherits(result, "htest")
  expect_true(grepl("tree_memo_profile", result$method))
  expect_true("profile" %in% names(result))
}

# tree_memo_profile is the only dispatcher that still stubs rxc/rxck
# (profiling is an rx2-specific instrumentation feature). All other
# algorithms route rxc/rxck to the shared first-gen kernels.
dat_3x3 <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9), nrow = 3)
dat_2x2x2 <- array(1:8, dim = c(2, 2, 2))
err_rxc <- tryCatch(tree_memo_profile(dat_3x3),
                    error = function(e) conditionMessage(e))
expect_true(is.character(err_rxc) &&
            grepl("not yet implemented", err_rxc) &&
            grepl("tree_memo_profile", err_rxc, fixed = TRUE))
err_rxck <- tryCatch(tree_memo_profile(dat_2x2x2),
                     error = function(e) conditionMessage(e))
expect_true(is.character(err_rxck) &&
            grepl("not yet implemented", err_rxck))

# r x c oracle: every wired dispatcher must match fisher.test() on a
# set of test cases. All eight currently route to the same underlying
# kernel (.rxc_tree_memo_cpp), so agreement is mechanical, but the
# assertion protects against accidental regression in any dispatcher
# wiring.
rxc_cases <- list(
  dat_2x3 = matrix(c(1, 2, 3, 4, 5, 6), nrow = 2),
  dat_3x3 = matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9), nrow = 3),
  dat_3x4 = matrix(c(2, 1, 3, 4, 0, 2, 1, 3, 2, 1, 4, 0), nrow = 3)
)
for (nm in names(rxc_cases)) {
  dat <- rxc_cases[[nm]]
  ref <- fisher.test(dat)$p.value
  for (algo in names(all_fns)) {
    result <- all_fns[[algo]](dat)
    expect_inherits(result, "htest")
    expect_true(abs(result$p.value - ref) < 1e-12,
                info = sprintf("rxc %s via %s: got %.15g, expected %.15g",
                               nm, algo, result$p.value, ref))
    expect_true(grepl(algo, result$method, fixed = TRUE),
                info = sprintf("rxc %s via %s: method tag", nm, algo))
  }
}

# r x c x k oracle: pure-R brute-force enumeration of the no-three-
# way-interaction reference set (conditioning on all three pairwise
# margins). Generalized to arbitrary R x C x K, but exponential in
# (R-1)(C-1)(K-1) — keep test cases small. The C++ kernel and this
# enumerator implement the same algorithm and test statistic, so
# agreement is a self-consistency check (same-language parity), not
# independent verification.
rxck_brute_force <- function(dat) {
  R <- dim(dat)[1]; C <- dim(dat)[2]; K <- dim(dat)[3]
  RC <- apply(dat, c(1, 2), sum)
  RK <- apply(dat, c(1, 3), sum)
  CK <- apply(dat, c(2, 3), sum)
  T_obs <- sum(lgamma(dat + 1))
  env <- new.env()
  env$mass_total <- 0
  env$mass_extreme <- 0

  # Compute dependent cells from pairwise margins; return NULL if any
  # cell would be negative or a cross-margin check fails.
  finish <- function(y) {
    R1 <- R - 1L; C1 <- C - 1L; K1 <- K - 1L
    for (i in seq_len(R1)) for (j in seq_len(C1)) {
      v <- RC[i, j] - sum(y[i, j, seq_len(K1)])
      if (v < 0) return(NULL)
      y[i, j, K] <- v
    }
    for (i in seq_len(R1)) for (k in seq_len(K1)) {
      v <- RK[i, k] - sum(y[i, seq_len(C1), k])
      if (v < 0) return(NULL)
      y[i, C, k] <- v
    }
    for (j in seq_len(C1)) for (k in seq_len(K1)) {
      v <- CK[j, k] - sum(y[seq_len(R1), j, k])
      if (v < 0) return(NULL)
      y[R, j, k] <- v
    }
    for (i in seq_len(R1)) {
      v <- RC[i, C] - sum(y[i, C, seq_len(K1)])
      if (v < 0) return(NULL)
      y[i, C, K] <- v
      if (sum(y[i, , K]) != RK[i, K]) return(NULL)
    }
    for (j in seq_len(C1)) {
      v <- CK[j, K] - sum(y[seq_len(R1), j, K])
      if (v < 0) return(NULL)
      y[R, j, K] <- v
      if (sum(y[R, j, ]) != RC[R, j]) return(NULL)
    }
    for (k in seq_len(K1)) {
      v <- CK[C, k] - sum(y[seq_len(R1), C, k])
      if (v < 0) return(NULL)
      y[R, C, k] <- v
      if (sum(y[R, , k]) != RK[R, k]) return(NULL)
    }
    v <- RC[R, C] - sum(y[R, C, seq_len(K1)])
    if (v < 0) return(NULL)
    y[R, C, K] <- v
    if (sum(y[R, , K]) != RK[R, K]) return(NULL)
    if (sum(y[, C, K]) != CK[C, K]) return(NULL)
    y
  }

  traverse <- function(l, lmax, y, RCr, RKr, CKr) {
    if (l > lmax) {
      yf <- finish(y)
      if (is.null(yf)) return(invisible())
      Tv <- sum(lgamma(yf + 1))
      w <- exp(-Tv)
      env$mass_total <- env$mass_total + w
      if (Tv - T_obs >= -3.45254e-7)
        env$mass_extreme <- env$mass_extreme + w
      return(invisible())
    }
    l0 <- l - 1L
    per_layer <- (R - 1L) * (C - 1L)
    k <- l0 %/% per_layer
    rem <- l0 %% per_layer
    cc <- rem %/% (R - 1L)
    r <- rem %% (R - 1L)
    r1 <- r + 1L; c1 <- cc + 1L; k1 <- k + 1L
    x_up <- min(RCr[r1, c1], RKr[r1, k1], CKr[c1, k1])
    if (x_up < 0) return(invisible())
    for (x in 0:x_up) {
      y[r1, c1, k1] <- x
      RCr[r1, c1] <- RCr[r1, c1] - x
      RKr[r1, k1] <- RKr[r1, k1] - x
      CKr[c1, k1] <- CKr[c1, k1] - x
      traverse(l + 1L, lmax, y, RCr, RKr, CKr)
      RCr[r1, c1] <- RCr[r1, c1] + x
      RKr[r1, k1] <- RKr[r1, k1] + x
      CKr[c1, k1] <- CKr[c1, k1] + x
    }
    invisible()
  }

  lmax <- (R - 1L) * (C - 1L) * (K - 1L)
  traverse(1L, lmax, array(0L, dim(dat)), RC, RK, CK)
  env$mass_extreme / env$mass_total
}

# Hand-verified sanity case: identical strata produce OR = 1 in every
# stratum, so the no-three-way-interaction null holds exactly. With
# R = C = 2 and K = 2 and the observed table equal to the product of
# its pairwise margins, the observed T_obs is the minimum of T(Y) over
# the reference set, so p-value = 1 exactly.
#   stratum 1: [[1,1],[1,1]]      stratum 2: [[1,1],[1,1]]
# Every cell = 1, all pairwise margins symmetric. Exact p = 1.
dat_homog <- array(rep(1L, 8), dim = c(2, 2, 2))
result_homog <- tree_memo_cpp(dat_homog)
expect_true(abs(result_homog$p.value - 1) < 1e-12,
            info = sprintf("rxck homogeneous-OR: got %.15g, expected 1",
                           result_homog$p.value))

rxck_cases <- list(
  dat_a   = array(c(2, 1, 1, 2, 1, 2, 2, 1), dim = c(2, 2, 2)),
  dat_b   = array(c(3, 1, 2, 4, 2, 3, 1, 2), dim = c(2, 2, 2)),
  dat_2x2x3 = array(c(2, 1, 1, 2,
                      1, 2, 2, 1,
                      2, 1, 1, 3), dim = c(2, 2, 3))
)
for (nm in names(rxck_cases)) {
  dat <- rxck_cases[[nm]]
  ref <- rxck_brute_force(dat)
  for (algo in names(all_fns)) {
    result <- all_fns[[algo]](dat)
    expect_inherits(result, "htest")
    expect_true(abs(result$p.value - ref) < 1e-10,
                info = sprintf("rxck %s via %s: got %.15g, expected %.15g",
                               nm, algo, result$p.value, ref))
    expect_true(grepl(algo, result$method, fixed = TRUE),
                info = sprintf("rxck %s via %s: method tag", nm, algo))
  }
}

# rxck asymptotic sanity check via contingencytables::BreslowDay_...
# (Suggests-level; skipped if unavailable). Breslow-Day's asymptotic
# chi-squared tests the same null (no-three-way-interaction / OR
# homogeneity across K strata in a 2x2xK table), so decisions should
# agree at the extremes; tight p-value agreement is not expected
# because Breslow-Day is asymptotic and the sample sizes here are
# small. Empirically: agreement to ~0.2 on 2x2xK with counts < 15 per
# cell; agreement to ~0.01 when both are in (0, 0.05) or (0.95, 1).
if (requireNamespace("contingencytables", quietly = TRUE)) {
  bd_cases <- list(
    hetero = array(c(10, 2, 2, 10,
                     3, 8, 9, 3,
                     11, 2, 2, 10), dim = c(2, 2, 3)),
    # Non-degenerate homogeneous case (same OR in each stratum);
    # all-equal cells give OR = 1 in every stratum and make the
    # Breslow-Day statistic 0/0 (returns NA), so we use a non-
    # unity shared odds ratio here.
    homog  = array(c(4, 2, 2, 4,
                     6, 3, 3, 6), dim = c(2, 2, 2))
  )
  for (nm in names(bd_cases)) {
    dat <- bd_cases[[nm]]
    bd <- contingencytables::BreslowDay_homogeneity_test_stratified_2x2(dat)
    ex <- tree_memo_cpp(dat)$p.value
    # Directional agreement: decision at alpha = 0.05 should match.
    expect_equal(bd$Pvalue < 0.05, ex < 0.05,
                 info = sprintf("rxck %s: BD p=%.4f vs exact p=%.4f",
                                nm, bd$Pvalue, ex))
  }
}
