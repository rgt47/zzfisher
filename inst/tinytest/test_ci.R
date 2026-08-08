library(tinytest)

# Stratified conditional-independence exact test (net_ci_cpp).

# With k = 1 the test reduces to the ordinary r x c Fisher exact
# test.
d1 <- matrix(c(1, 9, 11, 3), nrow = 2)
a1 <- array(d1, dim = c(2, 2, 1))
expect_true(
  abs(net_ci_cpp(a1)$p.value - fisher.test(d1)$p.value) < 1e-12,
  info = "k = 1 reduces to 2 x 2 fisher.test")

d2 <- matrix(c(2, 1, 0, 3, 1, 2, 1, 0, 2, 1, 1, 2), nrow = 3)
a2 <- array(d2, dim = c(3, 4, 1))
expect_true(
  abs(net_ci_cpp(a2)$p.value -
      fisher.test(d2, workspace = 2e6)$p.value) < 1e-12,
  info = "k = 1 reduces to 3 x 4 fisher.test")

# 2 x 2 x 2: brute force over the product of the two per-stratum
# fibers, from the definition.
brute_ci_22k <- function(arr, tol = 3.45254e-7) {
  k <- dim(arr)[3]
  fibers <- lapply(seq_len(k), function(l) {
    sl <- arr[, , l]
    r <- rowSums(sl); cc <- colSums(sl); n <- sum(sl)
    xs <- max(0, r[1] - cc[2]):min(r[1], cc[1])
    lp <- vapply(xs, function(x) {
      tab <- matrix(c(x, r[1] - x, cc[1] - x, n - r[1] - cc[1] + x),
                    nrow = 2)
      sum(lgamma(r + 1)) + sum(lgamma(cc + 1)) - lgamma(n + 1) -
        sum(lgamma(tab + 1))
    }, numeric(1))
    lp
  })
  obs_lp <- sum(vapply(seq_len(k), function(l) {
    sl <- arr[, , l]
    r <- rowSums(sl); cc <- colSums(sl); n <- sum(sl)
    sum(lgamma(r + 1)) + sum(lgamma(cc + 1)) - lgamma(n + 1) -
      sum(lgamma(sl + 1))
  }, numeric(1)))
  thresh <- obs_lp + log1p(tol)
  grids <- expand.grid(fibers)
  tot <- rowSums(as.matrix(grids))
  1 - sum(exp(tot)[tot > thresh])
}

set.seed(20260808)
for (rep in 1:5) {
  arr <- array(rpois(8, 3), dim = c(2, 2, 2))
  if (any(apply(arr, 3, function(s) any(rowSums(s) == 0) ||
                any(colSums(s) == 0)))) next
  expect_true(
    abs(net_ci_cpp(arr)$p.value - brute_ci_22k(arr)) < 1e-12,
    info = sprintf("2x2x2 brute force, rep %d", rep))
}

# Three strata.
arr3 <- array(c(5, 1, 2, 4, 4, 2, 3, 3, 1, 5, 4, 2),
              dim = c(2, 2, 3))
expect_true(
  abs(net_ci_cpp(arr3)$p.value - brute_ci_22k(arr3)) < 1e-12,
  info = "2x2x3 brute force")

# Shape guard.
expect_error(net_ci_cpp(matrix(1:4, 2)),
             info = "2D input rejected")
