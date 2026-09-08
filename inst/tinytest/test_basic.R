# Cross-kernel agreement and input validation. This file previously
# held a single expect_true(TRUE), which asserts nothing.

kern <- c("tree_dp", "tree_memo", "tree_memo_c", "tree_memo_cpp",
          "tree_dp_cpp", "tree_s4_c", "net_dp_cpp", "net_sweep_cpp",
          "net_vander_cpp")

# The package ships nine kernels for one quantity. Their whole point is
# that they agree, so that is what is asserted: against each other and
# against stats::fisher.test, over an exhaustive set of small tables.
tabs <- list()
for (a in 0:4) for (b in 0:4) for (cc in 0:4) for (dd in 0:4) {
  m <- matrix(c(a, b, cc, dd), 2, byrow = TRUE)
  if (all(rowSums(m) > 0) && all(colSums(m) > 0)) {
    tabs[[length(tabs) + 1L]] <- m
  }
}
for (m in tabs) {
  ref <- stats::fisher.test(m)$p.value
  ps <- vapply(kern, function(f) get(f)(m)$p.value, numeric(1))
  expect_true(max(abs(ps - ref)) < 1e-12,
    info = "every kernel matches stats::fisher.test on 2x2")
  expect_true(diff(range(ps)) < 1e-12,
    info = "the kernels agree with each other on 2x2")
}

# The same, for r x 2 and r x c shapes where the tree and network
# algorithms actually diverge from the 2x2 special case.
set.seed(101)
for (i in 1:40) {
  r <- sample(3:4, 1)
  repeat {
    m <- matrix(sample(0:5, 2 * r, TRUE), nrow = r)
    if (all(rowSums(m) > 0) && all(colSums(m) > 0)) break
  }
  ref <- stats::fisher.test(m)$p.value
  ps <- vapply(kern, function(f) get(f)(m)$p.value, numeric(1))
  expect_true(max(abs(ps - ref)) < 1e-10,
    info = "every kernel matches stats::fisher.test on r x 2")
}

# Malformed input must be refused, and refused the same way by all of
# them. A table containing -1 previously returned 0.25 from tree_dp and
# tree_memo_c, 0 from net_vander_cpp, and an error from net_sweep_cpp.
bad <- list(
  negative = matrix(c(-1, 2, 3, 4), 2),
  fractional = matrix(c(1.5, 2, 3, 4), 2),
  missing = matrix(c(NA, 2, 3, 4), 2),
  infinite = matrix(c(Inf, 2, 3, 4), 2),
  character = matrix(c("a", "b", "c", "d"), 2),
  one_column = matrix(c(1, 2), 2, 1)
)
for (nm in names(bad)) {
  for (f in kern) {
    expect_error(get(f)(bad[[nm]]),
      info = paste(f, "refuses", nm, "input"))
  }
}

# Valid but degenerate tables are still accepted: zero cells and an
# empty row are legitimate counts.
expect_true(is.numeric(tree_dp(matrix(c(0, 4, 4, 0), 2, byrow = TRUE))$p.value),
  info = "zero cells are accepted")
expect_true(is.numeric(tree_dp(matrix(c(0, 0, 3, 4), 2, byrow = TRUE))$p.value),
  info = "an empty row is accepted")

# net_ci_cpp documents that with k = 1 it reduces to the ordinary
# r x c Fisher exact test.
set.seed(102)
for (i in 1:20) {
  r <- sample(2:4, 1)
  cc <- sample(2:3, 1)
  repeat {
    m <- matrix(sample(0:6, r * cc, TRUE), r, cc)
    if (all(rowSums(m) > 0) && all(colSums(m) > 0)) break
  }
  expect_equal(net_ci_cpp(array(m, dim = c(r, cc, 1)))$p.value,
               stats::fisher.test(m)$p.value, tolerance = 1e-10,
    info = "net_ci_cpp with k = 1 is the r x c Fisher exact test")
}
expect_error(net_ci_cpp(array(c(-1, 1, 1, 3), dim = c(2, 2, 1))),
  info = "net_ci_cpp validates its counts too")
