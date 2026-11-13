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

# --- large-table regressions -------------------------------------
#
# lfact is a fixed array of MAX_N + 1 doubles and its fill loop runs to
# the table's total count, so a total above MAX_N wrote past the end of
# the struct. rx2_tree_dp.cpp had carried a bounds check all along;
# rx2_tree_memo_cpp.cpp had lost it. A 2x2 table of 6000s -- 24000
# observations, an unremarkable size for a clinical dataset, with every
# count a positive whole number -- segfaulted the session.
#
# Two of the network kernels returned NaN on the same table. NaN in a
# p-value slot is worse than an error: it compares FALSE against any
# threshold, so the test silently stops rejecting and nothing in the
# output says so.
big <- matrix(rep(6000, 4), 2)
for (f in kern) {
  res <- tryCatch(get(f)(big), error = function(e) e)
  if (inherits(res, "error")) {
    expect_true(nzchar(conditionMessage(res)),
      info = paste(f, "refuses a large table with a message"))
  } else {
    expect_true(is.finite(res$p.value),
      info = paste(f, "returns a finite p-value on a large table"))
    expect_true(res$p.value >= 0 && res$p.value <= 1,
      info = paste(f, "returns a probability on a large table"))
  }
}

# Whatever the kernel does, it must not hand back a non-finite value
# dressed as a p-value: make_htest() is the single point every kernel
# returns through. It is internal, and this file is shared between two
# packages with different names, so it is taken from the namespace of an
# exported kernel rather than through a named pkg::: call.
make_htest <- get("make_htest", envir = environment(tree_dp))
expect_error(make_htest(NaN, "d"), "rather than a p-value",
  info = "NaN is refused as a p-value")
expect_error(make_htest(NA_real_, "d"), "rather than a p-value",
  info = "NA is refused as a p-value")
expect_error(make_htest(Inf, "d"), "rather than a p-value",
  info = "Inf is refused as a p-value")
expect_error(make_htest(1.5, "d"), "not a probability",
  info = "a value above 1 is refused")
expect_equal(make_htest(0.25, "d")$p.value, 0.25,
  info = "an ordinary p-value passes through unchanged")

# Tables the kernels can handle must still be handled.
small <- matrix(c(30, 20, 20, 30), 2)
for (f in kern) {
  expect_equal(get(f)(small)$p.value, stats::fisher.test(small)$p.value,
    tolerance = 1e-10,
    info = paste(f, "still agrees with fisher.test on a moderate table"))
}

# net_dp_cpp() and net_vander_cpp() multiply binomial coefficients
# directly rather than working in logs, which is what makes the
# bulk-sum identity usable but caps the table size: choose(1020, 510)
# is 2.8e305 and choose(1030, 515) overflows a double. Past that the
# products became Inf, the ratios NaN, and a NaN came back as the
# p-value. The limit is now refused with its reason.
ok_size <- matrix(c(270, 230, 230, 270), 2)
for (f in c("net_dp_cpp", "net_vander_cpp")) {
  expect_equal(get(f)(ok_size)$p.value,
               stats::fisher.test(ok_size)$p.value, tolerance = 1e-10,
    info = paste(f, "is exact at a total of 1000"))
  expect_error(get(f)(matrix(rep(500, 4), 2)), "exceeds 1020",
    info = paste(f, "refuses a total beyond its overflow limit"))
}
# The kernels that work in logs carry on past it.
for (f in c("net_sweep_cpp", "tree_memo", "tree_memo_c", "tree_s4_c")) {
  expect_true(is.finite(get(f)(matrix(rep(6000, 4), 2))$p.value),
    info = paste(f, "handles a large table, working in logs"))
}
