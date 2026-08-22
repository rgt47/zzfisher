library(tinytest)

# fxpower(alternative = "two.sided") matches the pure-R reference
# (fisherexacttestrx2's fisher_power_fast(alternative = "two.sided"))
# to within eps-trimming tolerance. Note: this test suite cannot
# depend on that sibling package, so these values are hardcoded from a
# one-time cross-check (see fxpower_rcpp.cpp's cv_two_sided()
# comments for the algorithm and inst/benchmarks or the session that
# added this feature for the cross-validation).
two_sided_cases <- list(
  list(n1 = 40,  n2 = 40,  p1 = 0.3,  p2 = 0.1,  expected = 0.5339896594),
  list(n1 = 200, n2 = 100, p1 = 0.05, p2 = 0.30, expected = 0.9999025042),
  list(n1 = 5,   n2 = 5,   p1 = 0.01, p2 = 0.01, expected = 9.43478189101799e-08),
  list(n1 = 5,   n2 = 5,   p1 = 0.99, p2 = 0.99, expected = 9.43478189101803e-08)
)
for (cs in two_sided_cases) {
  res <- fxpower(cs$n1, cs$n2, cs$p1, cs$p2, eps = 0, alternative = "two.sided")
  expect_equal(
    res$power, cs$expected, tolerance = 1e-8,
    info = sprintf("n1=%d n2=%d p1=%.2f p2=%.2f", cs$n1, cs$n2, cs$p1, cs$p2)
  )
}

# alternative = "less" (the default) is completely unaffected: same
# results as calling fxpower() with no alternative argument at all.
default_res <- fxpower(100, 100, 0.05, 0.10, eps = 0)
less_res <- fxpower(100, 100, 0.05, 0.10, eps = 0, alternative = "less")
expect_identical(default_res$power, less_res$power)
expect_identical(default_res$alternative, "less")

# Regression test for a bug found and fixed while developing this
# feature: the two-sided path's eps-trimming walked each rejection
# block from its own OUTER boundary rather than from the alternative
# distribution's own mode (when the mode falls inside that block).
# Since the alternative-distribution density is itself unimodal, this
# made the walk's f values increase rather than decrease, tripping the
# f < flim early stop immediately and silently discarding whole blocks
# -- collapsing power to ~0 at large n where the trimming inevitably
# engages. n = 5000 triggered it; smaller n happened not to.
big_exact <- fxpower(5000, 5000, 0.05, 0.10, eps = 0, alternative = "two.sided")
big_trim  <- fxpower(5000, 5000, 0.05, 0.10, eps = 1e-6, alternative = "two.sided")
expect_true(
  big_exact$power > 0.99,
  info = sprintf("n=5000 two.sided exact power was %.6f, expected > 0.99",
                 big_exact$power)
)
expect_true(
  big_trim$power > 0.99,
  info = sprintf("n=5000 two.sided trimmed power was %.6f, expected > 0.99",
                 big_trim$power)
)

# The trimming error bound is a genuine upper bound on the actual
# discrepancy from the exact (eps = 0) result, across a grid including
# boundary-peak cases (extreme p1/p2 forcing the alternative's mode to
# either end of the support).
bound_cases <- list(
  list(n1 = 40,  n2 = 40,  p1 = 0.3,  p2 = 0.1),
  list(n1 = 30,  n2 = 60,  p1 = 0.05, p2 = 0.10),
  list(n1 = 200, n2 = 100, p1 = 0.05, p2 = 0.30),
  list(n1 = 5,   n2 = 5,   p1 = 0.01, p2 = 0.01),
  list(n1 = 5,   n2 = 5,   p1 = 0.99, p2 = 0.99),
  list(n1 = 500, n2 = 500, p1 = 0.05, p2 = 0.10),
  list(n1 = 1000, n2 = 500, p1 = 0.05, p2 = 0.10)
)
for (cs in bound_cases) {
  ex <- fxpower(cs$n1, cs$n2, cs$p1, cs$p2, eps = 0, alternative = "two.sided")
  tr <- fxpower(cs$n1, cs$n2, cs$p1, cs$p2, eps = 1e-6, alternative = "two.sided")
  expect_true(
    abs(ex$power - tr$power) <= tr$error + 1e-12,
    info = sprintf(
      "n1=%d n2=%d p1=%.2f p2=%.2f: actual diff %.3e exceeds claimed bound %.3e",
      cs$n1, cs$n2, cs$p1, cs$p2, abs(ex$power - tr$power), tr$error
    )
  )
}

# print.fxpower reports the correct sidedness.
expect_stdout(
  print(fxpower(40, 40, 0.3, 0.1, alternative = "two.sided")),
  "two-sided"
)
expect_stdout(
  print(fxpower(40, 40, 0.3, 0.1)),
  "one-sided"
)

# Invalid alternative is rejected.
expect_error(fxpower(40, 40, 0.3, 0.1, alternative = "greater"))
