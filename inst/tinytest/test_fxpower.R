library(tinytest)

# fxpower matches pure R reference for equal groups
cases <- list(
    list(n = 100, p1 = 0.05, p2 = 0.10, expected = 0.2781258),
    list(n = 200, p1 = 0.05, p2 = 0.10, expected = 0.5369213),
    list(n = 300, p1 = 0.05, p2 = 0.10, expected = 0.7107743),
    list(n = 100, p1 = 0.10, p2 = 0.20, expected = 0.5672980),
    list(n = 200, p1 = 0.10, p2 = 0.20, expected = 0.8541454),
    list(n = 100, p1 = 0.20, p2 = 0.30, expected = 0.4368245),
    list(n = 200, p1 = 0.20, p2 = 0.30, expected = 0.7135070)
)

for (tc in cases) {
    res <- fxpower(tc$n, tc$n, tc$p1, tc$p2, alpha = 0.05, eps = 0)
    expect_equal(
      res$power, tc$expected,
      tolerance = 1e-5,
      label = sprintf("n=%d, p1=%.2f, p2=%.2f", tc$n, tc$p1, tc$p2)
    )
}


# C++ matches pure R reference for unequal groups
params <- list(
    list(n1 = 50, n2 = 50, p1 = 0.10, p2 = 0.30),
    list(n1 = 100, n2 = 100, p1 = 0.05, p2 = 0.10),
    list(n1 = 200, n2 = 100, p1 = 0.05, p2 = 0.10),
    list(n1 = 100, n2 = 200, p1 = 0.10, p2 = 0.20),
    list(n1 = 150, n2 = 100, p1 = 0.20, p2 = 0.30),
    list(n1 = 50, n2 = 100, p1 = 0.05, p2 = 0.15)
)

for (p in params) {
    cpp <- fxpower(p$n1, p$n2, p$p1, p$p2, alpha = 0.05, eps = 0)
    r_ref <- fxpower_r(p$n1, p$n2, p$p1, p$p2, alpha = 0.05)
    expect_equal(
      cpp$power, r_ref,
      tolerance = 1e-10,
      label = sprintf("n1=%d, n2=%d, p1=%.2f, p2=%.2f",
                      p$n1, p$n2, p$p1, p$p2)
    )
}


# fxpower's one-sided test agrees with fisher_power(alternative = "less"):
# both reject for an unusually small x1 (few events in group 1), which is
# the direction fxpower's cv_unequal() accumulates from xmin_r upward.
for (p in params) {
  fx <- fxpower(p$n1, p$n2, p$p1, p$p2, alpha = 0.05, eps = 0)
  fp <- fisher_power(p$n1, p$n2, p1 = p$p1, p2 = p$p2, alpha = 0.05,
                      alternative = "less")
  expect_equal(
    fx$power, fp$power,
    tolerance = 1e-10,
    label = sprintf("fxpower vs fisher_power: n1=%d, n2=%d, p1=%.2f, p2=%.2f",
                    p$n1, p$n2, p$p1, p$p2)
  )
}


# more total subjects gives more power, at a fixed allocation ratio
# (n1 = 200, n2 = 100 has total N = 300 vs. N = 200 for n1 = n2 = 100)
equal <- fxpower(100, 100, 0.05, 0.10, alpha = 0.05)
larger_unequal <- fxpower(200, 100, 0.05, 0.10, alpha = 0.05)
expect_true(larger_unequal$power > equal$power)


# balanced allocation is more efficient than unequal allocation at
# a fixed total N
balanced <- fxpower(150, 150, 0.05, 0.10, alpha = 0.05)
unequal_same_n <- fxpower(200, 100, 0.05, 0.10, alpha = 0.05)
expect_true(
  balanced$power > unequal_same_n$power,
  info = "balanced n1=n2=150 (N=300) should beat n1=200,n2=100 (N=300)"
)


# power increases with sample size
p1 <- fxpower(100, 100, 0.10, 0.20, alpha = 0.05)
p2 <- fxpower(200, 200, 0.10, 0.20, alpha = 0.05)
p3 <- fxpower(300, 300, 0.10, 0.20, alpha = 0.05)
expect_true(p1$power < p2$power)
expect_true(p2$power < p3$power)


# alpha_star is conservative (leq nominal alpha)
res <- fxpower(100, 100, 0.05, 0.10, alpha = 0.05)
expect_true(res$alpha_star <= 0.05)


# fxpower rejects invalid inputs
expect_error(fxpower(0, 100, 0.05, 0.10))
expect_error(fxpower(100, 100, 0, 0.10))
expect_error(fxpower(100, 100, 0.05, 1.0))
expect_error(fxpower(100, 100, 0.05, 0.10, alpha = 0))
expect_error(fxpower(100, 100, 0.05, 0.10, eps = -1))


# eps trimming gives close approximation
exact <- fxpower(100, 100, 0.10, 0.20, alpha = 0.05, eps = 0)
approx <- fxpower(100, 100, 0.10, 0.20, alpha = 0.05, eps = 1e-6)
expect_equal(exact$power, approx$power, tolerance = 1e-4)
expect_true(approx$count < exact$count)


# error field: 0 when eps = 0 (no trimming), and a genuine upper
# bound on the trimming-induced power discrepancy when eps > 0
expect_equal(
  exact$error, 0,
  info = "error must be exactly 0 when eps = 0 (no trimming)"
)
expect_true(
  approx$error >= abs(approx$power - exact$power),
  info = "error must upper-bound the actual trimming discrepancy"
)
expect_true(
  approx$error > 0 && approx$error < 1,
  info = "error must be a nondegenerate bound when eps > 0"
)


# ss_fxpower finds correct sample size
ss <- ss_fxpower(p1 = 0.10, p2 = 0.20, alpha = 0.05,
                   target_power = 0.80, k = 1)
expect_true(ss$power >= 0.80)
below <- fxpower(ss$n1 - 1L, ss$n2 - 1L, 0.10, 0.20, alpha = 0.05)
expect_true(below$power < 0.80)


# component functions return valid values
cv <- zzfisher:::.cv_unequal_cpp(100L, 100L, 10L, 0.05)
expect_true(is.integer(cv) || is.numeric(cv))
expect_true(cv >= 0 || cv == -1L)

xfm <- zzfisher:::.xfmax_unequal_cpp(100L, 100L, 10L, 1.0)
expect_true(xfm >= 0 && xfm <= 10)

bj <- zzfisher:::.binom_joint_cpp(100L, 100L, 0.10, 0.20, 5L, 15L)
expect_true(bj > 0 && bj < 1)


# power_table returns correct structure
tbl <- power_table(
    n2_vec = c(50, 100),
    k_vec = c(1, 2),
    p1 = 0.05, p2 = 0.10
)
expect_inherits(tbl, "data.frame")
expect_equal(nrow(tbl), 4)
expect_true(all(c("n1", "n2", "k", "ntotal", "power",
                     "alpha_star") %in% names(tbl)))
expect_true(all(tbl$power > 0 & tbl$power < 1))


# print.fxpower produces output
res <- fxpower(100, 100, 0.05, 0.10)
expect_stdout(
  print(res), "Power of Fisher's Exact Test",
  info = "print.fxpower must report the test being described"
)
