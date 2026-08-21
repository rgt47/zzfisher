library(tinytest)

# Cross-checked against Exact::power.exact.test(method = "fisher"),
# n1 = n2 = 40, p1 = 0.3, p2 = 0.1.
pw_2s <- fisher_power(n1 = 40, n2 = 40, p1 = 0.3, p2 = 0.1,
                       alternative = "two.sided")
expect_inherits(pw_2s, "power.htest")
expect_equal(pw_2s$power, 0.5339897, tolerance = 1e-6)
expect_equal(pw_2s$n1, 40)
expect_equal(pw_2s$n2, 40)
expect_identical(pw_2s$alternative, "two.sided")

pw_g <- fisher_power(n1 = 40, n2 = 40, p1 = 0.3, p2 = 0.1,
                      alternative = "greater")
expect_equal(pw_g$power, 0.6554089, tolerance = 1e-6)

# n2 defaults to ceiling(ratio * n1) when not supplied.
pw_ratio <- fisher_power(n1 = 40, p1 = 0.3, p2 = 0.1, ratio = 1.5)
expect_equal(pw_ratio$n2, 60)

# Power increases monotonically with sample size (fixed effect, alpha).
pw_small <- fisher_power(n1 = 10, n2 = 10, p1 = 0.3, p2 = 0.1)
pw_large <- fisher_power(n1 = 100, n2 = 100, p1 = 0.3, p2 = 0.1)
expect_true(pw_large$power > pw_small$power)

# Power under the null (p1 == p2) is at most alpha (conservative test).
pw_null <- fisher_power(n1 = 30, n2 = 30, p1 = 0.2, p2 = 0.2, alpha = 0.05)
expect_true(pw_null$power <= 0.05)

# Argument validation.
expect_error(fisher_power(n1 = 10, n2 = 10, p1 = 0.1, p2 = 0.2, power = 0.8))
expect_error(fisher_power(p1 = 0.1, p2 = 0.2))
expect_error(fisher_power(n1 = 0, n2 = 10, p1 = 0.1, p2 = 0.2))
expect_error(fisher_power(n1 = 10, n2 = 10, p1 = 1.5, p2 = 0.2))
expect_error(fisher_power(n1 = 10, n2 = 10, p1 = 0.1, p2 = 0.2, alpha = 0))
expect_error(fisher_power(p1 = 0.1, p2 = 0.2, power = 1.5))
expect_error(fisher_power(p1 = 0.1, p2 = 0.2, power = 0.8, ratio = -1))

# Solving for n1 finds the smallest size reaching the target power, and
# the achieved power at that size matches a direct fisher_power() call
# with n1/n2 supplied.
sz <- fisher_power(p1 = 0.3, p2 = 0.1, power = 0.8, alternative = "two.sided")
expect_inherits(sz, "power.htest")
expect_true(sz$power >= 0.8)
expect_equal(
  fisher_power(sz$n1, sz$n2, p1 = 0.3, p2 = 0.1,
               alternative = "two.sided")$power,
  sz$power
)
expect_true(
  fisher_power(sz$n1 - 1, ceiling(sz$n1 - 1), p1 = 0.3, p2 = 0.1,
               alternative = "two.sided")$power < 0.8
)
