library(tinytest)

# fisher_power_fast() matches fisher_power() exactly (both optimizations
# are exact, not approximations) across direct power evaluation.
direct_cases <- list(
  list(n1 = 40,  n2 = 40,  p1 = 0.3,  p2 = 0.1,  alt = "less"),
  list(n1 = 40,  n2 = 40,  p1 = 0.3,  p2 = 0.1,  alt = "greater"),
  list(n1 = 40,  n2 = 40,  p1 = 0.3,  p2 = 0.1,  alt = "two.sided"),
  list(n1 = 30,  n2 = 60,  p1 = 0.05, p2 = 0.10, alt = "less"),
  list(n1 = 100, n2 = 100, p1 = 0.10, p2 = 0.20, alt = "greater"),
  list(n1 = 30,  n2 = 30,  p1 = 0.2,  p2 = 0.2,  alt = "two.sided")
)

for (cs in direct_cases) {
  slow <- fisher_power(cs$n1, cs$n2, p1 = cs$p1, p2 = cs$p2,
                        alternative = cs$alt)
  fast <- fisher_power_fast(cs$n1, cs$n2, p1 = cs$p1, p2 = cs$p2,
                             alternative = cs$alt)
  expect_identical(
    fast$power, slow$power,
    info = sprintf("n1=%d n2=%d p1=%.2f p2=%.2f alt=%s",
                   cs$n1, cs$n2, cs$p1, cs$p2, cs$alt)
  )
}

# Cross-checked against Exact::power.exact.test(method = "fisher"),
# same reference values used in test_power.R.
pw_2s <- fisher_power_fast(n1 = 40, n2 = 40, p1 = 0.3, p2 = 0.1,
                            alternative = "two.sided")
expect_equal(pw_2s$power, 0.5339897, tolerance = 1e-6)

pw_g <- fisher_power_fast(n1 = 40, n2 = 40, p1 = 0.3, p2 = 0.1,
                           alternative = "greater")
expect_equal(pw_g$power, 0.6554089, tolerance = 1e-6)

# n2 defaults to ceiling(ratio * n1) when not supplied, same as
# fisher_power().
pw_ratio <- fisher_power_fast(n1 = 40, p1 = 0.3, p2 = 0.1, ratio = 1.5)
expect_equal(pw_ratio$n2, 60)

# Solving for n1 (bisection) finds the same answer as fisher_power()'s
# linear scan, for both a default ratio = 1 search and an unequal
# allocation ratio, and for a fixed (non-derived) n2.
solve_cases <- list(
  list(p1 = 0.3, p2 = 0.1, power = 0.8, ratio = 1),
  list(p1 = 0.3, p2 = 0.1, power = 0.8, ratio = 2),
  list(p1 = 0.1, p2 = 0.3, power = 0.9, ratio = 1)
)
for (cs in solve_cases) {
  slow <- fisher_power(p1 = cs$p1, p2 = cs$p2, power = cs$power,
                        ratio = cs$ratio)
  fast <- fisher_power_fast(p1 = cs$p1, p2 = cs$p2, power = cs$power,
                             ratio = cs$ratio)
  expect_identical(
    fast$n1, slow$n1,
    info = sprintf("p1=%.2f p2=%.2f power=%.2f ratio=%.1f",
                   cs$p1, cs$p2, cs$power, cs$ratio)
  )
  expect_identical(fast$n2, slow$n2)
  expect_identical(fast$power, slow$power)
}

# Fixed n2 (not derived from ratio) while solving for n1: bisection
# still needs the correct n2 held constant across candidates.
slow_fixed_n2 <- fisher_power(n2 = 50, p1 = 0.3, p2 = 0.1, power = 0.8)
fast_fixed_n2 <- fisher_power_fast(n2 = 50, p1 = 0.3, p2 = 0.1, power = 0.8)
expect_identical(fast_fixed_n2$n1, slow_fixed_n2$n1)
expect_identical(fast_fixed_n2$n2, 50)
expect_identical(fast_fixed_n2$power, slow_fixed_n2$power)

# Argument validation matches fisher_power().
expect_error(fisher_power_fast(n1 = 10, n2 = 10, p1 = 0.1, p2 = 0.2, power = 0.8))
expect_error(fisher_power_fast(p1 = 0.1, p2 = 0.2))
expect_error(fisher_power_fast(n1 = 0, n2 = 10, p1 = 0.1, p2 = 0.2))
expect_error(fisher_power_fast(n1 = 10, n2 = 10, p1 = 1.5, p2 = 0.2))
expect_error(fisher_power_fast(n1 = 10, n2 = 10, p1 = 0.1, p2 = 0.2, alpha = 0))
expect_error(fisher_power_fast(p1 = 0.1, p2 = 0.2, power = 1.5))
expect_error(fisher_power_fast(p1 = 0.1, p2 = 0.2, power = 0.8, ratio = -1))
expect_error(fisher_power_fast(p1 = 0.99, p2 = 0.999, power = 0.999, n1_max = 5L))

# print method dispatches (inherited from "power.htest").
expect_stdout(
  print(fisher_power_fast(n1 = 40, n2 = 40, p1 = 0.3, p2 = 0.1)),
  "fast variant"
)

# Optimization 2 (doubling search + bisection) reaches the same
# answer as fisher_power()'s linear scan in a fraction of the time
# when the target n1 is large. n1 = 503 here; fisher_power_fast()
# should stay well under a second, while a plain linear scan (not
# exercised here to keep the test suite itself fast) took over 20s in
# manual verification.
system_time_fast <- system.time(
  large_n1 <- fisher_power_fast(p1 = 0.05, p2 = 0.10, power = 0.9,
                                 n1_max = 5000L, alternative = "less")
)[["elapsed"]]
expect_equal(large_n1$n1, 503)
expect_true(
  system_time_fast < 5,
  info = sprintf("fisher_power_fast() took %.2fs for a target n1 = 503 search",
                 system_time_fast)
)
