library(tinytest)

# ss_fxpower_fast() (doubling search + bisection) reaches the same
# answer as ss_fxpower()'s linear scan, across a range of target
# sizes spanning small (~hundreds) to large (~thousands) n.
solve_cases <- list(
  list(p1 = 0.05, p2 = 0.10, target_power = 0.9, k = 1),
  list(p1 = 0.05, p2 = 0.20, target_power = 0.8, k = 1),
  list(p1 = 0.10, p2 = 0.30, target_power = 0.8, k = 2),
  list(p1 = 0.30, p2 = 0.10, target_power = 0.8, k = 1)
)
for (cs in solve_cases) {
  slow <- ss_fxpower(p1 = cs$p1, p2 = cs$p2, target_power = cs$target_power,
                     k = cs$k, eps = 1e-6, n_max = 20000L)
  fast <- ss_fxpower_fast(p1 = cs$p1, p2 = cs$p2, target_power = cs$target_power,
                          k = cs$k, eps = 1e-6, n_max = 20000L)
  info <- sprintf("p1=%.2f p2=%.2f power=%.1f k=%.1f",
                  cs$p1, cs$p2, cs$target_power, cs$k)
  expect_identical(fast$n1, slow$n1, info = info)
  expect_identical(fast$n2, slow$n2, info = info)
  expect_identical(fast$power, slow$power, info = info)
  expect_identical(fast$n_total, slow$n_total, info = info)
}

# two.sided alternative also agrees between the two solvers.
slow_2s <- ss_fxpower(p1 = 0.05, p2 = 0.20, target_power = 0.8, eps = 1e-6)
# ss_fxpower() has no alternative argument (always "less"); confirm
# ss_fxpower_fast()'s default matches it, then check two.sided works.
fast_default <- ss_fxpower_fast(p1 = 0.05, p2 = 0.20, target_power = 0.8, eps = 1e-6)
expect_identical(fast_default$n1, slow_2s$n1)
expect_identical(fast_default$n2, slow_2s$n2)

fast_2s <- ss_fxpower_fast(p1 = 0.05, p2 = 0.20, target_power = 0.8, eps = 1e-6,
                           alternative = "two.sided")
expect_true(fast_2s$power >= 0.8)
expect_identical(fast_2s$alternative, "two.sided")

# n_max reached without success: both solvers warn rather than error,
# and the fast path still reports its best attempt at n_max.
expect_warning(
  ss_fxpower_fast(p1 = 0.5, p2 = 0.501, target_power = 0.999, n_max = 20L)
)

# Speed: the doubling search reaches a target requiring n ~ 500-9000
# in a small fraction of the linear scan's time. Not a strict timing
# assertion (flaky on shared/loaded machines) -- just confirms the
# fast path completes quickly on its own.
t_fast <- system.time(
  big <- ss_fxpower_fast(p1 = 0.05, p2 = 0.06, target_power = 0.9,
                         eps = 1e-6, n_max = 20000L)
)[["elapsed"]]
expect_equal(big$n1, 8916)
expect_true(
  t_fast < 2,
  info = sprintf("ss_fxpower_fast() took %.2fs for a target n ~ 8916 search",
                 t_fast)
)
