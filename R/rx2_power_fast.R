# rx2_power_fast.R
# A faster variant of fisher_power() (see R/rx2_power.R), incorporating
# two independent, exactness-preserving optimizations identified when
# comparing fisher_power() against fxpower()'s compiled, mode-centered
# kernel in zzfisher's benchmarking framework
# (inst/benchmarks/benchmark_power_comparison.R):
#
#   1. dhyper(xs, n1, n2, m) is only needed by the two.sided rejection
#      rule. fisher_power()'s inner loop computed it unconditionally
#      on every iteration of the outer m loop, wasting a full
#      vectorized hypergeometric call for the one-sided cases -- the
#      common case, and the only case fxpower() computes.
#   2. Solving for n1 (fisher_power(power = ...)) was a linear scan
#      over candidate n1 = 1, 2, 3, ... Power is monotonically
#      increasing in n1 at fixed p1, p2, alpha, ratio (asserted in
#      fisher_power()'s own test suite), so a bisection search reaches
#      the answer in O(log n1_max) power evaluations instead of
#      O(n1_max).
#
# Both changes are exact -- no eps-trimming, no approximation -- so
# fisher_power_fast() returns bitwise-identical power values to
# fisher_power() at every (n1, n2, p1, p2, alpha, alternative); only
# the n1-solving path length differs. Verified in
# inst/tinytest/test_power_fast.R.
#
# Caveat inherited from optimization 2: bisection requires power(n1)
# to be non-decreasing in n1. This holds for the exact test in every
# case checked so far, but discreteness can in principle produce a
# non-monotonic wiggle that a linear scan (which just takes the first
# n1 crossing the target, scanning in order) would not be tripped up
# by. If that mattered for a given design, fisher_power()'s linear
# scan remains the conservative choice.

.rx2_exact_power_fast <- function(n1, n2, p1, p2, alpha, alternative, tol) {
  b1 <- dbinom(0:n1, n1, p1)
  b2 <- dbinom(0:n2, n2, p2)

  power <- 0
  for (m in 0:(n1 + n2)) {
    lo <- max(0L, m - n2)
    hi <- min(n1, m)
    if (lo > hi) next
    xs <- lo:hi

    reject <- switch(
      alternative,
      greater = phyper(xs - 1, n1, n2, m, lower.tail = FALSE) <= alpha,
      less    = phyper(xs, n1, n2, m, lower.tail = TRUE) <= alpha,
      two.sided = {
        d <- dhyper(xs, n1, n2, m)
        vapply(
          seq_along(xs),
          function(i) sum(d[d <= d[i] * (1 + tol)]) <= alpha,
          logical(1)
        )
      }
    )

    if (any(reject))
      power <- power + sum(b1[xs + 1L] * b2[m - xs + 1L] * reject)
  }
  power
}

#' Power and sample size for Fisher's exact test (fast variant)
#'
#' A faster variant of \code{\link{fisher_power}} with the same exact
#' calling convention and numerically identical results, differing
#' only in cost: the \code{dhyper()} call wasted on one-sided
#' alternatives is removed, and solving for \code{n1} uses bisection
#' (which assumes power is non-decreasing in \code{n1}) instead of a
#' linear scan from 1. Added as a competing pure-R implementation
#' alongside \code{\link{fisher_power}} in zzfisher's power-function
#' benchmarking framework.
#'
#' @inheritParams fisher_power
#' @return An object of class \code{"power.htest"}, identical in
#'   structure to \code{\link{fisher_power}}'s return value.
#' @examples
#' fisher_power_fast(n1 = 40, n2 = 40, p1 = 0.3, p2 = 0.1)
#' fisher_power_fast(p1 = 0.1, p2 = 0.3, power = 0.8)
#' @export
fisher_power_fast <- function(n1 = NULL, n2 = NULL, p1, p2, alpha = 0.05,
                               power = NULL,
                               alternative = c('two.sided', 'greater', 'less'),
                               ratio = 1, n1_max = 1000L, tol = 3.45254e-7) {
  alternative <- match.arg(alternative)
  if (sum(is.null(n1), is.null(power)) != 1)
    stop("exactly one of 'n1' and 'power' must be NULL", call. = FALSE)
  if (p1 < 0 || p1 > 1 || p2 < 0 || p2 > 1)
    stop('p1 and p2 must be in [0, 1]', call. = FALSE)
  if (alpha <= 0 || alpha >= 1)
    stop('alpha must be in (0, 1)', call. = FALSE)
  if (!is.null(power) && (power <= 0 || power >= 1))
    stop('power must be in (0, 1)', call. = FALSE)
  if (ratio <= 0)
    stop('ratio must be positive', call. = FALSE)

  if (is.null(power)) {
    if (n1 < 1 || n1 != round(n1))
      stop('n1 must be a positive integer', call. = FALSE)
    if (is.null(n2)) {
      n2 <- ceiling(ratio * n1)
    } else if (n2 < 1 || n2 != round(n2)) {
      stop('n2 must be a positive integer', call. = FALSE)
    }
    power <- .rx2_exact_power_fast(n1, n2, p1, p2, alpha, alternative, tol)
  } else {
    n2_of <- function(cand_n1) if (is.null(n2)) ceiling(ratio * cand_n1) else n2
    power_of <- function(cand_n1) {
      .rx2_exact_power_fast(cand_n1, n2_of(cand_n1), p1, p2, alpha,
                             alternative, tol)
    }
    n1_max_int <- as.integer(n1_max)

    # Exponential (doubling) search for a bracket, THEN bisect inside
    # it. Each power_of() call costs O(n1 * n2), so a plain bisection
    # against n1_max would pay for one evaluation near n1_max on its
    # very first step even when the true answer is far smaller --
    # trading fewer evaluations for one arbitrarily expensive one.
    # Doubling keeps every evaluation close to the true answer's own
    # scale: lo starts at 0 (power_of(0) is treated as below target
    # without evaluating it) and hi doubles from 1 until it brackets
    # the target or reaches n1_max_int.
    lo <- 0L
    hi <- 1L
    repeat {
      if (hi >= n1_max_int) {
        hi <- n1_max_int
        if (power_of(hi) < power)
          stop(sprintf('target power not reached by n1_max = %d', n1_max),
               call. = FALSE)
        break
      }
      if (power_of(hi) >= power) break
      lo <- hi
      hi <- hi * 2L
    }

    while (lo + 1L < hi) {
      mid <- lo + (hi - lo) %/% 2L
      if (power_of(mid) >= power) {
        hi <- mid
      } else {
        lo <- mid
      }
    }

    n1 <- hi
    n2 <- n2_of(n1)
    power <- power_of(n1)
  }

  structure(
    list(
      n1 = n1, n2 = n2, p1 = p1, p2 = p2, alpha = alpha, power = power,
      alternative = alternative,
      method = paste("Exact power for Fisher's exact test",
                     "(Conlon & Thomas, 1993, Algorithm AS 280)",
                     "-- fast variant")
    ),
    class = 'power.htest'
  )
}
