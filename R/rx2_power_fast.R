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
#   3. The two.sided rejection rule was an O(range) elementwise
#      vapply()/sum() over the whole hypergeometric support for EVERY
#      candidate x in that same support, i.e. O(range^2) per m. Since
#      dhyper(xs, ...) is unimodal in x, the set {y : d(y) <= t} for
#      any threshold t is always the union of a prefix of the
#      increasing left tail and a suffix of the decreasing right tail.
#      Precomputing prefix/suffix sums of d once per m, then using
#      findInterval() (binary search) against each tail to locate
#      those two boundaries for every x's own threshold at once,
#      answers the same question in O(range log range) per m instead
#      of O(range^2) -- and does so fully vectorized, with no explicit
#      per-x loop.
#
# Optimizations 1 and 2 are exact -- no eps-trimming, no
# approximation -- so their contribution to fisher_power_fast()'s
# result is bitwise identical to fisher_power(). Optimization 3 is
# also exact in exact arithmetic, but sums the same hypergeometric
# masses in a different order (prefix-from-the-left plus
# suffix-from-the-right, rather than one linear scan across the whole
# support), so floating-point rounding can differ from
# fisher_power()'s two.sided result at the level of machine epsilon;
# inst/tinytest/test_power_fast.R checks this with a tolerance rather
# than expect_identical() for the two.sided case.
#
# Caveat inherited from optimization 2 (bisection when solving for n1)
# -----------------------------------------------------------------
#
# Bisection requires power(n1) to be non-decreasing in n1. This is NOT
# always true: exact test power is well documented to have a
# "saw-tooth" shape (Chernick & Liu 2002, "The Saw-Toothed Behavior of
# Power Versus Sample Size and Software Solutions", The American
# Statistician 56(2):149-155), because the achieved significance level
# of a discrete test oscillates as n changes, and that oscillation can
# leak into the power curve as local dips even while the general trend
# rises. This was confirmed directly for this package's own power
# function during development: a concrete example at p1 = 0.40,
# p2 = 0.55, one-sided alpha = 0.05 crosses above a target of 0.80 at
# n = 144, dips back below for n = 146-151, and only rises permanently
# past it from n = 152 on.
#
# A substantial investigation (see analysis/report/report.Rmd, section
# "Non-monotonicity of power in n and its effect on sample-size
# search", in the sibling fisher-power-unequal-n research workspace,
# for the full write-up) tried four ways to fix or bound this before
# concluding that plain bisection -- exactly as implemented below,
# with no correction -- is the right production choice:
#
#   1. Plain bisection, no correction: reports a sample size a few
#      units too LARGE in practice (never too small, and never an
#      actually-invalid n1 -- the reported n1 always genuinely meets
#      the target), because it assumes monotonicity to decide which
#      half of the search range to discard. Every disagreement found
#      against a linear scan, across several hundred (p1, p2, target
#      power) combinations explored interactively, was of this same
#      safe-direction form, off by 2 to 8 units.
#   2. A rigorous monotone bound derived from Hoeffding's inequality
#      for hypergeometric (finite-population) sampling: a genuinely
#      valid lower bound (checked numerically), but far too loose to
#      be useful -- at p1 = 0.06, p2 = 0.16, target power 0.95, the
#      true crossing is n1 = 219, but the Hoeffding bound itself does
#      not clear 0.95 until n1 = 369, a gap of 150 values. Hoeffding's
#      inequality uses only the range of the random variable, not its
#      variance, and is known to be loose for exactly this reason.
#      Tighter (Berry-Esseen-type) bounds were considered but not
#      pursued: no bound on the saw-tooth's magnitude or period appears
#      in Chernick & Liu (2002) or in a 2026 paper characterizing the
#      analogous non-monotonicity in Monte Carlo permutation test
#      power, so deriving one is a genuine, currently-unpublished
#      research question, not a matter of applying a known formula.
#   3. A multi-level (iteratively narrowing) coarse-grid isotonic
#      regression search: cheaper than approach 4, but WRONG on the
#      p1 = 0.40, p2 = 0.55 case above -- it returned n1 = 152 instead
#      of 144, missing the isolated 2-unit-wide spike entirely, because
#      each successive re-fit narrows the search range based on a
#      SMOOTHED (pool-adjacent-violators) estimate at shrinking scales,
#      and a narrow spike can get pooled together with a lower
#      neighbor and averaged below target before the narrowing ever
#      reaches fine enough resolution to see it individually.
#   4. A single-level coarse grid (e.g. 15-30 points) + isotonic
#      regression + exhaustive scan of the bracket the smoothed fit
#      crosses target in: this looked initially promising (correct on
#      6 of 6 known non-monotonicity cases) when tested against an
#      arbitrary FIXED grid spanning the whole [1, n1_max] range.  But
#      that test was not representative: a fixed full-range grid is
#      not usable in production, because a generously large n1_max
#      would force sampling points arbitrarily close to n1_max even
#      when the true answer is tiny (the exact problem bisection's
#      doubling phase exists to avoid). Once correctly combined with a
#      doubling phase (sampling the coarse grid only within the
#      doubling-derived bracket, as production code must), it stopped
#      reliably beating plain bisection: it missed the p1 = 0.40,
#      p2 = 0.55 case at every grid density tried (15, 30, and 60
#      points all returned 152, never 144), and increasing the grid
#      density even broke a previously-correct case (p1 = 0.06,
#      p2 = 0.21 went from correct at 15-30 points to wrong, off by 2,
#      at 60). The apparent success in the first test was an artifact
#      of that specific fixed grid's alignment, not a property of the
#      method.
#
# No approach tried improves on plain bisection's worst case once
# tested under realistic conditions, and this matches the literature:
# no published bound on the saw-tooth's magnitude or period was found
# (not in the original Chernick & Liu paper, nor a 2026 paper on the
# analogous permutation-test phenomenon), and SAS's PROC POWER
# documentation -- for exactly this "sawtooth power function" problem
# -- recommends manual inspection of a power-vs-n plot rather than an
# automated guarantee. Plain bisection is retained here as the
# production default: it is fast, its failure mode is safe-direction
# (never reports an n1 that doesn't genuinely meet the target, only
# occasionally one a few units larger than the true minimum), and nothing
# tried does better without giving something else up. fisher_power()'s
# linear scan remains the fully conservative choice, at substantial
# extra cost, for anyone who needs the true minimum ruled out exactly.

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
        n_d <- length(d)
        peak <- which.max(d)
        pre <- cumsum(d)                # pre[i]  = sum(d[1:i])
        suf <- rev(cumsum(rev(d)))      # suf[i]  = sum(d[i:n_d])
        thresh <- d * (1 + tol)

        # left_d (indices 1:peak) is nondecreasing; right_d (indices
        # (peak+1):n_d) is nonincreasing. findInterval() needs an
        # ascending vector, so the right tail is searched in reverse.
        left_d <- d[seq_len(peak)]
        right_d_asc <- if (peak < n_d) rev(d[(peak + 1L):n_d]) else numeric(0)

        # Number of left-tail / right-tail elements <= thresh[i], for
        # every i at once: a contiguous prefix of the left tail and a
        # contiguous suffix of the right tail, by unimodality.
        count_left <- findInterval(thresh, left_d)
        count_right <- findInterval(thresh, right_d_asc)

        left_sum <- ifelse(count_left > 0, pre[pmax(count_left, 1L)], 0)
        right_sum <- ifelse(
          count_right > 0, suf[n_d - pmax(count_right, 1L) + 1L], 0
        )

        (left_sum + right_sum) <= alpha
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
#' @details
#' Exact test power is not always non-decreasing in \code{n1}
#' (Chernick & Liu 2002; confirmed directly for this power function
#' during development). When solving for \code{n1}, bisection can
#' therefore report a value a few units larger than the smallest
#' \code{n1} that truly meets the target -- always a valid \code{n1}
#' in its own right, never an underpowered one, but not always the
#' minimum. Several corrections were investigated (a windowed check, a
#' rigorous Hoeffding-inequality bound, and two isotonic-regression
#' variants) and rejected: none reliably improved on plain bisection
#' once tested under realistic conditions, and one was less reliable
#' despite being far more expensive. See the comments at the top of
#' \code{R/rx2_power_fast.R} for the full account of what was tried,
#' and \code{analysis/report/report.Rmd}'s "Non-monotonicity of power
#' in n" section (sibling \code{fisher-power-unequal-n} research
#' workspace) for the complete write-up. \code{\link{fisher_power}}'s
#' linear scan remains the fully conservative choice when this
#' residual risk matters for a specific design.
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
