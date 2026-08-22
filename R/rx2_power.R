# rx2_power.R
# Power and sample size for the r x 2 (2 x 2) Fisher's exact test.
#
# Implements the power function defined in Conlon & Thomas (1993),
# "Algorithm AS 280: The Power Function for Fisher's Exact Test",
# J. R. Stat. Soc. Series C 42(1):258-260, and the sample size search
# built on it in Thomas & Conlon (1992), "Sample size determination
# based on Fisher's exact test for use in 2 x 2 comparative trials
# with low event rates", Controlled Clinical Trials 13(2):134-147.
#
# This is a direct re-derivation of the power function from its
# mathematical definition (exhaustive enumeration of the rejection
# region against the product-binomial alternative), not a port of the
# original Fortran/NAG source, which is not available here. Given
# fixed margins n1, n2, the rejection region under the null is found
# from the conditional hypergeometric distribution of x1 given the
# total event count m = x1 + x2 (the same conditioning fisher.test()
# uses); power is the probability mass of that rejection region under
# independent Binomial(n1, p1) and Binomial(n2, p2) draws.
#
# UX follows stats::power.t.test()/power.prop.test(): leave n1 or
# power NULL and the function solves for the missing one, returning a
# classed "power.htest" object -- the same base-R idiom this package
# already borrows via make_htest()/"htest" for its p-value functions.

.rx2_exact_power <- function(n1, n2, p1, p2, alpha, alternative, tol) {
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

#' Power and sample size for Fisher's exact test on a 2 x 2 table
#'
#' Computes the exact power of Fisher's exact test comparing two
#' independent binomial proportions, or the sample size needed to
#' reach a target power, following the power function of Conlon &
#' Thomas (1993, Algorithm AS 280) and the sample-size method of
#' Thomas & Conlon (1992).
#'
#' Exactly one of \code{n1} and \code{power} must be \code{NULL}: leave
#' \code{power} \code{NULL} to compute the power at a given \code{n1}
#' (and \code{n2}), or leave \code{n1} \code{NULL} to find the smallest
#' \code{n1} reaching a target \code{power}. This mirrors the
#' \code{stats::power.t.test()} / \code{stats::power.prop.test()}
#' calling convention.
#'
#' @param n1 Sample size of group 1, or \code{NULL} to solve for it.
#' @param n2 Sample size of group 2. If \code{NULL}, taken as
#'   \code{ceiling(ratio * n1)} (with \code{n1} the value supplied or
#'   the value being solved for).
#' @param p1 True event probability in group 1.
#' @param p2 True event probability in group 2.
#' @param alpha Nominal significance level.
#' @param power Target power, or \code{NULL} to solve for it.
#' @param alternative One of \code{"two.sided"}, \code{"greater"}
#'   (rejects for an excess of events in group 1), or \code{"less"}.
#' @param ratio Allocation ratio \code{n2 / n1}, used only when
#'   \code{n2} is not supplied.
#' @param n1_max Largest group 1 size to search when solving for
#'   \code{n1}.
#' @param tol Tolerance used to decide which tables belong to the
#'   two-sided rejection region (tables whose null probability is
#'   within \code{tol} of the observed table's are included). Matches
#'   the tolerance \code{3.45254e-7} used throughout this package.
#' @return An object of class \code{"power.htest"}, with \code{n1},
#'   \code{n2}, \code{p1}, \code{p2}, \code{alpha}, \code{power}, and
#'   \code{alternative} elements filled in (whichever of \code{n1}/
#'   \code{n2} or \code{power} was solved for).
#' @examples
#' fisher_power(n1 = 40, n2 = 40, p1 = 0.3, p2 = 0.1)
#' fisher_power(p1 = 0.1, p2 = 0.3, power = 0.8)
#' @importFrom stats dbinom dhyper phyper
#' @export
fisher_power <- function(n1 = NULL, n2 = NULL, p1, p2, alpha = 0.05,
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
    power <- .rx2_exact_power(n1, n2, p1, p2, alpha, alternative, tol)
  } else {
    found <- FALSE
    for (cand_n1 in seq_len(n1_max)) {
      cand_n2 <- if (is.null(n2)) ceiling(ratio * cand_n1) else n2
      achieved <- .rx2_exact_power(cand_n1, cand_n2, p1, p2, alpha,
                                    alternative, tol)
      if (achieved >= power) {
        n1 <- cand_n1
        n2 <- cand_n2
        power <- achieved
        found <- TRUE
        break
      }
    }
    if (!found)
      stop(sprintf('target power not reached by n1_max = %d', n1_max),
           call. = FALSE)
  }

  structure(
    list(
      n1 = n1, n2 = n2, p1 = p1, p2 = p2, alpha = alpha, power = power,
      alternative = alternative,
      method = "Exact power for Fisher's exact test (Conlon & Thomas, 1993, Algorithm AS 280)"
    ),
    class = 'power.htest'
  )
}
