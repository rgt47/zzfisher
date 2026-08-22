#' Power of Fisher's Exact Test with Unequal Group Sizes
#'
#' Compute the exact power of Fisher's conditional test for a
#' 2x2 table with independent binomial samples of sizes n1 and n2.
#' Generalizes Algorithm AS 280 (Conlon & Thomas, 1993) to
#' arbitrary allocation ratios.
#'
#' @param n1 Integer. Sample size in group 1 (treatment).
#' @param n2 Integer. Sample size in group 2 (control).
#' @param p1 Numeric. True event probability in group 1.
#' @param p2 Numeric. True event probability in group 2.
#' @param alpha Numeric. Type I error rate (one-sided if
#'   \code{alternative = "less"}, two-sided minimum-likelihood
#'   convention if \code{alternative = "two.sided"}).
#' @param eps Numeric. Trimming parameter controlling speed vs
#'   accuracy. Set to 0 for exact power (slower). Values around
#'   1e-5 to 1e-7 give excellent approximations with bounded error.
#' @param alternative One of \code{"less"} (the original one-sided
#'   test: rejects for an unusually small x1) or \code{"two.sided"}
#'   (Fisher's minimum-likelihood rule, matching
#'   \code{fisher_power()}'s \code{alternative = "two.sided"} and
#'   \code{fisher_power_fast()}'s -- see \code{tol} there for the
#'   tie tolerance, fixed here at \code{3.45254e-7}).
#'
#' @return A list of class \code{"fxpower"} with components:
#'   \item{power}{Computed power of the test.}
#'   \item{error}{Upper bound on the omitted probability mass due
#'     to eps-trimming: the smaller of (1 - power) and (flim times
#'     the number of support points not evaluated), 0 if eps = 0.
#'     Follows the ERROR calculation in the original Algorithm
#'     AS 280 FORTRAN (see docs/fxpower_original.f).}
#'   \item{count}{Number of sample space points evaluated.}
#'   \item{n1}{Group 1 sample size (echoed).}
#'   \item{n2}{Group 2 sample size (echoed).}
#'   \item{p1}{Group 1 event probability (echoed).}
#'   \item{p2}{Group 2 event probability (echoed).}
#'   \item{alpha}{Type I error rate (echoed).}
#'   \item{alpha_star}{Empirical type I error rate under H0.}
#'
#' @examples
#' fxpower(n1 = 100, n2 = 100, p1 = 0.05, p2 = 0.10, alpha = 0.05)
#' fxpower(n1 = 200, n2 = 100, p1 = 0.05, p2 = 0.10, alpha = 0.05)
#'
#' @export
fxpower <- function(n1, n2, p1, p2, alpha = 0.05, eps = 1e-6,
                    alternative = c("less", "two.sided")) {
  alternative <- match.arg(alternative)
  stopifnot(
    is.numeric(n1), length(n1) == 1, n1 > 0, n1 == as.integer(n1),
    is.numeric(n2), length(n2) == 1, n2 > 0, n2 == as.integer(n2),
    is.numeric(p1), length(p1) == 1, p1 > 0, p1 < 1,
    is.numeric(p2), length(p2) == 1, p2 > 0, p2 < 1,
    is.numeric(alpha), length(alpha) == 1, alpha > 0, alpha < 1,
    is.numeric(eps), length(eps) == 1, eps >= 0, eps < 1
  )

  two_sided <- identical(alternative, "two.sided")

  result <- .fxpower_cpp(
    as.integer(n1), as.integer(n2),
    as.double(p1), as.double(p2),
    as.double(alpha), as.double(eps),
    two_sided
  )

  alpha_star <- fxpower_alpha_star(n1, n2, p2, alpha, eps, alternative)

  structure(
    list(
      power = result$power,
      error = result$error,
      count = result$count,
      n1 = n1,
      n2 = n2,
      p1 = p1,
      p2 = p2,
      alpha = alpha,
      alpha_star = alpha_star,
      alternative = alternative
    ),
    class = "fxpower"
  )
}


#' Empirical Type I Error for Fisher's Exact Test
#'
#' Compute alpha* by evaluating the power function under H0
#' (p1 = p2 = p0).
#'
#' @param n1 Integer. Sample size in group 1.
#' @param n2 Integer. Sample size in group 2.
#' @param p0 Numeric. Common event probability under H0.
#' @param alpha Numeric. Nominal significance level.
#' @param eps Numeric. Trimming parameter.
#' @param alternative One of \code{"less"} or \code{"two.sided"};
#'   see \code{\link{fxpower}}.
#'
#' @return Numeric scalar: the empirical type I error rate.
#'
#' @keywords internal
fxpower_alpha_star <- function(n1, n2, p0, alpha, eps = 1e-6,
                               alternative = c("less", "two.sided")) {
  alternative <- match.arg(alternative)
  result <- .fxpower_cpp(
    as.integer(n1), as.integer(n2),
    as.double(p0), as.double(p0),
    as.double(alpha), as.double(eps),
    identical(alternative, "two.sided")
  )
  result$power
}


#' @export
print.fxpower <- function(x, ...) {
  side <- if (identical(x$alternative, "two.sided")) "two-sided" else "one-sided"
  cat(sprintf("Power of Fisher's Exact Test (%s)\n\n", side))
  cat("  n1 =", x$n1, " n2 =", x$n2,
      " (ratio =", sprintf("%.2f", x$n1 / x$n2), ")\n")
  cat("  p1 =", x$p1, " p2 =", x$p2, "\n")
  cat("  alpha =", x$alpha,
      " alpha* =", sprintf("%.5f", x$alpha_star), "\n")
  cat("  power =", sprintf("%.5f", x$power), "\n")
  cat("  error bound =", sprintf("%.2e", x$error), "\n")
  cat("  points evaluated =", x$count, "\n")
  invisible(x)
}


#' Sample Size for Fisher's Exact Test with Unequal Groups
#'
#' Find the minimum n2 (and corresponding n1 = ceiling(k * n2))
#' such that the power of Fisher's exact test reaches a target.
#'
#' @param p1 Numeric. True event probability in group 1.
#' @param p2 Numeric. True event probability in group 2.
#' @param alpha Numeric. One-sided type I error rate.
#' @param target_power Numeric. Desired power (default 0.80).
#' @param k Numeric. Allocation ratio n1/n2 (default 1).
#' @param eps Numeric. Trimming parameter for power calculation.
#' @param n_start Integer. Starting sample size for search.
#' @param n_max Integer. Maximum n2 to search.
#'
#' @return A list of class \code{"fxpower"} at the minimum
#'   sample size, plus component \code{n_total} = n1 + n2.
#'
#' @examples
#' ss_fxpower(p1 = 0.05, p2 = 0.10, alpha = 0.05, k = 1)
#' ss_fxpower(p1 = 0.05, p2 = 0.10, alpha = 0.05, k = 2)
#'
#' @export
ss_fxpower <- function(p1, p2, alpha = 0.05, target_power = 0.80,
                       k = 1, eps = 1e-6, n_start = 10L,
                       n_max = 100000L) {
  stopifnot(k > 0)

  n2 <- as.integer(n_start)

  repeat {
    n1 <- as.integer(ceiling(k * n2))
    res <- fxpower(n1 = n1, n2 = n2, p1 = p1, p2 = p2,
                   alpha = alpha, eps = eps)
    if (res$power >= target_power) break
    n2 <- n2 + 1L
    if (n2 > n_max) {
      warning("n_max reached without achieving target power")
      break
    }
  }

  res$n_total <- n1 + n2
  res
}


#' Power Table Across Allocation Ratios
#'
#' Compute power for a grid of sample sizes and allocation ratios,
#' producing a data frame suitable for tables in the paper.
#'
#' @param n2_vec Integer vector. Control group sizes.
#' @param k_vec Numeric vector. Allocation ratios n1/n2.
#' @param p1 Numeric. True event probability in group 1.
#' @param p2 Numeric. True event probability in group 2.
#' @param alpha Numeric. One-sided type I error rate.
#' @param eps Numeric. Trimming parameter.
#'
#' @return A data frame with columns n1, n2, k, ntotal, power,
#'   alpha_star, count.
#'
#' @examples
#' power_table(n2_vec = c(50, 100), k_vec = c(1, 1.5, 2),
#'             p1 = 0.05, p2 = 0.10)
#'
#' @export
power_table <- function(n2_vec, k_vec, p1, p2,
                        alpha = 0.05, eps = 1e-6) {
  grid <- expand.grid(n2 = n2_vec, k = k_vec)
  grid$n1 <- as.integer(ceiling(grid$k * grid$n2))
  grid$ntotal <- grid$n1 + grid$n2

  results <- mapply(
    function(n1, n2) {
      res <- fxpower(n1, n2, p1, p2, alpha, eps)
      c(power = res$power, alpha_star = res$alpha_star,
        count = res$count)
    },
    grid$n1, grid$n2,
    SIMPLIFY = TRUE
  )

  grid$power <- results["power", ]
  grid$alpha_star <- results["alpha_star", ]
  grid$count <- as.integer(results["count", ])

  grid[order(grid$k, grid$n2), c("n1", "n2", "k", "ntotal",
                                  "power", "alpha_star", "count")]
}


#' Sample Size Table Across Allocation Ratios
#'
#' Find the minimum sample size achieving target power for
#' each combination of parameter differences and allocation ratios.
#'
#' @param p2_vec Numeric vector. Control event probabilities.
#' @param delta_vec Numeric vector. Absolute differences p1 - p2
#'   (negative for lower treatment rates).
#' @param k_vec Numeric vector. Allocation ratios n1/n2.
#' @param alpha Numeric. One-sided type I error rate.
#' @param target_power Numeric. Desired power.
#' @param eps Numeric. Trimming parameter.
#'
#' @return A data frame with columns p2, delta, p1, k, n1, n2,
#'   ntotal, power.
#'
#' @examples
#' ss_table(p2_vec = 0.10, delta_vec = -0.05,
#'          k_vec = c(1, 2), target_power = 0.80)
#'
#' @export
ss_table <- function(p2_vec, delta_vec, k_vec,
                     alpha = 0.05, target_power = 0.80,
                     eps = 1e-6) {
  grid <- expand.grid(p2 = p2_vec, delta = delta_vec, k = k_vec)
  grid$p1 <- grid$p2 + grid$delta

  valid <- grid$p1 > 0 & grid$p1 < 1
  grid <- grid[valid, ]

  results <- mapply(
    function(p1, p2, k) {
      res <- ss_fxpower(p1, p2, alpha, target_power, k, eps)
      c(n1 = res$n1, n2 = res$n2, ntotal = res$n_total,
        power = res$power)
    },
    grid$p1, grid$p2, grid$k,
    SIMPLIFY = TRUE
  )

  grid$n1 <- as.integer(results["n1", ])
  grid$n2 <- as.integer(results["n2", ])
  grid$ntotal <- as.integer(results["ntotal", ])
  grid$power <- results["power", ]

  grid[order(grid$p2, grid$delta, grid$k),
       c("p2", "delta", "p1", "k", "n1", "n2", "ntotal", "power")]
}
