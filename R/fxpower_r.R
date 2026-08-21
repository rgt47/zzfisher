#' Pure R Reference Implementation of the Power Function
#'
#' Computes the exact power of Fisher's one-sided test by direct
#' enumeration. Slow but transparent -- used for validating the C
#' implementation.
#'
#' @param n1 Integer. Sample size in group 1.
#' @param n2 Integer. Sample size in group 2.
#' @param p1 Numeric. True event probability in group 1.
#' @param p2 Numeric. True event probability in group 2.
#' @param alpha Numeric. One-sided type I error rate.
#'
#' @return Numeric scalar: the exact power.
#'
#' @examples
#' fxpower_r(50, 50, p1 = 0.10, p2 = 0.30, alpha = 0.05)
#'
#' @export
fxpower_r <- function(n1, n2, p1, p2, alpha = 0.05) {
  ntot <- n1 + n2
  lfact <- c(0, cumsum(log(seq_len(ntot))))

  log_dhyper <- function(x, n1, n2, r) {
    y <- r - x
    lfact[n1 + 1] - lfact[x + 1] - lfact[n1 - x + 1] +
      lfact[n2 + 1] - lfact[y + 1] - lfact[n2 - y + 1] -
      lfact[ntot + 1] + lfact[r + 1] + lfact[ntot - r + 1]
  }

  log_dbinom_joint <- function(x, y) {
    lfact[n1 + 1] - lfact[x + 1] - lfact[n1 - x + 1] +
      lfact[n2 + 1] - lfact[y + 1] - lfact[n2 - y + 1] +
      x * log(p1) + (n1 - x) * log(1 - p1) +
      y * log(p2) + (n2 - y) * log(1 - p2)
  }

  powsum <- 0

  for (r in 0:ntot) {
    xmin_r <- max(0L, r - n2)
    xmax_r <- min(r, n1)
    if (xmax_r < xmin_r) next

    xs <- xmin_r:xmax_r
    hyp_probs <- exp(vapply(xs, log_dhyper,
                            double(1), n1 = n1, n2 = n2, r = r))
    cum_hyp <- cumsum(hyp_probs)
    xcv_idx <- max(which(cum_hyp <= alpha), 0L)
    if (xcv_idx == 0L) next
    xcv <- xs[xcv_idx]

    for (x in xmin_r:xcv) {
      y <- r - x
      powsum <- powsum + exp(log_dbinom_joint(x, y))
    }
  }

  powsum
}
