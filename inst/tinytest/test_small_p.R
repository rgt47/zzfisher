library(tinytest)

# Relative-error gate. An absolute-error check cannot see a p-value
# of 1e-20 reported as 0 or as 1e-14, so every comparison here is
# relative to an oracle that is itself accurate for tiny p: complete
# enumeration of the reference set, summed in log space.

.lse <- function(x) {
  if (!length(x)) return(-Inf)
  mx <- max(x)
  if (!is.finite(mx)) return(mx)
  mx + log(sum(exp(x - mx)))
}

# All tables with the margins of `dat`, as a vector of log P.
.enum_logp <- function(dat) {
  dat <- as.matrix(dat)
  r <- rowSums(dat)
  cs <- colSums(dat)
  n <- sum(dat)
  m <- nrow(dat)
  nc <- ncol(dat)
  lconst <- sum(lfactorial(r)) + sum(lfactorial(cs)) - lfactorial(n)
  out <- numeric(0)
  comps <- function(tot, cap) {
    if (length(cap) == 1L) {
      return(if (tot <= cap) matrix(tot, 1L) else matrix(0L, 0L, 1L))
    }
    res <- list()
    for (v in 0:min(tot, cap[1])) {
      sub <- comps(tot - v, cap[-1])
      if (nrow(sub)) res[[length(res) + 1L]] <- cbind(v, sub)
    }
    if (!length(res)) return(matrix(0L, 0L, length(cap)))
    do.call(rbind, res)
  }
  rec <- function(i, cap, acc) {
    if (i == m) {
      out <<- c(out, acc - sum(lfactorial(cap)))
      return(invisible())
    }
    cc <- comps(r[i], cap)
    for (j in seq_len(nrow(cc))) {
      rec(i + 1L, cap - cc[j, ], acc - sum(lfactorial(cc[j, ])))
    }
  }
  if (m == 1L) return(lconst - sum(lfactorial(dat)))
  rec(1L, cs, 0)
  lconst + out
}

.oracle_p <- function(dat, tol = 3.45254e-7) {
  lp <- .enum_logp(dat)
  lobs <- sum(lfactorial(rowSums(dat))) + sum(lfactorial(colSums(dat))) -
    lfactorial(sum(dat)) - sum(lfactorial(dat))
  min(1, exp(.lse(lp[lp <= lobs + log1p(tol)])))
}

.rel_ok <- function(got, ref, rtol = 1e-9) {
  is.finite(got) && abs(got - ref) <= rtol * ref
}

# The oracle agrees with fisher.test where FEXACT is accurate.
expect_true(.rel_ok(.oracle_p(matrix(c(1, 9, 11, 3), 2)),
                    fisher.test(matrix(c(1, 9, 11, 3), 2))$p.value),
            info = 'oracle matches fisher.test on a 2x2 table')

rx2_tabs <- list(
  tiny_2x2      = matrix(c(10, 0, 0, 10), 2),
  tiny_3x2      = matrix(c(20, 0, 0, 0, 20, 20), 3),
  tiny_3x2_b    = matrix(c(30, 0, 0, 0, 30, 30), 3),
  tiny_4x2      = matrix(c(40, 0, 0, 0, 0, 40, 40, 40), 4),
  tiny_5x2      = matrix(c(25, 2, 1, 0, 0, 1, 24, 25, 26, 25), 5),
  flipped_5x2   = matrix(c(1, 24, 25, 26, 25, 25, 2, 1, 0, 0), 5),
  moderate_4x2  = matrix(c(3, 5, 2, 7, 6, 2, 8, 1), 4),
  zero_rows_4x2 = matrix(c(0, 5, 9, 0, 7, 0, 1, 6), 4),
  p_one_2x2     = matrix(c(1, 1, 1, 1), 2),
  p_one_3x2     = matrix(c(2, 2, 2, 2, 2, 2), 3),
  two_rows      = matrix(c(12, 1, 0, 14), 2),
  # A row margin of 248 overflowed the sweep's leaf factor to
  # inf * 0 = NaN (Bedeian and Armenakis 1977, via Berry et al.).
  large_margin  = matrix(c(2, 2, 142, 0, 2, 106), 3),
  large_small_p = matrix(c(0, 6, 210, 6, 0, 170), 3)
)
rx2_fns <- list(
  tree_memo      = tree_memo,
  tree_dp        = tree_dp,
  tree_memo_c    = tree_memo_c,
  tree_s4_c      = tree_s4_c,
  tree_memo_cpp  = tree_memo_cpp,
  tree_dp_cpp    = tree_dp_cpp,
  net_vander_cpp = net_vander_cpp,
  net_dp_cpp     = net_dp_cpp,
  net_sweep_cpp  = net_sweep_cpp
)
for (tn in names(rx2_tabs)) {
  dat <- rx2_tabs[[tn]]
  ref <- .oracle_p(dat)
  for (fn in names(rx2_fns)) {
    got <- rx2_fns[[fn]](dat)$p.value
    expect_true(.rel_ok(got, ref),
                info = sprintf('%s on %s: got %.6g, oracle %.6g',
                               fn, tn, got, ref))
  }
}

rxc_tabs <- list(
  tiny_3x3     = matrix(c(8, 0, 0, 0, 8, 0, 0, 0, 8), 3),
  tiny_2x3     = matrix(c(15, 0, 0, 15, 0, 0), 2),
  moderate_3x3 = matrix(c(3, 1, 2, 0, 4, 1, 2, 2, 3), 3),
  p_one_2x3    = matrix(c(1, 1, 1, 1, 1, 1), 2),
  zero_col_3x3 = matrix(c(4, 0, 1, 0, 0, 3, 2, 5, 0), 3)
)
# Exported dispatchers only, so the file runs unchanged in any
# package that ships these kernels: tree_memo_c routes r x c to the
# r x c tree kernel and net_sweep_cpp to the r x c sweep.
rxc_fns <- list(
  rxc_tree_memo = function(d) tree_memo_c(d),
  rxc_net_sweep = function(d) net_sweep_cpp(d)
)
for (tn in names(rxc_tabs)) {
  dat <- rxc_tabs[[tn]]
  ref <- .oracle_p(dat)
  for (fn in names(rxc_fns)) {
    got <- rxc_fns[[fn]](dat)
    if (is.list(got)) got <- got$p.value
    expect_true(.rel_ok(got, ref),
                info = sprintf('%s on %s: got %.6g, oracle %.6g',
                               fn, tn, got, ref))
  }
}

# r x c x k, no-three-way-interaction kernel. For 2 x 2 x 2 the
# fiber with all two-way margins fixed is one-dimensional, so it
# is enumerated directly; weights are summed relative to the
# smallest T, as in the kernel.
.oracle_rxck_222 <- function(a, tol = 3.45254e-7) {
  ab <- apply(a, c(1, 2), sum)
  ac <- apply(a, c(1, 3), sum)
  bc <- apply(a, c(2, 3), sum)
  tt <- numeric(0)
  for (x in 0:sum(a)) {
    y <- array(0, c(2, 2, 2))
    y[1, 1, 1] <- x
    y[1, 1, 2] <- ab[1, 1] - x
    y[1, 2, 1] <- ac[1, 1] - x
    y[2, 1, 1] <- bc[1, 1] - x
    y[1, 2, 2] <- ab[1, 2] - y[1, 2, 1]
    y[2, 1, 2] <- ab[2, 1] - y[2, 1, 1]
    y[2, 2, 1] <- ac[2, 1] - y[2, 1, 1]
    y[2, 2, 2] <- ab[2, 2] - y[2, 2, 1]
    if (any(y < 0)) next
    if (!all(apply(y, c(1, 3), sum) == ac) ||
        !all(apply(y, c(2, 3), sum) == bc)) next
    tt <- c(tt, sum(lfactorial(y)))
  }
  t_obs <- sum(lfactorial(a))
  w <- exp(min(tt) - tt)
  sum(w[tt - t_obs >= -tol]) / sum(w)
}
rxck_tabs <- list(
  large_n = array(c(100, 95, 90, 105, 98, 97, 102, 99), c(2, 2, 2)),
  small_p = array(c(12, 0, 0, 12, 0, 11, 11, 0), c(2, 2, 2)),
  moderate = array(c(3, 5, 2, 4, 6, 1, 2, 7), c(2, 2, 2))
)
for (tn in names(rxck_tabs)) {
  a <- rxck_tabs[[tn]]
  storage.mode(a) <- 'integer'
  ref <- .oracle_rxck_222(a)
  got <- tryCatch(tree_memo_c(a)$p.value, error = function(e) NA_real_)
  expect_true(.rel_ok(got, ref),
              info = sprintf('rxck_tree_memo on %s: got %.6g, oracle %.6g',
                             tn, got, ref))
}

# Stratified conditional-independence kernel: the reference set is
# the product of the per-stratum fibers, enumerated directly.
.oracle_ci <- function(a, tol = 3.45254e-7) {
  lp <- 0
  lobs <- 0
  for (l in seq_len(dim(a)[3])) {
    s <- a[, , l]
    lp <- outer(lp, .enum_logp(s), `+`)
    lobs <- lobs + sum(lfactorial(rowSums(s))) +
      sum(lfactorial(colSums(s))) - lfactorial(sum(s)) -
      sum(lfactorial(s))
  }
  lp <- as.vector(lp)
  min(1, exp(.lse(lp[lp <= lobs + log1p(tol)])))
}
ci_tabs <- list(
  tiny_2strata = array(c(9, 0, 0, 9, 8, 0, 0, 8), c(2, 2, 2)),
  tiny_3x3x2 = array(c(5, 0, 0, 0, 5, 0, 0, 0, 5,
                       4, 0, 0, 0, 4, 0, 0, 0, 4), c(3, 3, 2)),
  moderate = array(c(3, 5, 2, 4, 6, 1, 2, 7), c(2, 2, 2))
)
for (tn in names(ci_tabs)) {
  a <- ci_tabs[[tn]]
  storage.mode(a) <- 'integer'
  ref <- .oracle_ci(a)
  got <- net_ci_cpp(a)$p.value
  expect_true(.rel_ok(got, ref),
              info = sprintf('net_ci on %s: got %.6g, oracle %.6g',
                             tn, got, ref))
}
