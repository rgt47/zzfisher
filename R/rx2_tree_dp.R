# tree_v5.R
# Version 5: Best pure-R tree implementation
#
# Base: tree_v4 (memoized find_max/find_min, closure-based scoping)
# Added from network v5/v6 (backported to tree paradigm):
#   (1) Precomputed log-factorial table — replaces lgamma(x+1) with
#       O(1) lookup, used in compute_prob
#   (2) Row ordering — sort rows by decreasing margin for tighter
#       early pruning
#   (3) Constrained suffix_max DP — cmax[k, c1] = max of
#       prod(choose(r[i], y[i])) for rows k..m subject to
#       sum(y[i]) = c1. Replaces the unconstrained product
#       prod(choose(r[i], floor(r[i]/2))) used in v3/v4.
#
# Matches tree_v5_opt.cpp algorithmically (minus C++-specific
# optimizations like fixed arrays and inline functions).

.rx2_tree_dp <- function(dat) {
  pval <- 1
  tol <- 3.45254e-7
  m <- nrow(dat)
  r_orig <- rowSums(dat)
  cc <- colSums(dat)
  n <- sum(r_orig)

  flipped <- cc[1] > cc[2]
  if (flipped) cc <- rev(cc)

  # (2) Row ordering: sort by decreasing margin
  ord <- order(r_orig, decreasing = TRUE)
  r <- r_orig[ord]

  y_obs <- if (flipped) dat[ord, 2] else dat[ord, 1]

  # (1) Precomputed log-factorial table
  lfact <- c(0, cumsum(log(seq_len(n))))

  # Suffix sums of row margins (1-indexed: suffix_r[k] = sum of r[k:m])
  suffix_r <- c(rev(cumsum(rev(r))), 0)

  log_const <- sum(lfact[r + 1]) + lfact[cc[1] + 1] + lfact[cc[2] + 1] - lfact[n + 1]
  p_obs <- exp(log_const - sum(lfact[y_obs + 1]) - sum(lfact[r - y_obs + 1]))

  # (3) Constrained suffix_max DP
  # cmax[k, c1] = max of prod(choose(r[i], y[i])) for rows k..m
  # subject to sum(y[i]) = c1.
  # Stored as matrix: cmax_mat[k, c1+1] (1-indexed k, 0-indexed c1)
  c0p1 <- cc[1] + 1
  cmax_mat <- matrix(0, nrow = m + 1, ncol = c0p1)
  cmax_mat[m + 1, 1] <- 1  # base: empty product = 1, c1 = 0

  for (k in m:1) {
    for (c1_idx in seq_len(c0p1)) {
      c1_val <- c1_idx - 1
      y_lo <- max(0, c1_val - suffix_r[k + 1])
      y_hi <- min(r[k], c1_val)
      if (y_lo > y_hi) next

      binom_val <- choose(r[k], y_lo)
      best <- 0

      for (y_k in y_lo:y_hi) {
        val <- binom_val * cmax_mat[k + 1, c1_val - y_k + 1]
        if (val > best) best <- val
        if (y_k < y_hi) {
          binom_val <- binom_val * (r[k] - y_k) / (y_k + 1)
        }
      }
      cmax_mat[k, c1_idx] <- best
    }
  }

  # Memoization caches
  memo_max <- vector("list", m)
  memo_min <- vector("list", m)
  for (i in seq_len(m)) {
    memo_max[[i]] <- vector("list", c0p1)
    memo_min[[i]] <- vector("list", c0p1)
  }

  y_max <- numeric(m)
  y_min <- numeric(m)
  p_max <- 0
  p_min <- 1

  compute_prob <- function(y) {
    exp(log_const - sum(lfact[y + 1]) - sum(lfact[r - y + 1]))
  }

  # --- find_min with memoization ---

  find_min <- function(c1, y, d) {
    k_cache <- d + 1
    cache_key <- c1 + 1

    if (!is.null(memo_min[[k_cache]][[cache_key]])) {
      cached <- memo_min[[k_cache]][[cache_key]]
      for (i in k_cache:m) y[i] <- cached$y[i]
      p_min <<- compute_prob(y)
      y_min <<- y
      return()
    }

    r_avail <- r
    if (d > 0) r_avail[1:d] <- NA
    n_rem <- sum(r[(d + 1):m])

    local_min_p <- Inf
    local_min_y <- y

    joe_min <- function(k, descend, r_avail, c1, n_rem, y) {
      if (k == m) {
        y[!is.na(r_avail)] <- c1
        prob <- compute_prob(y)
        if (prob < local_min_p) {
          local_min_p <<- prob
          local_min_y <<- y
        }
        if (prob < p_min) {
          y_min <<- y
          p_min <<- prob
        }
        return()
      }

      r_avail_save <- r_avail
      c1_save <- c1
      n_rem_save <- n_rem
      y_save <- y

      idx <- which.max(r_avail)
      y[idx] <- if (descend) min(r_avail[idx], c1) else r_avail[idx] - min(r_avail[idx], n_rem - c1)
      n_rem <- n_rem - r_avail[idx]
      r_avail[idx] <- NA
      c1 <- c1 - y[idx]

      joe_min(k + 1, TRUE, r_avail, c1, n_rem, y)

      if (descend) {
        joe_min(k, FALSE, r_avail_save, c1_save, n_rem_save, y_save)
      }
    }

    joe_min(d, TRUE, r_avail, c1, n_rem, y)

    memo_min[[k_cache]][[cache_key]] <<- list(p = local_min_p, y = local_min_y)
  }

  # --- find_max with memoization ---

  find_max <- function(k_start, c1, n_rem, y) {
    cache_key <- c1 + 1

    if (!is.null(memo_max[[k_start]][[cache_key]])) {
      cached <- memo_max[[k_start]][[cache_key]]
      for (i in k_start:m) y[i] <- cached$y[i]
      p_max <<- compute_prob(y)
      y_max <<- y
      return()
    }

    local_max_p <- 0
    local_max_y <- y

    stack <- list(list(k = k_start, c1 = c1, n_rem = n_rem, y = y))

    while (length(stack) > 0) {
      state <- stack[[length(stack)]]
      stack <- stack[-length(stack)]

      k <- state$k
      c1_rem <- state$c1
      n_rem_curr <- state$n_rem
      y_curr <- state$y

      if (k > m) {
        prob <- compute_prob(y_curr)
        if (prob > local_max_p) {
          local_max_p <- prob
          local_max_y <- y_curr
        }
        if (prob > p_max) {
          p_max <<- prob
          y_max <<- y_curr
        }
        next
      }

      d <- m - (k - 1)
      denom <- n_rem_curr + d
      if (denom == 0 || c1_rem == 0) {
        y_lo <- y_up <- 0
      } else {
        y_up <- floor((r[k] + 1) * (c1_rem + d - 1) / denom)
        y_lo <- ceiling(((r[k] + 1) * (c1_rem + 1) / denom) - 1)
      }
      if (c1_rem == 0) y_lo <- y_up <- 0

      n_after <- suffix_r[k + 1]
      y_lo <- max(y_lo, max(0, c1_rem - n_after))
      y_up <- min(y_up, min(r[k], c1_rem))

      if (y_lo > y_up) next

      for (y_k in y_up:y_lo) {
        y_new <- y_curr
        y_new[k] <- y_k
        stack[[length(stack) + 1]] <- list(
          k = k + 1,
          c1 = c1_rem - y_k,
          n_rem = n_rem_curr - r[k],
          y = y_new
        )
      }
    }

    memo_max[[k_start]][[cache_key]] <<- list(p = local_max_p, y = local_max_y)
  }

  # --- traverse with constrained suffix_max pruning ---

  traverse <- function(k, c1, y, prob_prefix, n_rem, descend, y_k, dir, mode) {
    if (k == m) {
      y[m] <- c1
      prob <- compute_prob(y)
      if (prob > p_obs * (1 + tol))
        pval <<- pval - prob
      return()
    }

    y_lo <- max(0, c1 - suffix_r[k + 1])
    y_hi <- min(r[k], c1)

    if (descend) {
      mode <- y_max[k]
      y_k <- mode
    }

    if (y_k < y_lo) return()
    if (y_k > y_hi) {
      traverse(k, c1, y, prob_prefix, n_rem, FALSE, mode - 1, -1, mode)
      return()
    }

    y[k] <- y_k
    hp <- dhyper(y_k, c1, n_rem - c1, r[k])
    new_prefix <- prob_prefix * hp
    new_c1 <- c1 - y_k
    new_n_rem <- suffix_r[k + 1]

    step_horizontal <- function() {
      traverse(k, c1, y, prob_prefix, n_rem, FALSE, y_k + dir, dir, mode)
    }

    if (m - k > 1) {
      # (3) Constrained suffix_max pruning
      smax <- cmax_mat[k + 1, new_c1 + 1]
      if (new_prefix * smax <= p_obs * (1 + tol)) {
        step_horizontal()
        return()
      }

      p_max <<- 0
      find_max(k + 1, new_c1, new_n_rem, y)
      p_min <<- 1
      find_min(new_c1, y, k)

      if (p_min > p_obs * (1 + tol)) {
        pval <<- pval - new_prefix
        step_horizontal()
        return()
      }

      if (p_max <= p_obs * (1 + tol)) {
        step_horizontal()
        return()
      }
    }

    traverse(k + 1, new_c1, y, new_prefix, new_n_rem, TRUE, 0, +1, 0)
    step_horizontal()
  }

  # ----- Main -----

  find_max(1, cc[1], n, y_obs)
  traverse(1, cc[1], y_obs, 1, n, TRUE, 0, +1, 0)

  pval
}
