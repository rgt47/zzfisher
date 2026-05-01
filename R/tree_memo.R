#' Memoized tree traversal (pure R)
#'
#' @param dat Integer matrix with 2 columns (r x 2 contingency table).
#' @return An htest object matching \code{fisher.test()} output.
#' @export
tree_memo <- function(dat) {
  pval <- 1
  tol <- 3.45254e-7
  m <- nrow(dat)
  r_orig <- rowSums(dat)
  cc <- colSums(dat)
  n <- sum(r_orig)

  flipped <- cc[1] > cc[2]
  if (flipped) cc <- rev(cc)

  ord <- order(r_orig, decreasing = FALSE)
  r <- r_orig[ord]
  y_obs <- if (flipped) dat[ord, 2] else dat[ord, 1]

  # suffix_r[k] = sum of r[k:m]; suffix_r[m+1] = 0
  suffix_r <- c(rev(cumsum(rev(r))), 0)

  lfact <- lgamma(seq_len(n + 2))
  log_const <- sum(lfact[r + 1]) +
    sum(lfact[cc + 1]) - lfact[n + 1]
  p_obs <- exp(log_const - sum(lfact[y_obs + 1]) -
    sum(lfact[r - y_obs + 1]))

  y_max <- numeric(m)
  p_max <- 0
  p_min <- 1

  memo_max <- new.env(hash = TRUE, parent = emptyenv())
  memo_min <- new.env(hash = TRUE, parent = emptyenv())

  compute_prob <- function(y) {
    exp(log_const -
      sum(lfact[y + 1]) - sum(lfact[r - y + 1]))
  }

  dhyper_lfact <- function(x, m_par, n_par, k_par) {
    if (x < 0 || x > k_par || x > m_par ||
        k_par - x > n_par) return(0)
    lp <- lfact[m_par + 1] - lfact[x + 1] -
      lfact[m_par - x + 1] +
      lfact[n_par + 1] - lfact[k_par - x + 1] -
      lfact[n_par - k_par + x + 1] -
      lfact[m_par + n_par + 1] + lfact[k_par + 1] +
      lfact[m_par + n_par - k_par + 1]
    exp(lp)
  }

  find_min <- function(c1, y, d) {
    key <- paste(d, c1)
    if (exists(key, envir = memo_min)) {
      cached_suffix <- get(key, envir = memo_min)
      y[d:m] <- cached_suffix
      p_min <<- compute_prob(y)
      return()
    }

    r_avail <- r
    if (d > 1) r_avail[1:(d - 1)] <- NA
    n_rem <- suffix_r[d]

    local_min_p <- Inf
    local_min_y <- y

    joe_min <- function(k, descend, r_avail,
                        c1, n_rem, y) {
      if (k > m) {
        y[!is.na(r_avail)] <- c1
        prob <- compute_prob(y)
        if (prob < local_min_p) {
          local_min_p <<- prob
          local_min_y <<- y
        }
        return()
      }
      r_save <- r_avail
      c1_save <- c1
      n_save <- n_rem
      y_save <- y
      idx <- which.max(r_avail)
      y[idx] <- if (descend) min(r_avail[idx], c1)
        else r_avail[idx] - min(r_avail[idx], n_rem - c1)
      n_rem <- n_rem - r_avail[idx]
      r_avail[idx] <- NA
      c1 <- c1 - y[idx]
      joe_min(k + 1, TRUE, r_avail, c1, n_rem, y)
      if (descend)
        joe_min(k, FALSE, r_save, c1_save, n_save, y_save)
    }

    joe_min(d, TRUE, r_avail, c1, n_rem, y)
    p_min <<- local_min_p
    assign(key, local_min_y[d:m], envir = memo_min)
  }

  find_max <- function(k_start, c1, y) {
    key <- paste(k_start, c1)
    if (exists(key, envir = memo_max)) {
      cached_suffix <- get(key, envir = memo_max)
      y[k_start:m] <- cached_suffix
      p_max <<- compute_prob(y)
      y_max <<- y
      return()
    }

    nr <- suffix_r[k_start]
    local_max_p <- 0
    local_max_y <- y

    # Incremental log-prob: compute prefix once
    log_p_prefix <- log_const
    if (k_start > 1) {
      for (i in 1:(k_start - 1))
        log_p_prefix <- log_p_prefix -
          lfact[y[i] + 1] - lfact[r[i] - y[i] + 1]
    }

    stack <- list(list(
      k = k_start, c1 = c1,
      n_rem = nr, y = y, log_p = log_p_prefix))

    while (length(stack) > 0) {
      cur <- stack[[length(stack)]]
      stack[[length(stack)]] <- NULL
      k <- cur$k
      c1r <- cur$c1
      nr <- cur$n_rem
      yc <- cur$y

      if (k > m) {
        prob <- exp(cur$log_p)
        if (prob > local_max_p) {
          local_max_p <- prob
          local_max_y <- yc
        }
        next
      }

      d <- m - k + 1
      denom <- nr + d
      if (denom == 0 || c1r == 0) {
        y_lo <- 0L
        y_up <- 0L
      } else {
        y_up <- floor((r[k] + 1) * (c1r + d - 1) / denom)
        y_lo <- ceiling((r[k] + 1) * (c1r + 1) /
          denom) - 1
      }
      n_after <- suffix_r[k + 1]
      y_lo <- max(y_lo, max(0, c1r - n_after))
      y_up <- min(y_up, min(r[k], c1r))
      if (y_lo > y_up) next

      for (yk in y_lo:y_up) {
        yn <- yc
        yn[k] <- yk
        stack[[length(stack) + 1]] <- list(
          k = k + 1, c1 = c1r - yk,
          n_rem = nr - r[k], y = yn,
          log_p = cur$log_p -
            lfact[yk + 1] - lfact[r[k] - yk + 1])
      }
    }

    p_max <<- local_max_p
    y_max <<- local_max_y
    assign(key, local_max_y[k_start:m], envir = memo_max)
  }

  # Inlined handler for bottom two rows (k_penult = m-1 in R)
  penult <- function(c1, n_rem, prob_prefix) {
    k <- m - 1L
    y_lo <- max(0, c1 - r[m])
    y_hi <- min(r[k], c1)
    if (y_lo > y_hi) return()

    c2 <- n_rem - c1
    thresh <- p_obs * (1 + tol)

    mode_k <- y_max[k]
    if (mode_k < y_lo) mode_k <- y_lo
    if (mode_k > y_hi) mode_k <- y_hi

    hp_mode <- dhyper_lfact(mode_k, c1, c2, r[k])
    prob <- prob_prefix * hp_mode
    if (prob > thresh) pval <<- pval - prob

    # Forward from mode
    prev_yk <- mode_k
    hp <- hp_mode
    if (mode_k < y_hi) {
      for (y_k in (mode_k + 1):y_hi) {
        hp <- hp * (c1 - prev_yk) *
          (r[k] - prev_yk) /
          ((prev_yk + 1) *
            (c2 - r[k] + prev_yk + 1))
        prev_yk <- y_k
        prob <- prob_prefix * hp
        if (prob > thresh) pval <<- pval - prob
      }
    }

    # Backward from mode
    prev_yk <- mode_k
    hp <- hp_mode
    if (mode_k > y_lo) {
      for (y_k in (mode_k - 1):y_lo) {
        hp <- hp * prev_yk *
          (c2 - r[k] + prev_yk) /
          ((c1 - prev_yk + 1) *
            (r[k] - prev_yk + 1))
        prev_yk <- y_k
        prob <- prob_prefix * hp
        if (prob > thresh) pval <<- pval - prob
      }
    }
  }

  traverse <- function(k, c1, y, prob_prefix, n_rem,
                       descend, y_k, dir, mode) {
    # Terminal: last row is forced
    if (k >= m) {
      if (k == m) {
        y[m] <- c1
        hp_last <- dhyper_lfact(
          c1, c1, n_rem - c1, r[m])
        prob <- prob_prefix * hp_last
        if (prob > p_obs * (1 + tol))
          pval <<- pval - prob
      }
      return()
    }

    # dhyper recurrence state
    hp <- 0
    prev_yk <- -1L
    need_full <- TRUE

    repeat {
      if (k >= m) {
        if (k == m) {
          y[m] <- c1
          hp_last <- dhyper_lfact(
            c1, c1, n_rem - c1, r[m])
          prob <- prob_prefix * hp_last
          if (prob > p_obs * (1 + tol))
            pval <<- pval - prob
        }
        return()
      }

      y_lo <- max(0, c1 - suffix_r[k + 1])
      y_hi <- min(r[k], c1)

      if (descend) {
        mode <- y_max[k]
        y_k <- mode
        need_full <- TRUE
      }
      if (y_k < y_lo) return()
      if (y_k > y_hi) {
        y_k <- mode - 1L
        dir <- -1L
        descend <- FALSE
        need_full <- TRUE
        next
      }

      y[k] <- y_k
      c2 <- n_rem - c1

      if (need_full) {
        hp <- dhyper_lfact(y_k, c1, c2, r[k])
        need_full <- FALSE
      } else if (dir == 1) {
        hp <- hp * (c1 - prev_yk) *
          (r[k] - prev_yk) /
          ((prev_yk + 1) *
            (c2 - r[k] + prev_yk + 1))
      } else {
        hp <- hp * prev_yk *
          (c2 - r[k] + prev_yk) /
          ((c1 - prev_yk + 1) *
            (r[k] - prev_yk + 1))
      }
      prev_yk <- y_k

      new_prefix <- prob_prefix * hp
      new_c1 <- c1 - y_k
      new_n_rem <- suffix_r[k + 1]

      # S3-before-S2 pruning (skip S2 if S3 prunes)
      if (m - k > 1) {
        p_max <<- 0
        find_max(k + 1, new_c1, y)
        if (p_max <= p_obs * (1 + tol)) {
          y_k <- y_k + dir
          descend <- FALSE
          next
        }
        p_min <<- 1
        find_min(new_c1, y, k + 1)
        if (p_min > p_obs * (1 + tol)) {
          pval <<- pval - new_prefix
          y_k <- y_k + dir
          descend <- FALSE
          next
        }
      }

      # Dispatch: penult for bottom two rows, recurse otherwise
      if (k == m - 2) {
        penult(new_c1, new_n_rem, new_prefix)
      } else {
        traverse(k + 1, new_c1, y, new_prefix,
          new_n_rem, TRUE, 0L, 1L, 0L)
      }
      y_k <- y_k + dir
      descend <- FALSE
    }
  }

  dname <- deparse(substitute(dat))
  find_max(1, cc[1], y_obs)
  traverse(1, cc[1], y_obs, 1, n, TRUE, 0L, 1L, 0L)
  make_htest(pval, dname, "[tree_memo]")
}
