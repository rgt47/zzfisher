# dispatch.R
# User-facing entry points. Each function captures the data name for
# the htest object, inspects dim(dat) to route to the rx2, rxc, or
# rxck kernel, and wraps the returned p-value in an htest S3 object.
#
# Internal kernels live in:
#   R/rx2_*.R          pure-R r x 2 implementations (one per algorithm)
#   src/rx2_*.{c,cpp}  compiled r x 2 implementations (one per algorithm)
#   src/rxc_tree_memo.cpp   first-gen r x c kernel (shared by all algos)
#   src/rxck_tree_memo.cpp  first-gen r x c x k kernel (shared by all)
#
# First-generation policy: every user-facing dispatcher (tree_memo,
# tree_dp, net_*, etc.) routes rxc and rxck shapes to the same pair
# of compiled kernels. The algorithm-specific rxc/rxck kernels for
# each paradigm (tree_dp, net_vander, net_dp) are second-generation
# work and will replace these fallback calls when written.
#
# tree_memo_profile is the exception: profiling is an rx2-specific
# feature (instrumentation of the memoized tree traversal) and stays
# stubbed for non-rx2 shapes.
#
# Compiled internals are auto-generated into R/RcppExports.R by
# Rcpp::compileAttributes().

.infer_shape <- function(dat) {
  d <- dim(dat)
  if (is.null(d) || !(length(d) %in% c(2L, 3L)))
    stop('dat must be a 2D matrix or 3D array', call. = FALSE)
  if (length(d) == 3L) return('rxck')
  if (d[2L] == 2L) 'rx2' else 'rxc'
}

.not_yet <- function(algo, shape) {
  stop(sprintf("'%s' is not yet implemented for %s tables",
               algo, shape),
       call. = FALSE)
}

# ----- Pure R -----

#' Memoized tree traversal (pure R)
#'
#' @param dat Integer matrix (r x 2 or r x c contingency table) or
#'   3D integer array (r x c x k).
#' @return An htest object matching \code{fisher.test()} output.
#' @export
tree_memo <- function(dat) {
  dname <- deparse(substitute(dat))
  pval <- switch(.infer_shape(dat),
    rx2  = .rx2_tree_memo(dat),
    rxc  = .rxc_tree_memo_cpp(dat),
    rxck = .rxck_tree_memo_cpp(dat)
  )
  make_htest(pval, dname, '[tree_memo]')
}

#' Tree traversal with constrained suffix DP (pure R)
#'
#' @param dat Integer matrix (r x 2 or r x c contingency table) or
#'   3D integer array (r x c x k).
#' @return An htest object matching \code{fisher.test()} output.
#' @export
tree_dp <- function(dat) {
  dname <- deparse(substitute(dat))
  pval <- switch(.infer_shape(dat),
    rx2  = .rx2_tree_dp(dat),
    rxc  = .rxc_tree_memo_cpp(dat),
    rxck = .rxck_tree_memo_cpp(dat)
  )
  make_htest(pval, dname, '[tree_dp]')
}

# ----- Compiled C -----

#' Memoized tree traversal (C implementation)
#'
#' @param dat Integer matrix (r x 2 or r x c contingency table) or
#'   3D integer array (r x c x k).
#' @return An htest object matching \code{fisher.test()} output.
#' @export
tree_memo_c <- function(dat) {
  dname <- deparse(substitute(dat))
  pval <- switch(.infer_shape(dat),
    rx2  = .rx2_tree_memo_c(dat),
    rxc  = .rxc_tree_memo_cpp(dat),
    rxck = .rxck_tree_memo_cpp(dat)
  )
  make_htest(pval, dname, '[tree_memo_c]')
}

#' Tree traversal with S4 mode-path dominance (C, no memoization)
#'
#' @param dat Integer matrix (r x 2 or r x c contingency table) or
#'   3D integer array (r x c x k).
#' @return An htest object matching \code{fisher.test()} output.
#' @export
tree_s4_c <- function(dat) {
  dname <- deparse(substitute(dat))
  pval <- switch(.infer_shape(dat),
    rx2  = .rx2_tree_s4_c(dat),
    rxc  = .rxc_tree_memo_cpp(dat),
    rxck = .rxck_tree_memo_cpp(dat)
  )
  make_htest(pval, dname, '[tree_s4_c]')
}

# ----- Compiled C++ -----

#' Memoized tree traversal (C++ implementation)
#'
#' @param dat Integer matrix (r x 2 or r x c contingency table) or
#'   3D integer array (r x c x k).
#' @return An htest object matching \code{fisher.test()} output.
#' @export
tree_memo_cpp <- function(dat) {
  dname <- deparse(substitute(dat))
  pval <- switch(.infer_shape(dat),
    rx2  = .rx2_tree_memo_cpp(dat),
    rxc  = .rxc_tree_memo_cpp(dat),
    rxck = .rxck_tree_memo_cpp(dat)
  )
  make_htest(pval, dname, '[tree_memo_cpp]')
}

#' Tree traversal with constrained suffix DP (C++)
#'
#' @param dat Integer matrix (r x 2 or r x c contingency table) or
#'   3D integer array (r x c x k).
#' @return An htest object matching \code{fisher.test()} output.
#' @export
tree_dp_cpp <- function(dat) {
  dname <- deparse(substitute(dat))
  pval <- switch(.infer_shape(dat),
    rx2  = .rx2_tree_dp_cpp(dat),
    rxc  = .rxc_tree_memo_cpp(dat),
    rxck = .rxck_tree_memo_cpp(dat)
  )
  make_htest(pval, dname, '[tree_dp_cpp]')
}

#' Network enumeration with Vandermonde bulk-sum pruning (C++)
#'
#' @param dat Integer matrix (r x 2 or r x c contingency table) or
#'   3D integer array (r x c x k).
#' @return An htest object matching \code{fisher.test()} output.
#' @export
net_vander_cpp <- function(dat) {
  dname <- deparse(substitute(dat))
  pval <- switch(.infer_shape(dat),
    rx2  = .rx2_net_vander_cpp(dat),
    rxc  = .rxc_tree_memo_cpp(dat),
    rxck = .rxck_tree_memo_cpp(dat)
  )
  make_htest(pval, dname, '[net_vander_cpp]')
}

#' Network enumeration with suffix DP and recurrence (C++)
#'
#' @param dat Integer matrix (r x 2 or r x c contingency table) or
#'   3D integer array (r x c x k).
#' @return An htest object matching \code{fisher.test()} output.
#' @export
net_dp_cpp <- function(dat) {
  dname <- deparse(substitute(dat))
  pval <- switch(.infer_shape(dat),
    rx2  = .rx2_net_dp_cpp(dat),
    rxc  = .rxc_tree_memo_cpp(dat),
    rxck = .rxck_tree_memo_cpp(dat)
  )
  make_htest(pval, dname, '[net_dp_cpp]')
}

#' Network state-sweep with prefix-class aggregation (C++)
#'
#' Stage-by-stage sweep over (stage, budget) states: paths are
#' aggregated into prefix log-mass classes per state, classified
#' against exact suffix extrema, and whole classes are subtracted
#' or dropped in O(1); at the last two rows each class contributes
#' one prefix-sum difference over shared leaf masses. The
#' Mehta-Patel network algorithm with the Clarkson-Fan-Joe past
#' device, specialized to two columns.
#'
#' @param dat Integer matrix (r x 2 or r x c contingency table) or
#'   3D integer array (r x c x k).
#' @return An htest object matching \code{fisher.test()} output.
#' @export
net_sweep_cpp <- function(dat) {
  dname <- deparse(substitute(dat))
  pval <- switch(.infer_shape(dat),
    rx2  = .rx2_net_sweep_cpp(dat),
    rxc  = .rxc_tree_memo_cpp(dat),
    rxck = .rxck_tree_memo_cpp(dat)
  )
  make_htest(pval, dname, '[net_sweep_cpp]')
}

# tree_memo_profile returns a list with $pvalue and profiling fields.
# The htest gets the profile attached on $profile.

#' Memoized tree traversal with profiling (C++)
#'
#' @param dat Integer matrix (r x 2 contingency table). Profiling is
#'   rx2-specific; r x c and r x c x k shapes are not yet supported.
#' @return An htest object with a \code{$profile} list attached.
#' @export
tree_memo_profile <- function(dat) {
  dname <- deparse(substitute(dat))
  raw <- switch(.infer_shape(dat),
    rx2  = .rx2_tree_memo_profile(dat),
    rxc  = .not_yet('tree_memo_profile', 'r x c'),
    rxck = .not_yet('tree_memo_profile', 'r x c x k')
  )
  result <- make_htest(raw$pvalue, dname, '[tree_memo_profile]')
  result$profile <- raw[setdiff(names(raw), 'pvalue')]
  result
}
