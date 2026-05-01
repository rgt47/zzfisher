#' Memoized tree traversal (C implementation)
#'
#' @param dat Integer matrix with 2 columns (r x 2 contingency table).
#' @return An htest object matching \code{fisher.test()} output.
#' @export
tree_memo_c <- function(dat) {
  dname <- deparse(substitute(dat))
  make_htest(.tree_memo_c(dat), dname, "[tree_memo_c]")
}

#' Memoized tree traversal with S4 mode-path dominance (C)
#'
#' @param dat Integer matrix with 2 columns (r x 2 contingency table).
#' @return An htest object matching \code{fisher.test()} output.
#' @export
tree_memo_s4_c <- function(dat) {
  dname <- deparse(substitute(dat))
  make_htest(.tree_memo_s4_c(dat), dname, "[tree_memo_s4_c]")
}

#' Memoized tree traversal (C++ implementation)
#'
#' @param dat Integer matrix with 2 columns (r x 2 contingency table).
#' @return An htest object matching \code{fisher.test()} output.
#' @export
tree_memo_cpp <- function(dat) {
  dname <- deparse(substitute(dat))
  make_htest(.tree_memo_cpp(dat), dname, "[tree_memo_cpp]")
}

#' Tree traversal with constrained suffix DP (C++)
#'
#' @param dat Integer matrix with 2 columns (r x 2 contingency table).
#' @return An htest object matching \code{fisher.test()} output.
#' @export
tree_dp_cpp <- function(dat) {
  dname <- deparse(substitute(dat))
  make_htest(.tree_dp_cpp(dat), dname, "[tree_dp_cpp]")
}

#' Network enumeration with Vandermonde bulk-sum pruning (C++)
#'
#' @param dat Integer matrix with 2 columns (r x 2 contingency table).
#' @return An htest object matching \code{fisher.test()} output.
#' @export
net_vander_cpp <- function(dat) {
  dname <- deparse(substitute(dat))
  make_htest(.net_vander_cpp(dat), dname, "[net_vander_cpp]")
}

#' Network enumeration with suffix DP and recurrence (C++)
#'
#' @param dat Integer matrix with 2 columns (r x 2 contingency table).
#' @return An htest object matching \code{fisher.test()} output.
#' @export
net_dp_cpp <- function(dat) {
  dname <- deparse(substitute(dat))
  make_htest(.net_dp_cpp(dat), dname, "[net_dp_cpp]")
}

#' Memoized tree traversal with profiling (C++)
#'
#' @param dat Integer matrix with 2 columns (r x 2 contingency table).
#' @return An htest object with a \code{$profile} list of counters.
#' @export
tree_memo_profile <- function(dat) {
  dname <- deparse(substitute(dat))
  raw <- .tree_memo_profile(dat)
  result <- make_htest(raw$pvalue, dname, "[tree_memo_profile]")
  result$profile <- raw[setdiff(names(raw), "pvalue")]
  result
}
