make_htest <- function(p.value, data_name,
                       method_detail = "") {
  # Every kernel returns through here, which makes it the one place to
  # insist that what comes back is a p-value. Two of the network
  # kernels returned NaN on a table whose total count overflowed their
  # arithmetic, and NaN in a p-value slot is worse than an error: it
  # compares FALSE against any threshold, so a test silently stops
  # rejecting and nothing in the output says so.
  if (length(p.value) != 1L || !is.numeric(p.value) ||
        !is.finite(p.value)) {
    stop("the kernel returned ", format(p.value),
         " rather than a p-value; this table is outside the range it ",
         "can handle. Try another algorithm, or reduce the table's ",
         "total count.", call. = FALSE)
  }
  if (p.value < -1e-9 || p.value > 1 + 1e-9) {
    stop("the kernel returned ", format(p.value),
         ", which is not a probability.", call. = FALSE)
  }
  p.value <- min(max(p.value, 0), 1)

  method <- if (nzchar(method_detail))
    paste("Fisher's Exact Test for Count Data",
          method_detail, sep = " ")
  else
    "Fisher's Exact Test for Count Data"
  structure(
    list(
      p.value     = p.value,
      alternative = "two.sided",
      method      = method,
      data.name   = data_name
    ),
    class = "htest"
  )
}
