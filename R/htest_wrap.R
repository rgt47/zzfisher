make_htest <- function(p.value, data_name,
                       method_detail = "") {
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
