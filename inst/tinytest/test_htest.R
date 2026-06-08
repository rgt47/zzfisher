library(tinytest)

dat_2x2 <- matrix(c(1, 9, 11, 3), nrow = 2)
dat_3x2 <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 3)
ref_2x2 <- fisher.test(dat_2x2)$p.value
ref_3x2 <- fisher.test(dat_3x2)$p.value

all_fns <- list(
  tree_memo      = tree_memo,
  tree_dp        = tree_dp,
  tree_memo_c    = tree_memo_c,
  tree_memo_s4_c = tree_memo_s4_c,
  tree_memo_cpp  = tree_memo_cpp,
  tree_dp_cpp    = tree_dp_cpp,
  net_vander_cpp = net_vander_cpp,
  net_dp_cpp     = net_dp_cpp
)

for (nm in names(all_fns)) {
  fn <- all_fns[[nm]]
  result <- fn(dat_2x2)

  expect_inherits(result, "htest",
    info = paste(nm, "returns htest class"))

  expect_true(abs(result$p.value - ref_2x2) < 1e-12,
    info = paste(nm, "p-value matches fisher.test (2x2)"))

  expect_true(abs(fn(dat_3x2)$p.value - ref_3x2) < 1e-12,
    info = paste(nm, "p-value matches fisher.test (3x2)"))

  expect_true(grepl(nm, result$method),
    info = paste(nm, "method contains function name"))

  expect_true(grepl("Fisher's Exact Test", result$method),
    info = paste(nm, "method contains Fisher's Exact Test"))

  expect_identical(result$alternative, "two.sided",
    info = paste(nm, "alternative is two.sided"))

  expect_identical(result$data.name, "dat_2x2",
    info = paste(nm, "data.name captures expression"))
}

# tree_memo_profile returns htest with profile (only when PROFILE_V4 compiled in)
run_profile <- tryCatch(
  { .tree_memo_profile(dat_2x2); TRUE },
  error = function(e) FALSE
)

if (run_profile) {
  result <- tree_memo_profile(dat_2x2)
  expect_inherits(result, "htest")
  expect_true(grepl("tree_memo_profile", result$method))
  expect_true("profile" %in% names(result))
}
