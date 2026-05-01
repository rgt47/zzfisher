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

test_that(paste(nm, "returns htest class"), {
    result <- fn(dat_2x2)
    expect_inherits(result, "htest")


test_that(paste(nm, "p-value matches fisher.test (2x2)"), {
    result <- fn(dat_2x2)
    expect_true(abs(result$p.value - ref_2x2) < 1e-12)


test_that(paste(nm, "p-value matches fisher.test (3x2)"), {
    result <- fn(dat_3x2)
    expect_true(abs(result$p.value - ref_3x2) < 1e-12)


test_that(paste(nm, "method contains function name"), {
    result <- fn(dat_2x2)
    expect_true(grepl(nm, result$method))
    expect_true(grepl("Fisher's Exact Test", result$method))


test_that(paste(nm, "alternative is two.sided"), {
    result <- fn(dat_2x2)
    expect_identical(result$alternative, "two.sided")


test_that(paste(nm, "data.name captures expression"), {
    result <- fn(dat_2x2)
    expect_identical(result$data.name, "dat_2x2")

}

# tree_memo_profile returns htest with profile
skip_if(
    tryCatch(
      { .tree_memo_profile(dat_2x2); FALSE },
      error = function(e) grepl("PROFILE_V4", e$message)
    ),
    "PROFILE_V4 not enabled"
)
result <- tree_memo_profile(dat_2x2)
expect_inherits(result, "htest")
expect_true(grepl("tree_memo_profile", result$method))
expect_true("profile" %in% names(result))

