// rx2_tree_memo_c_wrapper.cpp
// Rcpp wrapper for rx2_tree_memo_c_impl (pure C, r x 2 case)

#include <Rcpp.h>

extern "C" double rx2_tree_memo_c_impl(int *dat, int m);

// [[Rcpp::export(name = ".rx2_tree_memo_c")]]
double rx2_tree_memo_c(Rcpp::IntegerMatrix dat) {
    if (dat.ncol() != 2)
        Rcpp::stop("Input must be an r x 2 matrix");
    return rx2_tree_memo_c_impl(dat.begin(), dat.nrow());
}
