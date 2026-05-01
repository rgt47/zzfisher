// tree_memo_s4_wrapper.cpp
// Rcpp wrapper for tree_memo_s4_c_impl (S4 mode-path dominance)

#include <Rcpp.h>

extern "C" double tree_memo_s4_c_impl(int *dat, int m);

// [[Rcpp::export(name = ".tree_memo_s4_c")]]
double tree_memo_s4_c(Rcpp::IntegerMatrix dat) {
    if (dat.ncol() != 2)
        Rcpp::stop("Input must be an r x 2 matrix");
    return tree_memo_s4_c_impl(dat.begin(), dat.nrow());
}
