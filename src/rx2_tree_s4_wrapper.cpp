// rx2_tree_s4_wrapper.cpp
// Rcpp wrapper for rx2_tree_s4_c_impl (S4 mode-path dominance, r x 2)

#include <Rcpp.h>

extern "C" double rx2_tree_s4_c_impl(int *dat, int m);

// [[Rcpp::export(name = ".rx2_tree_s4_c")]]
double rx2_tree_s4_c(Rcpp::IntegerMatrix dat) {
    if (dat.ncol() != 2)
        Rcpp::stop("Input must be an r x 2 matrix");
    return rx2_tree_s4_c_impl(dat.begin(), dat.nrow());
}
