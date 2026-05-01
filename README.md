# zzfisher

Multiple algorithmic implementations of Fisher's exact test for
r x 2 contingency tables. Includes tree-based (mode-centric
traversal with pruning) and network-based (state-space enumeration
with Vandermonde bulk-sum) approaches in pure R, C, and C++.

## Installation

### From source (development version)

Requires a C/C++ compiler toolchain (Rtools on Windows, Xcode
Command Line Tools on macOS, `r-base-dev` on Debian/Ubuntu).

```r
# install.packages("pak")
pak::pak("path/to/zzfisher")
```

Or with devtools:

```r
# install.packages("devtools")
devtools::install("path/to/zzfisher")
```

Or from the command line:

```bash
R CMD build path/to/zzfisher
R CMD INSTALL zzfisher_0.1.0.tar.gz
```

### From GitHub

Once the repository is hosted on GitHub:

```r
pak::pak("username/zzfisher")
```

## Usage

All functions accept an integer matrix with two columns
(an r x 2 contingency table) and return a p-value.

```r
library(zzfisher)

dat <- matrix(c(1, 4, 3, 5,
                2, 3, 1, 4), ncol = 2)
dat
#>      [,1] [,2]
#> [1,]    1    2
#> [2,]    4    3
#> [3,]    3    1
#> [4,]    5    4

tree_memo(dat)
tree_dp(dat)
tree_dp_cpp(dat)
net_dp_cpp(dat)
```

### Exported functions

**Tree-based** (mode-centric traversal with memoized pruning):

- `tree_memo()` -- Pure R reference implementation
- `tree_memo_c()` -- C99 port with open-addressing hash tables
- `tree_memo_s4_c()` -- C99 with mode-path dominance pruning (S4)
- `tree_memo_cpp()` -- C++ port with `std::unordered_map`
- `tree_dp()` -- Pure R with constrained suffix-max DP
- `tree_dp_cpp()` -- C++ with constrained suffix-max DP

**Network-based** (state-space enumeration):

- `net_vander_cpp()` -- C++ with Vandermonde bulk-sum pruning
- `net_dp_cpp()` -- C++ with constrained suffix-max DP

**Diagnostics:**

- `tree_memo_profile()` -- C++ tree traversal returning a list
  with `$pvalue` and profiling counters (requires rebuilding
  with `PROFILE_V4` enabled in `src/Makevars`)

## Limitations

- Maximum 20 rows (`MAX_ROWS` in compiled code).
- Maximum total cell count of 10,000 (`MAX_N`).

## License

GPL (>= 3)
