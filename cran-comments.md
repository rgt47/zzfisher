## Submission summary

This is a new submission. zzfisher implements Fisher's exact test for r x 2 tables by several algorithms, with exact power and sample size.

This package contains compiled C and C++ code. Before submission it
should also be checked under the CRAN sanitizer builds (ASan/UBSan),
which has not yet been done on this machine.

## Test environments

* Local: macOS 26.6.2 (aarch64-apple-darwin25.4.0), R 4.6.1

Before submission this should also be checked on:

* win-builder (devel and release)
* macbuilder
* R-hub (a Linux and a Windows target)

## R CMD check results

0 errors | 0 warnings | 1 note

The note is:

```
* checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Ronald G. Thomas <rgthomas@ucsd.edu>'
  New submission
```

This note is expected for a first-time submission.

## Reverse dependencies

None. This is a new package with no reverse dependencies.

