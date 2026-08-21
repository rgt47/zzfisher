// fxpower_rcpp.cpp -- Rcpp wrapper around fxpower.c
//
// Provides .Call() interface to the C power function, avoiding
// the copy overhead of .C().  Also exposes individual components
// (cv, xfmax, binom_joint) for testing and paper illustrations.

#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <algorithm>

using namespace Rcpp;

// ------------------------------------------------------------------ //
// Log-factorial table                                                  //
// ------------------------------------------------------------------ //
static std::vector<double> make_lfact(int nmax) {
  std::vector<double> lf(nmax + 1);
  lf[0] = 0.0;
  for (int i = 1; i <= nmax; i++) {
    lf[i] = lf[i - 1] + std::log(static_cast<double>(i));
  }
  return lf;
}

// ------------------------------------------------------------------ //
// Joint binomial density                                               //
// ------------------------------------------------------------------ //
static double binom_joint(const std::vector<double> &lf,
                          int n1, int n2,
                          double lp1, double lp2,
                          double lmp1, double lmp2,
                          int x, int r) {
  int y = r - x;
  double logf = lf[n1] - lf[x] - lf[n1 - x] +
                lf[n2] - lf[y] - lf[n2 - y] +
                x * lp1 + (n1 - x) * lmp1 +
                y * lp2 + (n2 - y) * lmp2;
  return std::exp(logf);
}

// ------------------------------------------------------------------ //
// Mode of f(x, r) for given r (generalized quadratic)                 //
// ------------------------------------------------------------------ //
static int xfmax_unequal(int n1, int n2, int r, double theta) {
  int xmin_r = std::max(0, r - n2);
  int xmax_r = std::min(r, n1);
  int xm;

  if (std::fabs(theta - 1.0) > 1e-15) {
    double a = theta - 1.0;
    double b = -(static_cast<double>(n1 + r + 2) * theta +
                 static_cast<double>(n2 - r));
    double c = static_cast<double>(n1 + 1) *
               static_cast<double>(r + 1) * theta;
    double disc = b * b - 4.0 * a * c;
    if (disc < 0.0) disc = 0.0;
    double q = -0.5 * (b - std::sqrt(disc));
    xm = static_cast<int>(c / q);
  } else {
    double mode = static_cast<double>(r + 1) *
                  static_cast<double>(n1 + 1) /
                  static_cast<double>(n1 + n2 + 2);
    xm = static_cast<int>(mode);
  }

  xm = std::max(xm, xmin_r);
  xm = std::min(xm, xmax_r);
  return xm;
}

// ------------------------------------------------------------------ //
// Critical value for one-sided Fisher test                             //
// ------------------------------------------------------------------ //
static int cv_unequal(const std::vector<double> &lf,
                      int n1, int n2, int r, double alpha) {
  int ntot = n1 + n2;
  int xmin_r = std::max(0, r - n2);
  int xmax_r = std::min(r, n1);

  if (xmax_r < xmin_r) return -1;

  int y0 = r - xmin_r;
  double logp = lf[n1] - lf[xmin_r] - lf[n1 - xmin_r] +
                lf[n2] - lf[y0] - lf[n2 - y0] -
                lf[ntot] + lf[r] + lf[ntot - r];
  double prob = std::exp(logp);

  double cumprob = 0.0;
  int xcv = -1;

  for (int x = xmin_r; x <= xmax_r; x++) {
    cumprob += prob;
    if (cumprob <= alpha) {
      xcv = x;
    } else {
      break;
    }
    if (x < xmax_r) {
      prob *= static_cast<double>(n1 - x) *
              static_cast<double>(r - x) /
              (static_cast<double>(x + 1) *
               static_cast<double>(n2 - r + x + 1));
    }
  }

  return xcv;
}

// ------------------------------------------------------------------ //
// Sum f(x, r-x) along fixed r via adjacency recurrence                //
// ------------------------------------------------------------------ //
static double fsum_unequal(int from, int to, int inc, double f,
                           double flim, int n1, int n2,
                           double theta, int r, int &rcnt) {
  double fp = 0.0;
  rcnt = 0;
  int x = from;

  while (true) {
    if (f < flim) break;
    rcnt++;
    fp += f;

    int xnext = x + inc;
    if (inc > 0 && xnext > to) break;
    if (inc < 0 && xnext < to) break;

    if (inc > 0) {
      f *= theta * static_cast<double>(n1 - x) *
           static_cast<double>(r - x) /
           (static_cast<double>(x + 1) *
            static_cast<double>(n2 - r + x + 1));
    } else {
      f *= (1.0 / theta) * static_cast<double>(x) *
           static_cast<double>(n2 - r + x) /
           (static_cast<double>(n1 - x + 1) *
            static_cast<double>(r - x + 1));
    }
    x = xnext;
  }

  return fp;
}

// ------------------------------------------------------------------ //
// Main power function -- exported to R via .Call()                     //
// ------------------------------------------------------------------ //

// [[Rcpp::export(name = ".fxpower_cpp")]]
List fxpower_cpp(int n1, int n2, double p1, double p2,
                 double alpha, double eps) {
  if (n1 <= 0 || n2 <= 0 ||
      p1 <= 0.0 || p1 >= 1.0 ||
      p2 <= 0.0 || p2 >= 1.0 ||
      alpha <= 0.0 || alpha >= 1.0 ||
      eps < 0.0 || eps >= 1.0) {
    stop("Invalid arguments to fxpower_cpp");
  }

  int ntot = n1 + n2;
  std::vector<double> lf = make_lfact(ntot);

  double lp1  = std::log(p1);
  double lp2  = std::log(p2);
  double lmp1 = std::log(1.0 - p1);
  double lmp2 = std::log(1.0 - p2);
  double theta = p1 * (1.0 - p2) / (p2 * (1.0 - p1));

  int xmode = static_cast<int>((n1 + 1.0) * p1);
  int ymode = static_cast<int>((n2 + 1.0) * p2);
  int rmax  = xmode + ymode;

  double fmax = binom_joint(lf, n1, n2, lp1, lp2, lmp1, lmp2,
                            xmode, rmax);
  double flim = (eps > 0.0) ? fmax * eps : 0.0;

  double powsum = 0.0;
  int count = 0;
  bool finished = false;

  int r = rmax;
  int i = 0;
  int j = 1;

  while (!finished) {
    if (r >= 0 && r <= ntot) {
      int xmin_r = std::max(0, r - n2);
      int xcv = cv_unequal(lf, n1, n2, r, alpha);
      int xfm = xfmax_unequal(n1, n2, r, theta);

      if (xcv >= xmin_r) {
        int rcnt = 0;
        double f;

        if (xcv <= xfm) {
          f = binom_joint(lf, n1, n2, lp1, lp2, lmp1, lmp2,
                          xcv, r);
          powsum += fsum_unequal(xcv, xmin_r, -1, f, flim,
                                 n1, n2, theta, r, rcnt);
          if (rcnt == 0 && i > 0) { finished = true; break; }
          count += rcnt;
        } else {
          f = binom_joint(lf, n1, n2, lp1, lp2, lmp1, lmp2,
                          xfm, r);
          powsum += fsum_unequal(xfm, xmin_r, -1, f, flim,
                                 n1, n2, theta, r, rcnt);
          if (rcnt == 0 && i > 0) { finished = true; break; }
          count += rcnt;

          f = binom_joint(lf, n1, n2, lp1, lp2, lmp1, lmp2,
                          xfm + 1, r);
          powsum += fsum_unequal(xfm + 1, xcv, 1, f, flim,
                                 n1, n2, theta, r, rcnt);
          count += rcnt;
        }
      }

      if (i % 100 == 0) Rcpp::checkUserInterrupt();
    }

    i++;
    j = -j;
    r = r + i * j;
    if (r < 0 || r > ntot) {
      i++;
      j = -j;
      r = r + i * j;
      if (r < 0 || r > ntot) break;
    }
    if (i > 2 * ntot) break;
  }

  // Trimming-error bound, following the original Algorithm AS 280
  // ERROR calculation (docs/fxpower_original.f, lines 204-208):
  // the omitted probability mass is bounded above by flim times
  // the number of support points never evaluated, and is also
  // bounded above by 1 - powsum (the trivial complement bound).
  // "total" below is the exact count of (x, r) support points
  // over the full asymmetric hypergeometric support, computed in
  // closed form rather than the heuristic estimate used by the
  // 1993 balanced-case FORTRAN, since n1 may differ from n2.
  double error_bound = 0.0;
  if (eps > 0.0) {
    long long total = 0;
    for (int rr = 0; rr <= ntot; rr++) {
      int xmin_rr = std::max(0, rr - n2);
      int xmax_rr = std::min(rr, n1);
      total += (xmax_rr - xmin_rr + 1);
    }
    double flim_bound = flim * static_cast<double>(total - count);
    double complement_bound = 1.0 - powsum;
    error_bound = std::min(flim_bound, complement_bound);
    if (error_bound < 0.0) error_bound = 0.0;
  }

  return List::create(
    Named("power") = powsum,
    Named("error") = error_bound,
    Named("count") = count
  );
}

// ------------------------------------------------------------------ //
// Exposed helpers for testing / paper illustrations                    //
// ------------------------------------------------------------------ //

// [[Rcpp::export(name = ".cv_unequal_cpp")]]
int cv_unequal_r(int n1, int n2, int r, double alpha) {
  int ntot = n1 + n2;
  std::vector<double> lf = make_lfact(ntot);
  return cv_unequal(lf, n1, n2, r, alpha);
}

// [[Rcpp::export(name = ".xfmax_unequal_cpp")]]
int xfmax_unequal_r(int n1, int n2, int r, double theta) {
  return xfmax_unequal(n1, n2, r, theta);
}

// [[Rcpp::export(name = ".binom_joint_cpp")]]
double binom_joint_r(int n1, int n2, double p1, double p2,
                     int x, int r) {
  int ntot = n1 + n2;
  std::vector<double> lf = make_lfact(ntot);
  double lp1  = std::log(p1);
  double lp2  = std::log(p2);
  double lmp1 = std::log(1.0 - p1);
  double lmp2 = std::log(1.0 - p2);
  return binom_joint(lf, n1, n2, lp1, lp2, lmp1, lmp2, x, r);
}
