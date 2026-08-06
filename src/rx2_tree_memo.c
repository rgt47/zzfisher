// rx2_tree_memo.c
// Pure C99 tree-traversal kernel for Fisher's exact test on m x 2
// tables, redesigned around the c = 2 separability of the
// objective (docs/tree_memo_efficiency_review_2026-08-04.md).
//
// Because the table has two columns, the log-probability separates
// as log P(y) = log_const + sum_i h_i(y_i) with
//   h_i(y) = -lfact[y] - lfact[r_i - y]
// and each h_i strictly concave. Consequences used here:
//
//   1. Exact per-subtree bounds by dynamic programming. maxh[k][a]
//      (resp. minh[k][a]) is the exact maximum (minimum) of
//      sum_{i >= k} h_i(y_i) over completions with sum a. These
//      replace the former find_max / find_min searches and their
//      memoization hash tables; every bound query is one array
//      read. The maximum row is Requena and Martin Ciudad (2006),
//      Theorem 1, in dynamic-programming form.
//
//   2. Proven two-sided sibling cutoff. The child-attention value
//      F(y) = h_k(y) + maxh[k+1][a - y] is a sum of two concave
//      functions of y (concavity of the value function a -> maxh
//      is standard for separable concave maximization), hence
//      concave: the children whose subtrees are not entirely
//      inside the significance region S = {P <= P_obs} form a
//      contiguous interval (equivalently: the children whose
//      subtrees meet the complement region U = {P > P_obs},
//      which is convex, form the shadow of a convex set).
//      The walk locates the peak of F by local ascent and stops in
//      each direction at the first child with
//      pref_lp + F(y) <= log_thresh. This is the c = 2 proof of
//      the cascade that rx2_tree_s4.c conjectures (Requena and
//      Martin Ciudad 2006, Theorem 2).
//
//   3. penult tail cutoff. At the next-to-last level the tested
//      quantity is the exact one-dimensional hypergeometric pmf,
//      which is unimodal, so the leaf walk stops in each direction
//      at the first value not exceeding the threshold.
//
// All bound comparisons are performed in log space
// (log_thresh = log p_obs + log1p(tol)); the subtraction of
// complement mass remains in probability space via the conditional
// hypergeometric factorization (the dhyper factors of a subtree's
// completions sum to one, so the whole subtree's mass equals the
// prefix product).
//
// Row sort ascending by margin and column flip to the smaller
// first column are unchanged from the previous design.

#include <R.h>
#include <Rinternals.h>
#include <limits.h>
#include <math.h>
#include <string.h>

#define MAX_ROWS 20

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

typedef struct {
    int m;
    int n;
    int r[MAX_ROWS];
    int c[2];
    double tol;
    double log_const;
    double p_obs;
    double log_thresh;
    int suffix_r[MAX_ROWS + 1];
    double *lfact;
    double *maxh;   /* (m+1) x (c1+1), maxh[k*(c1+1) + a] */
    double *minh;
    int acap;       /* c1 + 1 */
    double pval;
} FisherState;

/* h_k(y) = -lfact[y] - lfact[r_k - y]; caller guarantees bounds. */
static double h_term(const FisherState *s, int k, int y) {
    return -s->lfact[y] - s->lfact[s->r[k] - y];
}

static double maxh_at(const FisherState *s, int k, int a) {
    return s->maxh[k * s->acap + a];
}
static double minh_at(const FisherState *s, int k, int a) {
    return s->minh[k * s->acap + a];
}

/* Exact suffix extrema of sum h_i under a budget, by backward DP.
 * maxh[m][0] = 0 and maxh[m][a>0] = -Inf; likewise minh with +Inf.
 * Feasible budgets at level k are 0..min(c1, suffix_r[k]). */
static void build_extrema(FisherState *s) {
    int m = s->m, acap = s->acap;
    int k, a, y;

    for (a = 0; a < acap; a++) {
        s->maxh[m * acap + a] = (a == 0) ? 0.0 : -INFINITY;
        s->minh[m * acap + a] = (a == 0) ? 0.0 : INFINITY;
    }
    for (k = m - 1; k >= 0; k--) {
        int amax = MIN(acap - 1, s->suffix_r[k]);
        for (a = 0; a < acap; a++) {
            double best_max = -INFINITY, best_min = INFINITY;
            if (a <= amax) {
                int y_lo = MAX(0, a - s->suffix_r[k + 1]);
                int y_hi = MIN(s->r[k], a);
                for (y = y_lo; y <= y_hi; y++) {
                    double nx = s->maxh[(k + 1) * acap + (a - y)];
                    double nn = s->minh[(k + 1) * acap + (a - y)];
                    double h = h_term(s, k, y);
                    if (isfinite(nx) && h + nx > best_max)
                        best_max = h + nx;
                    if (isfinite(nn) && h + nn < best_min)
                        best_min = h + nn;
                }
            }
            s->maxh[k * acap + a] = best_max;
            s->minh[k * acap + a] = best_min;
        }
    }
}

/* Hypergeometric pmf. Parameter names deliberately avoid the CM
 * symbols m (rows), n (grand total), and k (traversal level):
 * m_white/n_black/k_draw are the white-ball, black-ball, and draw
 * counts of the classical urn parameterization. */
static double dhyper_lf(int x, int m_white, int n_black, int k_draw,
                          const double *lfact) {
    if (x < 0 || x > k_draw || x > m_white || k_draw - x > n_black)
        return 0.0;
    double lp = lfact[m_white] - lfact[x] - lfact[m_white - x]
              + lfact[n_black] - lfact[k_draw - x]
              - lfact[n_black - k_draw + x]
              - lfact[m_white + n_black] + lfact[k_draw]
              + lfact[m_white + n_black - k_draw];
    return exp(lp);
}

/* Ratio steps for the dhyper recurrence along y -> y + 1. */
static double hp_step_up(double hp, int prev, int c1, int c2, int rk) {
    return hp * (double)(c1 - prev) * (double)(rk - prev) /
           ((double)(prev + 1) * (double)(c2 - rk + prev + 1));
}
static double hp_step_down(double hp, int prev, int c1, int c2, int rk) {
    return hp * (double)prev * (double)(c2 - rk + prev) /
           ((double)(c1 - prev + 1) * (double)(rk - prev + 1));
}

/* Next-to-last level: y at level m-2 fixes the leaf. The tested
 * quantity is prob_prefix times the exact 1-D hypergeometric pmf,
 * unimodal in y, so the walk breaks at the first value at or below
 * the threshold in each direction. */
static void penult(int a, int n_rem, double prob_prefix,
                   FisherState *s) {
    int k = s->m - 2;
    int y_lo = MAX(0, a - s->r[s->m - 1]);
    int y_hi = MIN(s->r[k], a);
    if (y_lo > y_hi) return;

    int c2 = n_rem - a;
    double thresh = s->p_obs * (1 + s->tol);

    /* Locate the pmf mode by the closed form, clamped, then refined
     * by local ascent (guards the floor formula's edge cases). */
    int mode_k = (int)((double)(s->r[k] + 1) * (a + 1) /
                       (n_rem + 2));
    if (mode_k < y_lo) mode_k = y_lo;
    if (mode_k > y_hi) mode_k = y_hi;
    double hp_mode = dhyper_lf(mode_k, a, c2, s->r[k], s->lfact);
    while (mode_k < y_hi) {
        double up = dhyper_lf(mode_k + 1, a, c2, s->r[k], s->lfact);
        if (up <= hp_mode) break;
        mode_k++; hp_mode = up;
    }
    while (mode_k > y_lo) {
        double dn = dhyper_lf(mode_k - 1, a, c2, s->r[k], s->lfact);
        if (dn <= hp_mode) break;
        mode_k--; hp_mode = dn;
    }

    double prob = prob_prefix * hp_mode;
    if (prob > thresh)
        s->pval -= prob;
    else
        return;              /* peak fails: nothing qualifies */

    int prev, y_k;
    double hp;

    prev = mode_k; hp = hp_mode;
    for (y_k = mode_k + 1; y_k <= y_hi; y_k++) {
        hp = hp_step_up(hp, prev, a, c2, s->r[k]);
        prev = y_k;
        prob = prob_prefix * hp;
        if (prob > thresh) s->pval -= prob;
        else break;          /* unimodal: all further fail */
    }
    prev = mode_k; hp = hp_mode;
    for (y_k = mode_k - 1; y_k >= y_lo; y_k--) {
        hp = hp_step_down(hp, prev, a, c2, s->r[k]);
        prev = y_k;
        prob = prob_prefix * hp;
        if (prob > thresh) s->pval -= prob;
        else break;
    }
}

/* pref_lp + F(y) with F(y) = h_k(y) + maxh[k+1][a - y]; -Inf when
 * the child is infeasible. */
static double child_max_lp(const FisherState *s, int k, int a,
                           double pref_lp, int y) {
    double nx = maxh_at(s, k + 1, a - y);
    if (!isfinite(nx)) return -INFINITY;
    return pref_lp + h_term(s, k, y) + nx;
}

static void traverse(int k, int a, double prob_prefix,
                     double pref_lp, int n_rem, FisherState *s);

/* Handle one child: bulk-subtract if its subtree lies entirely in
 * the complement, otherwise descend. The caller has already
 * established that the child is not entirely in the acceptance
 * region. */
static void handle_child(int k, int a, int y_k, double hp,
                         double prob_prefix, double pref_lp,
                         FisherState *s) {
    double new_prefix = prob_prefix * hp;
    int new_a = a - y_k;
    double nn = minh_at(s, k + 1, new_a);
    if (isfinite(nn) &&
        pref_lp + h_term(s, k, y_k) + nn > s->log_thresh) {
        s->pval -= new_prefix;   /* whole subtree in complement */
        return;
    }
    double new_pref_lp = pref_lp + h_term(s, k, y_k);
    if (k == s->m - 3)
        penult(new_a, s->suffix_r[k + 1], new_prefix, s);
    else
        traverse(k + 1, new_a, new_prefix, new_pref_lp,
                 s->suffix_r[k + 1], s);
}

static void traverse(int k, int a, double prob_prefix,
                     double pref_lp, int n_rem, FisherState *s) {
    if (k >= s->m - 1) {
        if (k == s->m - 1) {
            double hp_last = dhyper_lf(a, a, n_rem - a,
                                       s->r[s->m - 1], s->lfact);
            double prob = prob_prefix * hp_last;
            if (prob > s->p_obs * (1 + s->tol))
                s->pval -= prob;
        }
        return;
    }

    int y_lo = MAX(0, a - s->suffix_r[k + 1]);
    int y_hi = MIN(s->r[k], a);
    if (y_lo > y_hi) return;

    /* Locate the peak of the concave attention value F by local
     * ascent from the 1-D mode guess. */
    int y_pk = (int)((double)(s->r[k] + 1) * (a + 1) / (n_rem + 2));
    if (y_pk < y_lo) y_pk = y_lo;
    if (y_pk > y_hi) y_pk = y_hi;
    {
        double f = child_max_lp(s, k, a, pref_lp, y_pk);
        while (y_pk < y_hi) {
            double up = child_max_lp(s, k, a, pref_lp, y_pk + 1);
            if (up <= f) break;
            y_pk++; f = up;
        }
        while (y_pk > y_lo) {
            double dn = child_max_lp(s, k, a, pref_lp, y_pk - 1);
            if (dn <= f) break;
            y_pk--; f = dn;
        }
        if (f <= s->log_thresh)
            return;   /* peak fails: every child is in S */
    }

    int c2 = n_rem - a;
    int prev, y_k;
    double hp;

    /* Upward from the peak; concavity of F licenses the break. */
    hp = dhyper_lf(y_pk, a, c2, s->r[k], s->lfact);
    handle_child(k, a, y_pk, hp, prob_prefix, pref_lp, s);
    prev = y_pk;
    for (y_k = y_pk + 1; y_k <= y_hi; y_k++) {
        if (child_max_lp(s, k, a, pref_lp, y_k) <= s->log_thresh)
            break;
        hp = hp_step_up(hp, prev, a, c2, s->r[k]);
        prev = y_k;
        handle_child(k, a, y_k, hp, prob_prefix, pref_lp, s);
    }
    /* Downward from the peak. */
    hp = dhyper_lf(y_pk, a, c2, s->r[k], s->lfact);
    prev = y_pk;
    for (y_k = y_pk - 1; y_k >= y_lo; y_k--) {
        if (child_max_lp(s, k, a, pref_lp, y_k) <= s->log_thresh)
            break;
        hp = hp_step_down(hp, prev, a, c2, s->r[k]);
        prev = y_k;
        handle_child(k, a, y_k, hp, prob_prefix, pref_lp, s);
    }
}

double rx2_tree_memo_c_impl(int *dat, int m) {
    if (m < 2) return 1.0;
    if (m > MAX_ROWS)
        Rf_error("Too many rows (max %d)", MAX_ROWS);

    FisherState s;
    memset(&s, 0, sizeof(s));
    s.m = m;
    s.tol = 3.45254e-7;
    s.pval = 1.0;

    int r_unsorted[MAX_ROWS];
    int i;
    for (i = 0; i < m; i++) {
        r_unsorted[i] = dat[i] + dat[i + m];
        s.n += r_unsorted[i];
        s.c[0] += dat[i];
        s.c[1] += dat[i + m];
    }

    int flipped = s.c[0] > s.c[1];
    if (flipped) {
        int tmp = s.c[0];
        s.c[0] = s.c[1];
        s.c[1] = tmp;
    }

    s.acap = s.c[0] + 1;
    if ((size_t)(m + 1) * (size_t)s.acap >
        (size_t)INT_MAX / (2 * sizeof(double)))
        Rf_error("Table too large for bound arrays");

    s.lfact = (double *)R_alloc(s.n + 1, sizeof(double));
    s.maxh = (double *)R_alloc((size_t)(m + 1) * s.acap,
                               sizeof(double));
    s.minh = (double *)R_alloc((size_t)(m + 1) * s.acap,
                               sizeof(double));

    int order[MAX_ROWS];
    for (i = 0; i < m; i++) order[i] = i;
    int j;
    for (i = 1; i < m; i++) {
        int tmp_ord = order[i];
        int tmp_val = r_unsorted[tmp_ord];
        j = i - 1;
        while (j >= 0 &&
               r_unsorted[order[j]] > tmp_val) {
            order[j + 1] = order[j];
            j--;
        }
        order[j + 1] = tmp_ord;
    }

    int y_obs[MAX_ROWS];
    for (i = 0; i < m; i++) {
        int orig = order[i];
        s.r[i] = r_unsorted[orig];
        y_obs[i] = flipped ? dat[orig + m] : dat[orig];
    }

    s.suffix_r[m] = 0;
    for (i = m - 1; i >= 0; i--)
        s.suffix_r[i] = s.suffix_r[i + 1] + s.r[i];

    s.lfact[0] = 0.0;
    for (i = 1; i <= s.n; i++)
        s.lfact[i] = s.lfact[i - 1] + log((double)i);

    s.log_const = 0.0;
    for (i = 0; i < m; i++)
        s.log_const += s.lfact[s.r[i]];
    s.log_const += s.lfact[s.c[0]] +
                   s.lfact[s.c[1]] -
                   s.lfact[s.n];

    double log_p_obs = s.log_const;
    for (i = 0; i < m; i++) {
        log_p_obs -= s.lfact[y_obs[i]];
        log_p_obs -= s.lfact[s.r[i] - y_obs[i]];
    }
    s.p_obs = exp(log_p_obs);
    s.log_thresh = log_p_obs + log1p(s.tol);

    build_extrema(&s);

    if (m == 2) {
        penult(s.c[0], s.n, 1.0, &s);
    } else {
        traverse(0, s.c[0], 1.0, s.log_const, s.n, &s);
    }

    return s.pval;
}
