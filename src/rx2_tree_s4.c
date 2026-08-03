// tree_s4.c
// Tree traversal + S4 mode-path dominance pruning.
// No memoization: find_max/find_min recompute on every call.
// Uses R_alloc for interrupt-safe memory management.
//
// S4: When traversing on the mode path and find_max of the
// mode-path subtree <= P_obs, all sibling subtrees at that
// level also have find_max <= P_obs (Mode-Path Dominance
// Theorem). The entire level can be pruned with a single
// find_max call instead of O(w) calls.
//
// The traverse function tracks whether the current node is
// on the mode path. When S3 fires on a mode-path node, S4
// signals the caller to skip all remaining siblings.

#include <R.h>
#include <Rinternals.h>
#include <math.h>
#include <string.h>

#define MAX_ROWS 20
#define FIND_MAX_STACK 4096

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

// S4 return code: traverse returns 1 if S4 cascade fires
#define S4_NONE 0
#define S4_CASCADE 1

typedef struct {
    int m;
    int n;
    int r[MAX_ROWS];
    int c[2];
    double tol;
    double log_const;
    double p_obs;
    double *cmax;
    int c0p1;
    int suffix_r[MAX_ROWS + 1];
    double *lfact;
    double pval;
    double p_max;
    double p_min;
    int y_max[MAX_ROWS];
    int mode_path[MAX_ROWS];
} FisherState;

static double compute_prob(const int *y,
                                const FisherState *s) {
    double log_prob = s->log_const;
    int i;
    for (i = 0; i < s->m; i++) {
        log_prob -= s->lfact[y[i]];
        log_prob -= s->lfact[s->r[i] - y[i]];
    }
    return exp(log_prob);
}

static void joe_min_impl(int k, int descend,
    int *r_avail, int c1, int n_rem, int *y,
    FisherState *s,
    double *local_min_p, int *local_min_y) {
    int i;
    if (k == s->m) {
        for (i = 0; i < s->m; i++) {
            if (r_avail[i] >= 0) y[i] = c1;
        }
        double prob = compute_prob(y, s);
        if (prob < *local_min_p) {
            *local_min_p = prob;
            memcpy(local_min_y, y, s->m * sizeof(int));
        }
        return;
    }

    int r_avail_save[MAX_ROWS], y_save[MAX_ROWS];
    memcpy(r_avail_save, r_avail, s->m * sizeof(int));
    memcpy(y_save, y, s->m * sizeof(int));
    int c1_save = c1, n_rem_save = n_rem;

    int idx = 0, max_r = -1;
    for (i = 0; i < s->m; i++) {
        if (r_avail[i] > max_r) {
            max_r = r_avail[i]; idx = i;
        }
    }

    if (descend)
        y[idx] = MIN(r_avail[idx], c1);
    else
        y[idx] = r_avail[idx] -
                 MIN(r_avail[idx], n_rem - c1);

    n_rem -= r_avail[idx];
    r_avail[idx] = -1;
    c1 -= y[idx];

    joe_min_impl(k + 1, 1, r_avail, c1, n_rem, y,
                     s, local_min_p, local_min_y);

    if (descend) {
        memcpy(r_avail, r_avail_save,
               s->m * sizeof(int));
        memcpy(y, y_save, s->m * sizeof(int));
        joe_min_impl(k, 0, r_avail, c1_save,
                         n_rem_save, y, s,
                         local_min_p, local_min_y);
    }
}

static void find_min(int c1, int *y, int k_start,
                          FisherState *s) {
    int i;
    int r_avail[MAX_ROWS];
    int n_rem = 0;
    for (i = 0; i < s->m; i++) {
        if (i < k_start) r_avail[i] = -1;
        else { r_avail[i] = s->r[i]; n_rem += s->r[i]; }
    }

    double local_min_p = 1e300;
    int local_min_y[MAX_ROWS];
    memcpy(local_min_y, y, s->m * sizeof(int));

    joe_min_impl(k_start, 1, r_avail, c1, n_rem, y, s,
                     &local_min_p, local_min_y);

    s->p_min = local_min_p;
}

typedef struct {
    int k, c1, n_rem;
    double log_p;
    int y[MAX_ROWS];
} FindMaxEntryC;

static void find_max(int k_start, int c1,
                          int *y, FisherState *s) {
    int n_rem = s->suffix_r[k_start];
    double local_max_p = 0;
    int local_max_y[MAX_ROWS];
    memcpy(local_max_y, y, s->m * sizeof(int));

    double log_p_prefix = s->log_const;
    int i;
    for (i = 0; i < k_start; i++)
        log_p_prefix -= s->lfact[y[i]]
                      + s->lfact[s->r[i] - y[i]];

    FindMaxEntryC stack[FIND_MAX_STACK];
    int sp = 0;

    stack[sp].k = k_start;
    stack[sp].c1 = c1;
    stack[sp].n_rem = n_rem;
    stack[sp].log_p = log_p_prefix;
    memcpy(stack[sp].y, y, s->m * sizeof(int));
    sp++;

    while (sp > 0) {
        sp--;
        FindMaxEntryC curr = stack[sp];
        int k = curr.k;
        int c1_rem = curr.c1;
        int n_rem_curr = curr.n_rem;

        if (k >= s->m) {
            double prob = exp(curr.log_p);
            if (prob > local_max_p) {
                local_max_p = prob;
                memcpy(local_max_y, curr.y,
                       s->m * sizeof(int));
            }
            continue;
        }

        int d_rem = s->m - k;
        int denom = n_rem_curr + d_rem;
        int y_lo, y_up;

        if (denom == 0 || c1_rem == 0) {
            y_lo = y_up = 0;
        } else {
            y_up = (int)floor(
                (double)(s->r[k] + 1) *
                (c1_rem + d_rem - 1) / denom);
            y_lo = (int)ceil(
                (double)(s->r[k] + 1) *
                (c1_rem + 1) / denom) - 1;
        }

        int n_after = s->suffix_r[k + 1];
        y_lo = MAX(y_lo, MAX(0, c1_rem - n_after));
        y_up = MIN(y_up, MIN(s->r[k], c1_rem));

        if (y_lo > y_up) continue;

        int y_k;
        for (y_k = y_lo; y_k <= y_up; y_k++) {
            if (sp >= FIND_MAX_STACK)
                Rf_error("find_max stack overflow "
                         "(capacity %d_rem)", FIND_MAX_STACK);
            stack[sp].k = k + 1;
            stack[sp].c1 = c1_rem - y_k;
            stack[sp].n_rem = n_rem_curr - s->r[k];
            stack[sp].log_p = curr.log_p
                - s->lfact[y_k]
                - s->lfact[s->r[k] - y_k];
            memcpy(stack[sp].y, curr.y,
                   s->m * sizeof(int));
            stack[sp].y[k] = y_k;
            sp++;
        }
    }

    s->p_max = local_max_p;
    memcpy(s->y_max, local_max_y, s->m * sizeof(int));
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

// Returns S4_CASCADE if mode-path S3 fired (caller should
// skip remaining siblings), S4_NONE otherwise.
static int traverse(int k, int c1, int *y,
    double prob_prefix, int n_rem,
    int descend, int y_k, int dir, int mode,
    int on_mode_path, FisherState *s) {

    double hp = 0.0;
    int prev_yk = -1;
    int need_full_dhyper = 1;

    for (;;) {
        if (k == s->m - 1) {
            y[s->m - 1] = c1;
            double hp_last = dhyper_lf(
                c1, c1, n_rem - c1,
                s->r[s->m - 1], s->lfact);
            double prob = prob_prefix * hp_last;
            if (prob > s->p_obs * (1 + s->tol))
                s->pval -= prob;
            return S4_NONE;
        }

        int y_lo = MAX(0, c1 - s->suffix_r[k + 1]);
        int y_hi = MIN(s->r[k], c1);

        if (descend) {
            mode = s->y_max[k];
            y_k = mode;
            need_full_dhyper = 1;
        }
        if (y_k < y_lo) return S4_NONE;
        if (y_k > y_hi) {
            y_k = mode - 1; dir = -1;
            descend = 0;
            need_full_dhyper = 1;
            continue;
        }

        y[k] = y_k;
        int c2 = n_rem - c1;

        if (need_full_dhyper) {
            hp = dhyper_lf(y_k, c1, c2,
                            s->r[k], s->lfact);
            need_full_dhyper = 0;
        } else {
            if (dir == +1) {
                hp = hp * (double)(c1 - prev_yk) *
                     (double)(s->r[k] - prev_yk) /
                     ((double)(prev_yk + 1) *
                      (double)(c2 - s->r[k] +
                               prev_yk + 1));
            } else {
                hp = hp * (double)prev_yk *
                     (double)(c2 - s->r[k] +
                              prev_yk) /
                     ((double)(c1 - prev_yk + 1) *
                      (double)(s->r[k] - prev_yk +
                               1));
            }
        }
        prev_yk = y_k;

        double new_prefix = prob_prefix * hp;
        int new_c1 = c1 - y_k;
        int new_n_rem = s->suffix_r[k + 1];

        int is_mode_node = on_mode_path &&
                           (y_k == s->mode_path[k]);

        if (s->m - k > 2) {
            // S1: constrained suffix-max bound
            if (new_prefix *
                s->cmax[(k + 1) * s->c0p1 + new_c1] <=
                s->p_obs * (1 + s->tol)) {
                on_mode_path = on_mode_path &&
                               !is_mode_node;
                y_k += dir; descend = 0;
                continue;
            }

            s->p_max = 0;
            find_max(k + 1, new_c1, y, s);
            s->p_min = 1;
            find_min(new_c1, y, k + 1, s);

            // S2: find_min bulk subtraction
            if (s->p_min > s->p_obs * (1 + s->tol)) {
                s->pval -= new_prefix;
                on_mode_path = on_mode_path &&
                               !is_mode_node;
                y_k += dir; descend = 0;
                continue;
            }

            // S3: find_max exclusion
            if (s->p_max <= s->p_obs * (1 + s->tol)) {
                // S4: mode-path dominance cascade
                if (is_mode_node)
                    return S4_CASCADE;
                on_mode_path = 0;
                y_k += dir; descend = 0;
                continue;
            }
        }

        traverse(k + 1, new_c1, y, new_prefix,
                     new_n_rem, 1, 0, +1, 0,
                     is_mode_node, s);

        on_mode_path = on_mode_path && !is_mode_node;
        y_k += dir; descend = 0;
    }
}

static double choose_round(int n_arg, int k_arg,
                            const double *lfact) {
    if (k_arg < 0 || k_arg > n_arg) return 0.0;
    if (k_arg == 0 || k_arg == n_arg) return 1.0;
    return floor(exp(lfact[n_arg] - lfact[k_arg]
                     - lfact[n_arg - k_arg]) + 0.5);
}

static void precompute_cmax(FisherState *s) {
    int c0p1 = s->c[0] + 1;
    int i, k, c1, y_k;
    s->c0p1 = c0p1;
    s->cmax = (double *)R_alloc(
        (s->m + 1) * c0p1, sizeof(double));
    memset(s->cmax, 0,
           (size_t)(s->m + 1) * c0p1 * sizeof(double));

    s->cmax[s->m * c0p1 + 0] = 1.0;

    for (k = s->m - 1; k >= 0; k--) {
        for (c1 = 0; c1 < c0p1; c1++) {
            int y_lo = MAX(0, c1 - s->suffix_r[k + 1]);
            int y_hi = MIN(s->r[k], c1);
            if (y_lo > y_hi) continue;

            double best = 0.0;
            double binom = choose_round(
                s->r[k], y_lo, s->lfact);
            int rk = s->r[k];

            for (y_k = y_lo; y_k <= y_hi; y_k++) {
                double val = binom *
                    s->cmax[(k + 1) * c0p1 +
                            (c1 - y_k)];
                if (val > best) best = val;
                if (y_k < y_hi)
                    binom = binom * (rk - y_k) /
                            (y_k + 1);
            }
            s->cmax[k * c0p1 + c1] = best;
        }
    }
}

double rx2_tree_s4_c_impl(int *dat, int m) {
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

    s.lfact = (double *)R_alloc(s.n + 1, sizeof(double));

    int flipped = s.c[0] > s.c[1];
    if (flipped) {
        int tmp = s.c[0];
        s.c[0] = s.c[1];
        s.c[1] = tmp;
    }

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

    precompute_cmax(&s);

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

    for (i = 0; i < m; i++) s.y_max[i] = 0;
    s.p_max = 0;
    s.p_min = 1;

    find_max(0, s.c[0], y_obs, &s);

    memcpy(s.mode_path, s.y_max, s.m * sizeof(int));

    traverse(0, s.c[0], y_obs, 1.0, s.n,
                 1, 0, +1, 0, 1, &s);

    return s.pval;
}
