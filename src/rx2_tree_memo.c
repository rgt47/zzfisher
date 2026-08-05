// tree_memo.c
// Pure C99 port of tree_memo algorithm.
// Open-addressing hash table with dynamically sized capacity.
// Uses R_alloc for interrupt-safe memory management.

#include <R.h>
#include <Rinternals.h>
#include <limits.h>
#include <math.h>
#include <string.h>

#define MAX_ROWS 20
#define FIND_MAX_STACK 4096

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

typedef struct {
    int y[MAX_ROWS];
} CacheEntry;

typedef struct {
    int key;
    int occupied;
    CacheEntry entry;
} HashBucket;

typedef struct {
    int m;
    int n;
    int r[MAX_ROWS];
    int c[2];
    int key_mult;
    int hash_cap;
    double tol;
    double log_const;
    double p_obs;
    int suffix_r[MAX_ROWS + 1];
    double *lfact;
    double pval;
    double p_max;
    double p_min;
    int y_max[MAX_ROWS];
    HashBucket *memo_max;
    HashBucket *memo_min;
} FisherState;

static int cache_key(int k, int c1, int mult) {
    return k * mult + c1;
}

static int hash_index(int key, int cap) {
    unsigned int h = (unsigned int)key;
    h = ((h >> 16) ^ h) * 0x45d9f3b;
    h = ((h >> 16) ^ h) * 0x45d9f3b;
    h = (h >> 16) ^ h;
    return (int)(h & (cap - 1));
}

static CacheEntry *cache_lookup(HashBucket *table,
                                  int key, int cap) {
    int idx = hash_index(key, cap);
    int i;
    for (i = 0; i < cap; i++) {
        int probe = (idx + i) & (cap - 1);
        if (!table[probe].occupied)
            return NULL;
        if (table[probe].key == key)
            return &table[probe].entry;
    }
    return NULL;
}

static void cache_insert(HashBucket *table, int key,
                          const int *y, int from, int to,
                          int cap) {
    int idx = hash_index(key, cap);
    int i;
    for (i = 0; i < cap; i++) {
        int probe = (idx + i) & (cap - 1);
        if (!table[probe].occupied ||
            table[probe].key == key) {
            table[probe].occupied = 1;
            table[probe].key = key;
            memcpy(table[probe].entry.y + from, y + from,
                   (to - from) * sizeof(int));
            return;
        }
    }
}

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
    int key = cache_key(k_start, c1, s->key_mult);
    int i;

    CacheEntry *cached = cache_lookup(s->memo_min, key,
                                        s->hash_cap);
    if (cached) {
        memcpy(y + k_start, cached->y + k_start,
               (s->m - k_start) * sizeof(int));
        s->p_min = compute_prob(y, s);
        return;
    }

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
    cache_insert(s->memo_min, key, local_min_y, k_start, s->m,
                 s->hash_cap);
}

typedef struct {
    int k, c1, n_rem;
    int y[MAX_ROWS];
    double log_p;
} FindMaxEntry;

static void find_max(int k_start, int c1,
                          int *y, FisherState *s) {
    int key = cache_key(k_start, c1, s->key_mult);

    CacheEntry *cached = cache_lookup(s->memo_max, key,
                                        s->hash_cap);
    if (cached) {
        memcpy(y + k_start, cached->y + k_start,
               (s->m - k_start) * sizeof(int));
        s->p_max = compute_prob(y, s);
        memcpy(s->y_max, y, s->m * sizeof(int));
        return;
    }

    int i;
    int n_rem = s->suffix_r[k_start];
    double local_max_p = 0;
    int local_max_y[MAX_ROWS];
    memcpy(local_max_y, y, s->m * sizeof(int));

    double log_p_prefix = s->log_const;
    for (i = 0; i < k_start; i++)
        log_p_prefix -= s->lfact[y[i]]
                      + s->lfact[s->r[i] - y[i]];

    FindMaxEntry stack[FIND_MAX_STACK];
    int sp = 0;

    stack[sp].k = k_start;
    stack[sp].c1 = c1;
    stack[sp].n_rem = n_rem;
    memcpy(stack[sp].y, y, s->m * sizeof(int));
    stack[sp].log_p = log_p_prefix;
    sp++;

    while (sp > 0) {
        sp--;
        FindMaxEntry curr = stack[sp];
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
                         "(capacity %d)", FIND_MAX_STACK);
            stack[sp].k = k + 1;
            stack[sp].c1 = c1_rem - y_k;
            stack[sp].n_rem = n_rem_curr - s->r[k];
            memcpy(stack[sp].y, curr.y,
                   s->m * sizeof(int));
            stack[sp].y[k] = y_k;
            stack[sp].log_p = curr.log_p
                - s->lfact[y_k]
                - s->lfact[s->r[k] - y_k];
            sp++;
        }
    }

    s->p_max = local_max_p;
    memcpy(s->y_max, local_max_y, s->m * sizeof(int));

    cache_insert(s->memo_max, key, local_max_y,
                 k_start, s->m, s->hash_cap);
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

static void penult(int c1, int n_rem,
    double prob_prefix, FisherState *s) {
    int k = s->m - 2;
    int y_lo = MAX(0, c1 - s->r[s->m - 1]);
    int y_hi = MIN(s->r[k], c1);
    if (y_lo > y_hi) return;

    int c2 = n_rem - c1;
    double thresh = s->p_obs * (1 + s->tol);

    int mode_k = s->y_max[k];
    if (mode_k < y_lo) mode_k = y_lo;
    if (mode_k > y_hi) mode_k = y_hi;

    double hp_mode = dhyper_lf(mode_k, c1, c2,
                                 s->r[k], s->lfact);
    double prob = prob_prefix * hp_mode;
    if (prob > thresh)
        s->pval -= prob;

    int prev_yk;
    double hp;
    int y_k;

    prev_yk = mode_k;
    hp = hp_mode;
    for (y_k = mode_k + 1; y_k <= y_hi; y_k++) {
        hp *= (double)(c1 - prev_yk) *
              (double)(s->r[k] - prev_yk) /
              ((double)(prev_yk + 1) *
               (double)(c2 - s->r[k] + prev_yk + 1));
        prev_yk = y_k;
        prob = prob_prefix * hp;
        if (prob > thresh)
            s->pval -= prob;
    }

    prev_yk = mode_k;
    hp = hp_mode;
    for (y_k = mode_k - 1; y_k >= y_lo; y_k--) {
        hp *= (double)prev_yk *
              (double)(c2 - s->r[k] + prev_yk) /
              ((double)(c1 - prev_yk + 1) *
               (double)(s->r[k] - prev_yk + 1));
        prev_yk = y_k;
        prob = prob_prefix * hp;
        if (prob > thresh)
            s->pval -= prob;
    }
}

static void traverse(int k, int c1, int *y,
    double prob_prefix, int n_rem, int descend,
    int y_k, int dir, int mode, FisherState *s) {

    double hp = 0.0;
    int prev_yk = -1;
    int need_full = 1;

    for (;;) {
        if (k >= s->m - 1) {
            if (k == s->m - 1) {
                y[s->m - 1] = c1;
                double hp_last = dhyper_lf(
                    c1, c1, n_rem - c1,
                    s->r[s->m - 1], s->lfact);
                double prob = prob_prefix * hp_last;
                if (prob > s->p_obs * (1 + s->tol))
                    s->pval -= prob;
            }
            return;
        }

        int y_lo = MAX(0, c1 - s->suffix_r[k + 1]);
        int y_hi = MIN(s->r[k], c1);

        if (descend) {
            mode = s->y_max[k]; y_k = mode;
            need_full = 1;
        }
        if (y_k < y_lo) return;
        if (y_k > y_hi) {
            y_k = mode - 1; dir = -1;
            descend = 0;
            need_full = 1;
            continue;
        }

        y[k] = y_k;
        int c2 = n_rem - c1;
        if (need_full) {
            hp = dhyper_lf(y_k, c1, c2,
                            s->r[k], s->lfact);
            need_full = 0;
        } else if (dir == +1) {
            hp *= (double)(c1 - prev_yk) *
                  (double)(s->r[k] - prev_yk) /
                  ((double)(prev_yk + 1) *
                   (double)(c2 - s->r[k] + prev_yk + 1));
        } else {
            hp *= (double)prev_yk *
                  (double)(c2 - s->r[k] + prev_yk) /
                  ((double)(c1 - prev_yk + 1) *
                   (double)(s->r[k] - prev_yk + 1));
        }
        prev_yk = y_k;

        double new_prefix = prob_prefix * hp;
        int new_c1 = c1 - y_k;
        int new_n_rem = s->suffix_r[k + 1];

        if (s->m - k > 2) {
            s->p_max = 0;
            find_max(k + 1, new_c1, y, s);
            if (s->p_max <= s->p_obs * (1 + s->tol)) {
                y_k += dir; descend = 0;
                continue;
            }
            s->p_min = 1;
            find_min(new_c1, y, k + 1, s);
            if (s->p_min > s->p_obs * (1 + s->tol)) {
                s->pval -= new_prefix;
                y_k += dir; descend = 0;
                continue;
            }
        }

        if (k == s->m - 3) {
            penult(new_c1, new_n_rem, new_prefix, s);
        } else {
            traverse(k + 1, new_c1, y, new_prefix,
                         new_n_rem, 1, 0, +1, 0, s);
        }
        y_k += dir; descend = 0;
    }
}

static int next_power_of_two(int v) {
    v--;
    v |= v >> 1; v |= v >> 2;
    v |= v >> 4; v |= v >> 8;
    v |= v >> 16;
    return v + 1;
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

    s.key_mult = s.c[0] + 1;
    size_t n_states = (size_t)(m + 1) * (size_t)s.key_mult;
    if (n_states > (size_t)INT_MAX / 2)
        Rf_error("Table too large for hash table "
                 "(n_states = %zu)", n_states);
    s.hash_cap = next_power_of_two(
        (int)(n_states < 64 ? 64 : 2 * n_states));

    s.lfact = (double *)R_alloc(s.n + 1, sizeof(double));
    s.memo_max = (HashBucket *)R_alloc(s.hash_cap,
        sizeof(HashBucket));
    s.memo_min = (HashBucket *)R_alloc(s.hash_cap,
        sizeof(HashBucket));
    memset(s.memo_max, 0,
           s.hash_cap * sizeof(HashBucket));
    memset(s.memo_min, 0,
           s.hash_cap * sizeof(HashBucket));

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

    for (i = 0; i < m; i++) s.y_max[i] = 0;
    s.p_max = 0;
    s.p_min = 1;

    find_max(0, s.c[0], y_obs, &s);
    traverse(0, s.c[0], y_obs, 1.0, s.n,
                 1, 0, +1, 0, &s);

    return s.pval;
}
