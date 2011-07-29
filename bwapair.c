#include "bwapair.h"

#include "ksort.h"
#include <math.h>

typedef struct {
    int o_n;
    int subo_n;
    int cnt_chg;
    int max_len;
    /* TODO: make these into arrays of pointers */
    position_t last_pos[2][2];
    position_t o_pos[2];
    uint64_t subo_score, o_score;
} pairing_internals_t;

extern int g_log_n[256]; // in bwase.c

static int position_lt(const position_t a, const position_t b) {
    return a.remapped_pos < b.remapped_pos;
}
KSORT_INIT(position, position_t, position_lt);

static inline uint64_t hash_64(uint64_t key) {
    key += ~(key << 32);
    key ^= (key >> 22);
    key += ~(key << 13);
    key ^= (key >> 8);
    key += (key << 3);
    key ^= (key >> 15);
    key += ~(key << 27);
    key ^= (key >> 31);
    return key;
}

static inline int mappings_overlap(const position_t *a, const position_t* b, const alngrp_t aln[2]) {
    const alignment_t *aln1;
    const alignment_t *aln2;

    if (a->pos == (uint64_t)-1 || b->pos == (uint64_t)-1)
        return 0;

    aln1 = &__aln(*a, aln);
    aln2 = &__aln(*b, aln);

    if ((a->remapped_pos == b->remapped_pos) &&
        (a->idx_and_end&1) == (b->idx_and_end&1)
    ) {
        return 1;
    }
    return 0;
}

static inline position_t *select_mapping(const alngrp_t aln[2], const pos_arr_t *arr, int begin, int end, int *n_optimal) {
    int i;
    position_t *best = &arr->a[begin];

    for (i = begin+1; i <= end; ++i) {
        position_t *p = &arr->a[i];
        int this_score = __aln(*p, aln).aln.score;
        int best_score = __aln(*best, aln).aln.score;
        if (__aln(*p, aln).aln.score < __aln(*best, aln).aln.score) {
            *n_optimal = 1;
            best = p;
        } else if (this_score == best_score) {
            ++(*n_optimal);
            /* randomly pick one here? */
        }
    }

    return best;
}

static void pairing_aux(const pairing_param_t *param, pairing_internals_t *pint, const position_t *u, const position_t *v, int n_optimal) {
    bwa_seq_t **p = param->p;
    const alngrp_t *aln = param->aln;
    const pe_opt_t *opt = param->opt;
    const isize_info_t *ii = param->ii;

    // here v>=u. When ii is set, we check insert size with ii; otherwise with opt->max_isize
    bwtint_t l = v->remapped_pos + p[__aln_end(*v)]->len - u->remapped_pos;
    if (u->remapped_pos != (uint64_t)-1 && v->remapped_pos > u->remapped_pos
        && l >= pint->max_len && ((ii->high && l <= ii->high_bayesian)
        || (ii->high == 0 && l <= opt->max_isize)))
    {
        uint64_t s = __aln(*v, aln).aln.score + __aln(*u, aln).aln.score;
        s *= 10;
        /* penalize for deviation from the average insert size */
        if (ii->high)
            s += (int)(-4.343 * log(.5 * erfc(M_SQRT1_2 * fabs(l - ii->avg) / ii->std)) + .499);

        s = s<<32 | (uint32_t)hash_64(u->remapped_pos<<32 | v->remapped_pos);

        if (s>>32 == pint->o_score>>32) {
            pint->o_n += n_optimal;
        } else if (s>>32 < pint->o_score>>32) {
            pint->subo_n += pint->o_n;
            pint->o_n = n_optimal;
        } else {
            pint->subo_n += n_optimal;
        }

        if (s < pint->o_score) {
            pint->subo_score = pint->o_score;
            pint->o_score = s;
            pint->o_pos[__aln_end(*u)] = *u;
            pint->o_pos[__aln_end(*v)] = *v;
        } else if (s < pint->subo_score) {
            pint->subo_score = s;
        }
    }
}

static void pairing_aux2(const pairing_param_t *param, pairing_internals_t *pint, bwa_seq_t *read, const position_t *pos) {
    const alngrp_t *aln = param->aln;

    const bwt_aln1_t *r = &__aln(*pos, aln).aln;
    read->extra_flag |= SAM_FPP;
    if (read->pos != pos->pos || read->strand != r->a) {
        read->n_mm = r->n_mm; read->n_gapo = r->n_gapo; read->n_gape = r->n_gape; read->strand = r->a;
        read->score = r->score;
        read->pos = pos->pos;
        read->remapped_pos = pos->remapped_pos;
        if (read->mapQ > 0) ++pint->cnt_chg;
    }
}

int find_optimal_pair(const pairing_param_t *param) {
    bwa_seq_t **p = param->p;
    const pos_arr_t *arr = param->arr;
    const alngrp_t *aln = param->aln;
    const pe_opt_t *opt = param->opt;
    int s_mm = param->s_mm;

    int i, j;
    pairing_internals_t pint = { 0 };
    pint.max_len = p[0]->full_len;
    if (pint.max_len < p[1]->full_len)
        pint.max_len = p[1]->full_len;

    pint.o_score = pint.subo_score = (uint64_t)-1;
    pint.o_n = pint.subo_n = 0;
    ks_introsort(position, arr->n, arr->a);
    for (j = 0; j < 2; ++j) {
        pint.last_pos[j][0].pos = pint.last_pos[j][1].pos = (uint64_t)-1;
        pint.last_pos[j][0].remapped_pos = pint.last_pos[j][1].remapped_pos = (uint64_t)-1;
    }
    if (opt->type == BWA_PET_STD) {
        i = 0;
        while (i < arr->n) {
            position_t *pos = &arr->a[i];
            int strand = aln[pos->idx_and_end&1].a[pos->idx_and_end>>1].aln.a;
            int n_optimal = 1;

            if (i < arr->n-1) {
                int k = i;
                while (mappings_overlap(pos, &arr->a[k+1], aln))
                    k++;
                if (k > i) {
                    pos = select_mapping(aln, arr, i, k, &n_optimal);
                    i = k;
                } else ++i;
            } else {
                ++i;
            }

            if (strand == 1) { // reverse strand, then check
                int y = 1 - (pos->idx_and_end&1);
                pairing_aux(param, &pint, &pint.last_pos[y][1], pos, n_optimal);
                pairing_aux(param, &pint, &pint.last_pos[y][0], pos, n_optimal);
            } else { // forward strand, then push
                pint.last_pos[pos->idx_and_end&1][0] = pint.last_pos[pos->idx_and_end&1][1];
                pint.last_pos[pos->idx_and_end&1][1] = *pos;
            }
        }
/*
    } else if (opt->type == BWA_PET_SOLID) {
        for (i = 0; i < arr->n; ++i) {
            uint64_t x = arr->a[i];
            int strand = aln[x&1].a[(uint32_t)x>>1].aln.a;
            if ((strand^x)&1) { // push
                int y = 1 - (x&1);
                __pairing_aux(pint.last_pos[y][1], x);
                __pairing_aux(pint.last_pos[y][0], x);
            } else { // check
                pint.last_pos[x&1][0] = pint.last_pos[x&1][1];
                pint.last_pos[x&1][1] = x;
            }
        }
*/
    } else {
        fprintf(stderr, "[paring] not implemented yet!\n");
        exit(1);
    }
    // set pairing
    if (pint.o_score != (uint64_t)-1) {
        int mapQ_p = 0; // this is the maximum mapping quality when one end is moved
        int rr[2];
        if (pint.o_n == 1) {
            if (pint.subo_score == (uint64_t)-1) {
                mapQ_p = 29; // no sub-optimal pair
            } else if ((pint.subo_score>>32) - (pint.o_score>>32) > s_mm * 10) {
                 mapQ_p = 23; // poor sub-optimal pair
            } else {
                int n = pint.subo_n > 255? 255 : pint.subo_n;
                mapQ_p = ((pint.subo_score>>32) - (pint.o_score>>32)) / 2 - g_log_n[n];
                if (mapQ_p < 0) mapQ_p = 0;
            }
        }
        rr[0] = __aln(pint.o_pos[0], aln).aln.a;
        rr[1] = __aln(pint.o_pos[1], aln).aln.a;
        if ((p[0]->pos == pint.o_pos[0].pos && p[0]->strand == rr[0]) &&
            (p[1]->pos == pint.o_pos[1].pos && p[1]->strand == rr[1]))
        { // both ends not moved
            if (p[0]->mapQ > 0 && p[1]->mapQ > 0) {
                int mapQ = p[0]->mapQ + p[1]->mapQ;
                if (mapQ > 60) mapQ = 60;
                p[0]->mapQ = p[1]->mapQ = mapQ;
            } else {
                if (p[0]->mapQ == 0)
                    p[0]->mapQ = (mapQ_p + 7 < p[1]->mapQ)? mapQ_p + 7 : p[1]->mapQ;
                if (p[1]->mapQ == 0)
                    p[1]->mapQ = (mapQ_p + 7 < p[0]->mapQ)? mapQ_p + 7 : p[0]->mapQ;
            }
        } else if (p[0]->pos == pint.o_pos[0].pos && p[0]->strand == rr[0]) { // [1] moved
            p[1]->seQ = 0; p[1]->mapQ = p[0]->mapQ;
            if (p[1]->mapQ > mapQ_p) p[1]->mapQ = mapQ_p;
        } else if (p[1]->pos == pint.o_pos[1].pos && p[1]->strand == rr[1]) { // [0] moved
            p[0]->seQ = 0; p[0]->mapQ = p[1]->mapQ;
            if (p[0]->mapQ > mapQ_p) p[0]->mapQ = mapQ_p;
        } else { // both ends moved
            p[0]->seQ = p[1]->seQ = 0;
            mapQ_p -= 20;
            if (mapQ_p < 0) mapQ_p = 0;
            p[0]->mapQ = p[1]->mapQ = mapQ_p;
        }
        pairing_aux2(param, &pint, p[0], &pint.o_pos[0]);
        pairing_aux2(param, &pint, p[1], &pint.o_pos[1]);
    }
    return pint.cnt_chg;
}
