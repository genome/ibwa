#include "saiset.h"

#include "ksort.h"
#include "utils.h"

#include <limits.h>

static int __cmp_alignment_scores(const alignment_t a, const alignment_t b) {
    return a.aln.score < b.aln.score;
}
KSORT_INIT(alignment, alignment_t, __cmp_alignment_scores);



saiset_t *saiset_create(int n, const char **files[]) {
    int i, j;
    saiset_t *s = calloc(1, sizeof(saiset_t));
    s->count = n;
    s->fp[0] = calloc(n, sizeof(FILE*));
    s->fp[1] = calloc(n, sizeof(FILE*));

    for (i = 0; i < 2; ++i) {
        for (j = 0; j < n; ++j) {
            if (!(s->fp[i][j] = xopen(files[j][i], "r"))) {
                fprintf(stderr, "[saiset_create] failed to open sai file %s\n", files[j][i]);
            }
            /* TODO: verify opt records match! */
            fread(&s->opt[i], sizeof(gap_opt_t), 1, s->fp[i][j]);
        }
    }
    return s;
}

void saiset_destroy(saiset_t *s) {
    int i, j;
    for (i = 0; i < 2; ++i) {
        for (j = 0; j < s->count; ++j) {
            fclose(s->fp[i][j]);
        }
        free(s->fp[i]);
    }
    free(s);
}

alngrp_t *alngrp_create(const dbset_t *dbs, const saiset_t *s, int which) {
    int i, j;
    uint32_t count = 0;
    alngrp_t *ag = (alngrp_t*)calloc(1, sizeof(alngrp_t));
    int best_score;

    for (i = 0; i < s->count; ++i) {
        fread(&count, 4, 1, s->fp[which][i]);
        kv_resize(alignment_t, *ag, ag->n+count);
        for (j = 0; j < count; ++j) {
            fread(&ag->a[ag->n].aln, sizeof(bwt_aln1_t), 1, s->fp[which][i]);
            ag->a[ag->n].db = dbs->db[i];
            ag->a[ag->n].dbidx = i;
            ag->a[ag->n].remapped = 0;
            ++ag->n;
        }
    }

    /* only need to do this if we have multiple sai streams */
    if (s->count > 1 && ag->n > 0) {
        ks_introsort(alignment, ag->n, ag->a);
        best_score = ag->a[0].aln.score;
        
        /* filter out bad scores */
        for (i = 0; i < ag->n; ++i) {
            if (ag->a[i].aln.score > best_score + s->opt[0].s_mm) {
                ag->n = i;
                break;
            }
        }
    }

    return ag;
}

void alngrp_destroy(alngrp_t *ag) {
    free(ag->a);
    free(ag);
}

void select_sai(const alngrp_t *ag, bwa_seq_t *s, int *main_idx) {
    int i, cnt, best;
    if (ag->n == 0) {
        s->type = BWA_TYPE_NO_MATCH;
        s->c1 = s->c2 = 0;
        return;
    }

    if (main_idx) {
        best = ag->a[0].aln.score;
        for (i = cnt = 0; i < ag->n; ++i) {
            const bwt_aln1_t *p = &ag->a[i].aln;
            if (p->score > best) break;
            if (drand48() * (p->l - p->k + 1 + cnt) > (double)cnt) {
                *main_idx = i;
                s->n_mm = p->n_mm; s->n_gapo = p->n_gapo; s->n_gape = p->n_gape; s->strand = p->a;
                s->score = p->score;
                s->sa = p->k + (bwtint_t)((p->l - p->k + 1) * drand48());
            }
            cnt += p->l - p->k + 1;
        }
        s->c1 = cnt;
        for (; i < ag->n; ++i) cnt += ag->a[i].aln.l - ag->a[i].aln.k + 1;
        s->c2 = cnt - s->c1;
        s->type = s->c1 > 1? BWA_TYPE_REPEAT : BWA_TYPE_UNIQUE;
    }
}

void select_sai_multi(const alngrp_t *ag, bwa_seq_t *s, int n_multi) {
    int k, rest, n_occ, z = 0;
    for (k = n_occ = 0; k < ag->n; ++k) {
        const bwt_aln1_t *q = &ag->a[k].aln ;
        n_occ += q->l - q->k + 1;
    }
    if (s->multi) free(s->multi); 
    if (n_occ > n_multi + 1) { // if there are too many hits, generate none of them
        s->multi = 0; s->n_multi = 0;
        return;
    }
    /* The following code is more flexible than what is required
     * here. In principle, due to the requirement above, we can
     * simply output all hits, but the following samples "rest"
     * number of random hits. */
    rest = n_occ > n_multi + 1? n_multi + 1 : n_occ; // find one additional for ->sa
    s->multi = calloc(rest, sizeof(bwt_multi1_t));
    for (k = 0; k < ag->n; ++k) {
        const bwtdb_t *db = ag->a[k].db;
        const bwt_aln1_t *q = &ag->a[k].aln;
        if (q->l - q->k + 1 <= rest) {
            bwtint_t l;
            for (l = q->k; l <= q->l; ++l) {
                s->multi[z].pos = bwtdb_sa2seq(db, q->a, l, s->len);
                s->multi[z].gap = q->n_gapo + q->n_gape;
                s->multi[z].mm = q->n_mm;
                s->multi[z++].strand = q->a;
            }
            rest -= q->l - q->k + 1;
        } else { // Random sampling (http://code.activestate.com/recipes/272884/). In fact, we never come here. 
            int j, i, k;
            for (j = rest, i = q->l - q->k + 1, k = 0; j > 0; --j) {
                double p = 1.0, x = drand48();
                while (x < p) p -= p * j / (i--);
                s->multi[z].pos = bwtdb_sa2seq(db, q->a, q->l-1, s->len);
                s->multi[z].gap = q->n_gapo + q->n_gape;
                s->multi[z].mm = q->n_mm;
                s->multi[z++].strand = q->a;
            }
            rest = 0;
            break;
        }
    }
    s->n_multi = z;
    for (k = z = 0; k < s->n_multi; ++k)
        if (s->multi[k].pos != s->pos)
            s->multi[z++] = s->multi[k];
    s->n_multi = z < n_multi? z : n_multi;
}


