#include "bwasw.h"
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <string.h>

#include "threadblock.h"

#define SW_MIN_MATCH_LEN 20
#define SW_MIN_MAPQ 17

typedef struct {
    uint64_t n_tot[2];
    uint64_t n_mapped[2];
} bwa_paired_sw_out_t;

typedef struct {
    dbset_t *dbs;
    int n_seqs;
    bwa_seq_t *seqs[2];
    const pe_opt_t *popt;
    const isize_info_t *ii;
    bwa_paired_sw_out_t *out;
} bwa_paired_sw_data_t;


// cnt = n_mm<<16 | n_gapo<<8 | n_gape
bwa_cigar_t *bwa_sw_core(const dbset_t *dbs, int len, const ubyte_t *seq, int64_t *beg, int reglen,
                      int *n_cigar, uint32_t *_cnt)
{
    bwa_cigar_t *cigar = 0;
    ubyte_t *ref_seq;
    bwtint_t k, x, y, l;
    int path_len, ret;
    AlnParam ap = aln_param_bwa;
    path_t *path, *p;

    // check whether there are too many N's
    if (reglen < SW_MIN_MATCH_LEN || (int64_t)dbs->l_pac - *beg < len) return 0;
    for (k = 0, x = 0; k < len; ++k)
        if (seq[k] >= 4) ++x;
    if ((float)x/len >= 0.25 || len - x < SW_MIN_MATCH_LEN) return 0;

    // get reference subsequence
    ref_seq = (ubyte_t*)calloc(reglen, 1);
    l = dbset_extract_sequence(dbs, dbs->bns, ref_seq, *beg, reglen);
    path = (path_t*)calloc(l+len, sizeof(path_t));

    // do alignment
    ret = aln_local_core(ref_seq, l, (ubyte_t*)seq, len, &ap, path, &path_len, 1, 0);
    if (ret < 0) {
        free(path); free(cigar); free(ref_seq); *n_cigar = 0;
        return 0;
    }
    cigar = bwa_aln_path2cigar(path, path_len, n_cigar);

    // check whether the alignment is good enough
    for (k = 0, x = y = 0; k < *n_cigar; ++k) {
        bwa_cigar_t c = cigar[k];
        if (__cigar_op(c) == FROM_M) x += __cigar_len(c), y += __cigar_len(c);
        else if (__cigar_op(c) == FROM_D) x += __cigar_len(c);
        else y += __cigar_len(c);
    }
    if (x < SW_MIN_MATCH_LEN || y < SW_MIN_MATCH_LEN) { // not good enough
        free(path); free(cigar); free(ref_seq);
        *n_cigar = 0;
        return 0;
    }

    { // update cigar and coordinate;
        int start, end;
        p = path + path_len - 1;
        *beg += (p->i? p->i : 1) - 1;
        start = (p->j? p->j : 1) - 1;
        end = path->j;
        cigar = (bwa_cigar_t*)realloc(cigar, sizeof(bwa_cigar_t) * (*n_cigar + 2));
        if (start) {
            memmove(cigar + 1, cigar, sizeof(bwa_cigar_t) * (*n_cigar));
            cigar[0] = __cigar_create(3, start);
            ++(*n_cigar);
        }
        if (end < len) {
            /*cigar[*n_cigar] = 3<<14 | (len - end);*/
            cigar[*n_cigar] = __cigar_create(3, (len - end));
            ++(*n_cigar);
        }
    }

    { // set *cnt
        int n_mm, n_gapo, n_gape;
        n_mm = n_gapo = n_gape = 0;
        p = path + path_len - 1;
        x = p->i? p->i - 1 : 0; y = p->j? p->j - 1 : 0;
        for (k = 0; k < *n_cigar; ++k) {
            bwa_cigar_t c = cigar[k];
            if (__cigar_op(c) == FROM_M) {
                for (l = 0; l < (__cigar_len(c)); ++l)
                    if (ref_seq[x+l] < 4 && seq[y+l] < 4 && ref_seq[x+l] != seq[y+l]) ++n_mm;
                x += __cigar_len(c), y += __cigar_len(c);
            } else if (__cigar_op(c) == FROM_D) {
                x += __cigar_len(c), ++n_gapo, n_gape += (__cigar_len(c)) - 1;
            } else if (__cigar_op(c) == FROM_I) {
                y += __cigar_len(c), ++n_gapo, n_gape += (__cigar_len(c)) - 1;
            }
        }
        *_cnt = (uint32_t)n_mm<<16 | n_gapo<<8 | n_gape;
    }
    
    free(ref_seq); free(path);
    return cigar;
}

static void set_right_coordinate(
    int64_t *beg, int64_t *end,
    bwa_seq_t *ref, bwa_seq_t *mate,
    const isize_info_t *ii,
    uint64_t l_pac
    )
{
    *beg = (int64_t)ref->remapped_pos + ii->avg - 3*ii->std - mate->len*1.5;
    *end = *beg + 6*ii->std + 2*mate->len;

    if (*beg < (int64_t)ref->remapped_pos + ref->len)
        *beg = ref->remapped_pos + ref->len;
    if (*end > l_pac)
        *end = l_pac;
}

static void set_left_coordinate(
    int64_t *beg, int64_t *end,
    bwa_seq_t *ref, bwa_seq_t *mate,
    const isize_info_t *ii
    )
{
    *beg = (int64_t)ref->remapped_pos + ref->len - ii->avg - 3*ii->std - mate->len*0.5;
    *end = *beg + 6*ii->std + 2*mate->len;

    if (*beg < 0)
        *beg = 0;
    if (*end > ref->remapped_pos)
        *end = ref->remapped_pos;
}

static void bwa_paired_sw_thread(uint32_t idx, uint32_t size, void* data)
{
    bwa_paired_sw_data_t *d = (bwa_paired_sw_data_t*)data;
    dbset_t *dbs = d->dbs;
    int n_seqs = d->n_seqs;
    bwa_seq_t *seqs[2] = { d->seqs[0], d->seqs[1] };
    const pe_opt_t *popt = d->popt;;
    const isize_info_t *ii = d->ii;;
    bwa_paired_sw_out_t *out = &d->out[idx];
    int i;

    // perform mate alignment
    out->n_tot[0] = out->n_tot[1] = out->n_mapped[0] = out->n_mapped[1] = 0;
    for (i = idx; i < n_seqs; i += size) {
        bwa_seq_t *p[2];
        p[0] = seqs[0] + i; p[1] = seqs[1] + i;
        if ((p[0]->mapQ >= SW_MIN_MAPQ || p[1]->mapQ >= SW_MIN_MAPQ) && (p[0]->extra_flag&SAM_FPP) == 0) { // unpaired and one read has high mapQ
            int k, n_cigar[2], is_singleton, mapQ = 0, mq_adjust[2];
            int64_t beg[2], end[2];
            bwa_cigar_t *cigar[2];
            uint32_t cnt[2];

            /* In the following, _pref points to the reference read
             * which must be aligned; _pmate points to its mate which is
             * considered to be modified. */

#define __set_fixed(_pref, _pmate, _beg, _cnt) do {                        \
                _pmate->type = BWA_TYPE_MATESW;                            \
                _pmate->pos = _beg;                                        \
                _pmate->remapped_pos = _beg;                                        \
                _pmate->dbidx = 0;                                        \
                _pmate->remapped_dbidx = 0;                                        \
                _pmate->seQ = _pref->seQ;                                \
                _pmate->strand = (popt->type == BWA_PET_STD)? 1 - _pref->strand : _pref->strand; \
                _pmate->n_mm = _cnt>>16; _pmate->n_gapo = _cnt>>8&0xff; _pmate->n_gape = _cnt&0xff; \
                _pmate->extra_flag |= SAM_FPP;                            \
                _pref->extra_flag |= SAM_FPP;                            \
            } while (0)

            mq_adjust[0] = mq_adjust[1] = 255; // not effective
            is_singleton = (p[0]->type == BWA_TYPE_NO_MATCH || p[1]->type == BWA_TYPE_NO_MATCH)? 1 : 0;

            ++out->n_tot[is_singleton];
            cigar[0] = cigar[1] = 0;
            n_cigar[0] = n_cigar[1] = 0;
            if (popt->type != BWA_PET_STD && popt->type != BWA_PET_SOLID)
                continue; // other types of pairing is not considered
            for (k = 0; k < 2; ++k) { // p[1-k] is the reference read and p[k] is the read considered to be modified
                ubyte_t *seq;
                if (p[1-k]->type == BWA_TYPE_NO_MATCH) continue; // if p[1-k] is unmapped, skip
                if (popt->type == BWA_PET_STD) {
                    if (p[1-k]->strand == 0) { // then the mate is on the reverse strand and has larger coordinate
                        set_right_coordinate(beg+k, end+k, p[1-k], p[k], ii, dbs->l_pac);
                        seq = p[k]->rseq;
                    } else { // then the mate is on forward stand and has smaller coordinate
                        set_left_coordinate(beg+k, end+k, p[1-k], p[k], ii);
                        seq = p[k]->seq;
                        seq_reverse(p[k]->len, seq, 0); // because ->seq is reversed; this will reversed back shortly
                    }
                } else { // BWA_PET_SOLID
                    if (p[1-k]->strand == 0) { // R3-F3 pairing
                        if (k == 0) 
                            set_left_coordinate(beg+k, end+k, p[1-k], p[k], ii); // p[k] is R3
                        else
                            set_right_coordinate(beg+k, end+k, p[1-k], p[k], ii, dbs->l_pac); // p[k] is F3
                        seq = p[k]->rseq;
                        seq_reverse(p[k]->len, seq, 0); // because ->seq is reversed
                    } else { // F3-R3 pairing
                        if (k == 0)
                            set_right_coordinate(beg+k, end+k, p[1-k], p[k], ii, dbs->l_pac); // p[k] is R3
                        else
                            set_left_coordinate(beg+k, end+k, p[1-k], p[k], ii); // p[k] is F3
                        seq = p[k]->seq;
                    }
                }
                // perform SW alignment
                cigar[k] = bwa_sw_core(dbs, p[k]->len, seq, &beg[k], end[k] - beg[k], &n_cigar[k], &cnt[k]);
                if (cigar[k] && p[k]->type != BWA_TYPE_NO_MATCH) { // re-evaluate cigar[k]
                    int s_old, clip = 0, s_new;
                    if (__cigar_op(cigar[k][0]) == 3) clip += __cigar_len(cigar[k][0]);
                    if (__cigar_op(cigar[k][n_cigar[k]-1]) == 3) clip += __cigar_len(cigar[k][n_cigar[k]-1]);
                    s_old = (int)((p[k]->n_mm * 9 + p[k]->n_gapo * 13 + p[k]->n_gape * 2) / 3. * 8. + .499);
                    s_new = (int)(((cnt[k]>>16) * 9 + (cnt[k]>>8&0xff) * 13 + (cnt[k]&0xff) * 2 + clip * 3) / 3. * 8. + .499);
                    s_old += -4.343 * log(ii->ap_prior / dbs->l_pac);
                    s_new += (int)(-4.343 * log(.5 * erfc(M_SQRT1_2 * 1.5) + .499)); // assume the mapped isize is 1.5\sigma
                    if (s_old < s_new) { // reject SW alignment
                        mq_adjust[k] = s_new - s_old;
                        free(cigar[k]); cigar[k] = 0; n_cigar[k] = 0;
                    } else mq_adjust[k] = s_old - s_new;
                }
                // now revserse sequence back such that p[*]->seq looks untouched
                if (popt->type == BWA_PET_STD) {
                    if (p[1-k]->strand == 1) seq_reverse(p[k]->len, seq, 0);
                } else {
                    if (p[1-k]->strand == 0) seq_reverse(p[k]->len, seq, 0);
                }
            }
            k = -1; // no read to be changed
            if (cigar[0] && cigar[1]) {
                k = p[0]->mapQ < p[1]->mapQ? 0 : 1; // p[k] to be fixed
                mapQ = abs(p[1]->mapQ - p[0]->mapQ);
            } else if (cigar[0]) k = 0, mapQ = p[1]->mapQ;
            else if (cigar[1]) k = 1, mapQ = p[0]->mapQ;
            if (k >= 0 && p[k]->pos != beg[k]) {
                ++out->n_mapped[is_singleton];
                { // recalculate mapping quality
                    int tmp = (int)p[1-k]->mapQ - p[k]->mapQ/2 - 8;
                    if (tmp <= 0) tmp = 1;
                    if (mapQ > tmp) mapQ = tmp;
                    p[k]->mapQ = p[1-k]->mapQ = mapQ;
                    p[k]->seQ = p[1-k]->seQ = p[1-k]->seQ < mapQ? p[1-k]->seQ : mapQ;
                    if (p[k]->mapQ > mq_adjust[k]) p[k]->mapQ = mq_adjust[k];
                    if (p[k]->seQ > mq_adjust[k]) p[k]->seQ = mq_adjust[k];
                }
                // update CIGAR
                free(p[k]->cigar); p[k]->cigar = cigar[k]; cigar[k] = 0;
                p[k]->n_cigar = n_cigar[k];
                // update the rest of information
                __set_fixed(p[1-k], p[k], beg[k], cnt[k]);
            }
            free(cigar[0]); free(cigar[1]);
        }
    }
}

void bwa_paired_sw(dbset_t *dbs, int n_seqs, bwa_seq_t *seqs[2], const pe_opt_t *popt, const isize_info_t *ii)
{
    int i;
    uint64_t n_tot[2] = {0,0};
    uint64_t n_mapped[2] = {0,0};
    bwa_paired_sw_data_t td;

    dbset_load_pac(dbs);

    if (popt->is_sw && ii->avg >= 0.0) {
        td.dbs = dbs;
        td.n_seqs = n_seqs;
        td.seqs[0] = seqs[0];
        td.seqs[1] = seqs[1];
        td.popt = popt;
        td.ii = ii;
        td.out = calloc(popt->n_threads, sizeof(bwa_paired_sw_out_t));
        threadblock_exec(popt->n_threads, &bwa_paired_sw_thread, &td);

        for (i = 0; i < popt->n_threads; ++i) {
            n_tot[0] += td.out[i].n_tot[0];
            n_tot[1] += td.out[i].n_tot[1];
            n_mapped[0] += td.out[i].n_mapped[0];
            n_mapped[1] += td.out[i].n_mapped[1];
        }
        free(td.out);

        fprintf(stderr, "[bwa_paired_sw] %lld out of %lld Q%d singletons are mated.\n",
                (long long)n_mapped[1], (long long)n_tot[1], SW_MIN_MAPQ);
        fprintf(stderr, "[bwa_paired_sw] %lld out of %lld Q%d discordant pairs are fixed.\n",
                (long long)n_mapped[0], (long long)n_tot[0], SW_MIN_MAPQ);
    }

    dbset_unload_pac(dbs);
}


