#include <unistd.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <string.h>
#include "bwtaln.h"
#include "kvec.h"
#include "bntseq.h"
#include "utils.h"
#include "stdaln.h"
#include "bwtcache.h"
#include "dbset.h"
#include "saiset.h"

#include "threadblock.h"

typedef struct {
	int count;
	const char *fq[2];
	kvec_t(const char*) prefixes;
	kvec_t(const char**) sai_pair;
} pe_inputs_t;

typedef struct {
	double avg, std, ap_prior;
	bwtint_t low, high, high_bayesian;
} isize_info_t;

typedef kvec_t(uint64_t) pos_arr_t;
typedef kvec_t(bwt_aln1_t) aln_buf_t;

typedef struct {
	dbset_t *dbs;
	const alngrp_t *buf[2];
	int n_seqs;
	bwa_seq_t *seqs[2];
	isize_info_t *ii;
	const pe_opt_t *opt;
	const gap_opt_t *gopt;
	int *cnt_chg;
} cal_pac_pos_params_t;

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

#include "ksort.h"
KSORT_INIT_GENERIC(uint64_t)

#define MIN_HASH_WIDTH 1000

extern int g_log_n[256]; // in bwase.c

void bwase_initialize();
void bwa_aln2seq_core(int n_aln, const bwt_aln1_t *aln, bwa_seq_t *s, int set_main, int n_multi);
void bwa_aln2seq(int n_aln, const bwt_aln1_t *aln, bwa_seq_t *s);
int bwa_approx_mapQ(const bwa_seq_t *p, int mm);
void bwa_print_sam1(const dbset_t *dbs, bwa_seq_t *p, const bwa_seq_t *mate, int mode, int max_top2);
void bwa_refine_gapped(dbset_t *dbs, int n_seqs, bwa_seq_t *seqs);
bntseq_t *bwa_open_nt(const char *prefix);
void bwa_print_sam_PG();

pe_opt_t *bwa_init_pe_opt()
{
	pe_opt_t *po;
	po = (pe_opt_t*)calloc(1, sizeof(pe_opt_t));
	po->max_isize = 500;
	po->force_isize = 0;
	po->max_occ = 100000;
	po->n_multi = 3;
	po->N_multi = 10;
	po->type = BWA_PET_STD;
	po->is_sw = 1;
	po->ap_prior = 1e-5;
	po->n_threads = 1;
	return po;
}

static inline uint64_t hash_64(uint64_t key)
{
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
/*
static double ierfc(double x) // inverse erfc(); iphi(x) = M_SQRT2 *ierfc(2 * x);
{
	const double a = 0.140012;
	double b, c;
	b = log(x * (2 - x));
	c = 2./M_PI/a + b / 2.;
	return sqrt(sqrt(c * c - b / a) - c);
}
*/

// for normal distribution, this is about 3std
#define OUTLIER_BOUND 2.0

static int infer_isize(int n_seqs, bwa_seq_t *seqs[2], isize_info_t *ii, double ap_prior, int64_t L)
{
	uint64_t x, *isizes, n_ap = 0;
	int n, i, tot, p25, p75, p50, max_len = 1, tmp;
	double skewness = 0.0, kurtosis = 0.0, y;

	ii->avg = ii->std = -1.0;
	ii->low = ii->high = ii->high_bayesian = 0;
	isizes = (uint64_t*)calloc(n_seqs, 8);
	for (i = 0, tot = 0; i != n_seqs; ++i) {
		bwa_seq_t *p[2];
		p[0] = seqs[0] + i; p[1] = seqs[1] + i;
		if (p[0]->mapQ >= 20 && p[1]->mapQ >= 20) {
			x = (p[0]->pos < p[1]->pos)? p[1]->pos + p[1]->len - p[0]->pos : p[0]->pos + p[0]->len - p[1]->pos;
			if (x < 100000) isizes[tot++] = x;
		}
		if (p[0]->len > max_len) max_len = p[0]->len;
		if (p[1]->len > max_len) max_len = p[1]->len;
	}
	if (tot < 20) {
		fprintf(stderr, "[infer_isize] fail to infer insert size: too few good pairs\n");
		free(isizes);
		return -1;
	}
	ks_introsort(uint64_t, tot, isizes);
	p25 = isizes[(int)(tot*0.25 + 0.5)];
	p50 = isizes[(int)(tot*0.50 + 0.5)];
	p75 = isizes[(int)(tot*0.75 + 0.5)];
	tmp  = (int)(p25 - OUTLIER_BOUND * (p75 - p25) + .499);
	ii->low = tmp > max_len? tmp : max_len; // ii->low is unsigned
	ii->high = (int)(p75 + OUTLIER_BOUND * (p75 - p25) + .499);
	for (i = 0, x = n = 0; i < tot; ++i)
		if (isizes[i] >= ii->low && isizes[i] <= ii->high)
			++n, x += isizes[i];
	ii->avg = (double)x / n;
	for (i = 0; i < tot; ++i) {
		if (isizes[i] >= ii->low && isizes[i] <= ii->high) {
			double tmp = (isizes[i] - ii->avg) * (isizes[i] - ii->avg);
			ii->std += tmp;
			skewness += tmp * (isizes[i] - ii->avg);
			kurtosis += tmp * tmp;
		}
	}
	kurtosis = kurtosis/n / (ii->std / n * ii->std / n) - 3;
	ii->std = sqrt(ii->std / n); // it would be better as n-1, but n is usually very large
	skewness = skewness / n / (ii->std * ii->std * ii->std);
	for (y = 1.0; y < 10.0; y += 0.01)
		if (.5 * erfc(y / M_SQRT2) < ap_prior / L * (y * ii->std + ii->avg)) break;
	ii->high_bayesian = (bwtint_t)(y * ii->std + ii->avg + .499);
	for (i = 0; i < tot; ++i)
		if (isizes[i] > ii->high_bayesian) ++n_ap;
	ii->ap_prior = .01 * (n_ap + .01) / tot;
	if (ii->ap_prior < ap_prior) ii->ap_prior = ap_prior;
	free(isizes);
	fprintf(stderr, "[infer_isize] (25, 50, 75) percentile: (%d, %d, %d)\n", p25, p50, p75);
	if (isnan(ii->std) || p75 > 100000) {
		ii->low = ii->high = ii->high_bayesian = 0; ii->avg = ii->std = -1.0;
		fprintf(stderr, "[infer_isize] fail to infer insert size: weird pairing\n");
		return -1;
	}
	for (y = 1.0; y < 10.0; y += 0.01)
		if (.5 * erfc(y / M_SQRT2) < ap_prior / L * (y * ii->std + ii->avg)) break;
	ii->high_bayesian = (bwtint_t)(y * ii->std + ii->avg + .499);
	fprintf(stderr, "[infer_isize] low and high boundaries: %d and %d for estimating avg and std\n", ii->low, ii->high);
	fprintf(stderr, "[infer_isize] inferred external isize from %d pairs: %.3lf +/- %.3lf\n", n, ii->avg, ii->std);
	fprintf(stderr, "[infer_isize] skewness: %.3lf; kurtosis: %.3lf; ap_prior: %.2e\n", skewness, kurtosis, ii->ap_prior);
	fprintf(stderr, "[infer_isize] inferred maximum insert size: %d (%.2lf sigma)\n", ii->high_bayesian, y);
	return 0;
}

static int pairing(bwa_seq_t *p[2], const pos_arr_t *arr, const alngrp_t aln[2], const pe_opt_t *opt, int s_mm, const isize_info_t *ii)
{
	int i, j, o_n, subo_n, cnt_chg = 0, low_bound = ii->low, max_len;
	uint64_t last_pos[2][2], o_pos[2], subo_score, o_score;
	max_len = p[0]->full_len;
	if (max_len < p[1]->full_len) max_len = p[1]->full_len;
	if (low_bound < max_len) low_bound = max_len;

	// here v>=u. When ii is set, we check insert size with ii; otherwise with opt->max_isize
#define __pairing_aux(u,v) do {											\
		bwtint_t l = ((v)>>32) + p[(v)&1]->len - ((u)>>32);				\
		if ((u) != (uint64_t)-1 && (v)>>32 > (u)>>32 && l >= max_len	\
			&& ((ii->high && l <= ii->high_bayesian) || (ii->high == 0 && l <= opt->max_isize))) \
		{																\
			uint64_t s = aln[(v)&1].a[(uint32_t)(v)>>1].aln.score + aln[(u)&1].a[(uint32_t)(u)>>1].aln.score; \
			s *= 10;													\
			if (ii->high) s += (int)(-4.343 * log(.5 * erfc(M_SQRT1_2 * fabs(l - ii->avg) / ii->std)) + .499); \
			s = s<<32 | (uint32_t)hash_64((u)>>32<<32 | (v)>>32);		\
			if (s>>32 == o_score>>32) ++o_n;							\
			else if (s>>32 < o_score>>32) { subo_n += o_n; o_n = 1; }	\
			else ++subo_n;												\
			if (s < o_score) subo_score = o_score, o_score = s, o_pos[(u)&1] = (u), o_pos[(v)&1] = (v); \
			else if (s < subo_score) subo_score = s;					\
		}																\
	} while (0)

#define __pairing_aux2(q, w) do {										\
		const bwt_aln1_t *r = &aln[(w)&1].a[((uint32_t)(w)>>1)].aln;	\
		(q)->extra_flag |= SAM_FPP;										\
		if ((q)->pos != (w)>>32 || (q)->strand != r->a) {				\
			(q)->n_mm = r->n_mm; (q)->n_gapo = r->n_gapo; (q)->n_gape = r->n_gape; (q)->strand = r->a; \
			(q)->score = r->score;										\
			(q)->pos = (w)>>32;											\
			if ((q)->mapQ > 0) ++cnt_chg;								\
		}																\
	} while (0)

	o_score = subo_score = (uint64_t)-1;
	o_n = subo_n = 0;
	ks_introsort(uint64_t, arr->n, arr->a);
	for (j = 0; j < 2; ++j) last_pos[j][0] = last_pos[j][1] = (uint64_t)-1;
	if (opt->type == BWA_PET_STD) {
		for (i = 0; i < arr->n; ++i) {
			uint64_t x = arr->a[i];
			int strand = aln[x&1].a[(uint32_t)x>>1].aln.a;
			if (strand == 1) { // reverse strand, then check
				int y = 1 - (x&1);
				__pairing_aux(last_pos[y][1], x);
				__pairing_aux(last_pos[y][0], x);
			} else { // forward strand, then push
				last_pos[x&1][0] = last_pos[x&1][1];
				last_pos[x&1][1] = x;
			}
		}
	} else if (opt->type == BWA_PET_SOLID) {
		for (i = 0; i < arr->n; ++i) {
			uint64_t x = arr->a[i];
			int strand = aln[x&1].a[(uint32_t)x>>1].aln.a;
			if ((strand^x)&1) { // push
				int y = 1 - (x&1);
				__pairing_aux(last_pos[y][1], x);
				__pairing_aux(last_pos[y][0], x);
			} else { // check
				last_pos[x&1][0] = last_pos[x&1][1];
				last_pos[x&1][1] = x;
			}
		}
	} else {
		fprintf(stderr, "[paring] not implemented yet!\n");
		exit(1);
	}
	// set pairing
	//fprintf(stderr, "[%d, %d, %d, %d]\n", arr->n, (int)(o_score>>32), (int)(subo_score>>32), o_n);
	if (o_score != (uint64_t)-1) {
		int mapQ_p = 0; // this is the maximum mapping quality when one end is moved
		int rr[2];
		//fprintf(stderr, "%d, %d\n", o_n, subo_n);
		if (o_n == 1) {
			if (subo_score == (uint64_t)-1) mapQ_p = 29; // no sub-optimal pair
			else if ((subo_score>>32) - (o_score>>32) > s_mm * 10) mapQ_p = 23; // poor sub-optimal pair
			else {
				int n = subo_n > 255? 255 : subo_n;
				mapQ_p = ((subo_score>>32) - (o_score>>32)) / 2 - g_log_n[n];
				if (mapQ_p < 0) mapQ_p = 0;
			}
		}
		rr[0] = aln[o_pos[0]&1].a[(uint32_t)o_pos[0]>>1].aln.a;
		rr[1] = aln[o_pos[1]&1].a[(uint32_t)o_pos[1]>>1].aln.a;
		if ((p[0]->pos == o_pos[0]>>32 && p[0]->strand == rr[0]) && (p[1]->pos == o_pos[1]>>32 && p[1]->strand == rr[1])) { // both ends not moved
			if (p[0]->mapQ > 0 && p[1]->mapQ > 0) {
				int mapQ = p[0]->mapQ + p[1]->mapQ;
				if (mapQ > 60) mapQ = 60;
				p[0]->mapQ = p[1]->mapQ = mapQ;
			} else {
				if (p[0]->mapQ == 0) p[0]->mapQ = (mapQ_p + 7 < p[1]->mapQ)? mapQ_p + 7 : p[1]->mapQ;
				if (p[1]->mapQ == 0) p[1]->mapQ = (mapQ_p + 7 < p[0]->mapQ)? mapQ_p + 7 : p[0]->mapQ;
			}
		} else if (p[0]->pos == o_pos[0]>>32 && p[0]->strand == rr[0]) { // [1] moved
			p[1]->seQ = 0; p[1]->mapQ = p[0]->mapQ;
			if (p[1]->mapQ > mapQ_p) p[1]->mapQ = mapQ_p;
		} else if (p[1]->pos == o_pos[1]>>32 && p[1]->strand == rr[1]) { // [0] moved
			p[0]->seQ = 0; p[0]->mapQ = p[1]->mapQ;
			if (p[0]->mapQ > mapQ_p) p[0]->mapQ = mapQ_p;
		} else { // both ends moved
			p[0]->seQ = p[1]->seQ = 0;
			mapQ_p -= 20;
			if (mapQ_p < 0) mapQ_p = 0;
			p[0]->mapQ = p[1]->mapQ = mapQ_p;
		}
		__pairing_aux2(p[0], o_pos[0]);
		__pairing_aux2(p[1], o_pos[1]);
	}
	return cnt_chg;
}

static void bwa_cal_pac_pos_pe_thread(uint32_t idx, uint32_t size, void *data)
{
	cal_pac_pos_params_t const *tdata = (cal_pac_pos_params_t*)data;
	const alngrp_t *buf[2] = {tdata->buf[0], tdata->buf[1]};
	int n_seqs = tdata->n_seqs;
	bwa_seq_t *seqs[2] = {tdata->seqs[0], tdata->seqs[1]};
	isize_info_t *ii = tdata->ii;
	const pe_opt_t *opt = tdata->opt;
	const gap_opt_t *gopt = tdata->gopt;

	int i,j;
	alngrp_t aln[2] = {{0,0}, {0,0}};
	pos_arr_t arr = {0};

	tdata->cnt_chg[idx] = 0;
	for (i = idx; i < n_seqs; i += size) {
		bwa_seq_t *p[2];
		for (j = 0; j < 2; ++j) {
			p[j] = seqs[j] + i;
			aln[j] = buf[j][i];
		}
		if ((p[0]->type == BWA_TYPE_UNIQUE || p[0]->type == BWA_TYPE_REPEAT)
			&& (p[1]->type == BWA_TYPE_UNIQUE || p[1]->type == BWA_TYPE_REPEAT))
		{ // only when both ends mapped
			uint64_t x;
			int j, k, n_occ[2];
			for (j = 0; j < 2; ++j) {
				n_occ[j] = 0;
				for (k = 0; k < aln[j].n; ++k)
					n_occ[j] += aln[j].a[k].aln.l - aln[j].a[k].aln.k + 1;
			}
			if (n_occ[0] > opt->max_occ || n_occ[1] > opt->max_occ) continue;
			arr.n = 0;
			for (j = 0; j < 2; ++j) {
				for (k = 0; k < aln[j].n; ++k) {
					alignment_t *ar = &aln[j].a[k];
					bwtint_t l;
					if (ar->aln.l - ar->aln.k + 1 >= MIN_HASH_WIDTH) { // then check hash table
						poslist_t pos = bwtdb_cached_sa2seq(ar->db, &ar->aln, p[j]->len);
						for (l = 0; l < pos.n; ++l) {
							x  = pos.a[l];
							x = x<<32 | k<<1 | j;
							kv_push(uint64_t, arr, x);
						}
					} else { // then calculate on the fly
						for (l = ar->aln.k; l <= ar->aln.l; ++l) {
							x = bwtdb_sa2seq(ar->db, ar->aln.a, l, p[j]->len);
							x = x<<32 | k<<1 | j;
							kv_push(uint64_t, arr, x);
						}
					}

				}
			}
			tdata->cnt_chg[idx] += pairing(p, &arr, aln, opt, gopt->s_mm, ii);
		}

		if (opt->N_multi || opt->n_multi) {
			for (j = 0; j < 2; ++j) {
				if (p[j]->type != BWA_TYPE_NO_MATCH) {
					int max_multi = opt->n_multi;
					if (!(p[j]->extra_flag&SAM_FPP) && p[1-j]->type != BWA_TYPE_NO_MATCH)
						max_multi = p[j]->c1+p[j]->c2-1 > opt->N_multi? opt->n_multi : opt->N_multi;
					select_sai_multi(&aln[j], p[j], max_multi);
				}
			}
		}
	}
	kv_destroy(arr);
}

int bwa_cal_pac_pos_pe(dbset_t *dbs, int n_seqs, bwa_seq_t *seqs[2], saiset_t *saiset, isize_info_t *ii,
					   const pe_opt_t *opt, const gap_opt_t *gopt, const isize_info_t *last_ii)
{
	int i, j, cnt_chg = 0;
	alngrp_t *aln[2];
	alngrp_t *buf[2];
	cal_pac_pos_params_t tp;

	buf[0] = (alngrp_t*)calloc(n_seqs, sizeof(alngrp_t));
	buf[1] = (alngrp_t*)calloc(n_seqs, sizeof(alngrp_t));

	// SE
	for (i = 0; i != n_seqs; ++i) {
		bwa_seq_t *p[2];
		for (j = 0; j < 2; ++j) {
			int main_idx = 0;
			p[j] = seqs[j] + i;
			p[j]->n_multi = 0;
			p[j]->extra_flag |= SAM_FPD | (j == 0? SAM_FR1 : SAM_FR2);
			aln[j] = alngrp_create(dbs, saiset, j);
			kv_copy(alignment_t, buf[j][i], *aln[j]); // backup aln[j]

			// generate SE alignment and mapping quality
			select_sai(aln[j], p[j], &main_idx);

			if (p[j]->type == BWA_TYPE_UNIQUE || p[j]->type == BWA_TYPE_REPEAT) {
				alignment_t *main_aln = &aln[j]->a[main_idx];
				int max_diff = gopt->fnr > 0.0? bwa_cal_maxdiff(p[j]->len, BWA_AVG_ERR, gopt->fnr) : gopt->max_diff;
				p[j]->pos = bwtdb_sa2seq(main_aln->db, p[j]->strand, p[j]->sa, p[j]->len);
				p[j]->seQ = p[j]->mapQ = bwa_approx_mapQ(p[j], max_diff);
			}
		}
	}

	// infer isize
	infer_isize(n_seqs, seqs, ii, opt->ap_prior, dbs->total_bwt_seq_len[0]);
	if (ii->avg < 0.0 && last_ii->avg > 0.0) *ii = *last_ii;
	if (opt->force_isize) {
		fprintf(stderr, "[%s] discard insert size estimate as user's request.\n", __func__);
		ii->low = ii->high = 0; ii->avg = ii->std = -1.0;
	}

	// PE
	for (i = 0; i < 2; ++i) {
		tp.dbs = dbs;
		tp.buf[i] = buf[i];
		tp.seqs[i] = seqs[i];
	}
	tp.n_seqs = n_seqs;
	tp.ii = ii;
	tp.opt = opt;
	tp.gopt = gopt;
	tp.cnt_chg = calloc(opt->n_threads, sizeof(int));
	threadblock_exec(opt->n_threads, &bwa_cal_pac_pos_pe_thread, &tp);

	for (i = 0; i < opt->n_threads; ++i)
		cnt_chg += tp.cnt_chg[i];

	// free
	free(tp.cnt_chg);
	alngrp_destroy(aln[0]);
	alngrp_destroy(aln[1]);
	for (i = 0; i < n_seqs; ++i) {
		kv_destroy(buf[0][i]);
		kv_destroy(buf[1][i]);
	}
	free(buf[0]); free(buf[1]);
	return cnt_chg;
}

#define SW_MIN_MATCH_LEN 20
#define SW_MIN_MAPQ 17

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

	dbset_load_pac(dbs);

	// perform mate alignment
	/* TODO: make out->n_mapped|tot thread safe! */
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

#define __set_rght_coor(_a, _b, _pref, _pmate) do {						\
				(_a) = (int64_t)_pref->pos + ii->avg - 3 * ii->std - _pmate->len * 1.5; \
				(_b) = (_a) + 6 * ii->std + 2 * _pmate->len;			\
				if ((_a) < (int64_t)_pref->pos + _pref->len) (_a) = _pref->pos + _pref->len; \
				if ((_b) > dbs->l_pac) (_b) = dbs->l_pac;				\
			} while (0)

#define __set_left_coor(_a, _b, _pref, _pmate) do {						\
				(_a) = (int64_t)_pref->pos + _pref->len - ii->avg - 3 * ii->std - _pmate->len * 0.5; \
				(_b) = (_a) + 6 * ii->std + 2 * _pmate->len;			\
				if ((_a) < 0) (_a) = 0;									\
				if ((_b) > _pref->pos) (_b) = _pref->pos;				\
			} while (0)
			
#define __set_fixed(_pref, _pmate, _beg, _cnt) do {						\
				_pmate->type = BWA_TYPE_MATESW;							\
				_pmate->pos = _beg;										\
				_pmate->seQ = _pref->seQ;								\
				_pmate->strand = (popt->type == BWA_PET_STD)? 1 - _pref->strand : _pref->strand; \
				_pmate->n_mm = _cnt>>16; _pmate->n_gapo = _cnt>>8&0xff; _pmate->n_gape = _cnt&0xff; \
				_pmate->extra_flag |= SAM_FPP;							\
				_pref->extra_flag |= SAM_FPP;							\
			} while (0)

			mq_adjust[0] = mq_adjust[1] = 255; // not effective
			is_singleton = (p[0]->type == BWA_TYPE_NO_MATCH || p[1]->type == BWA_TYPE_NO_MATCH)? 1 : 0;

			++out->n_tot[is_singleton];
			cigar[0] = cigar[1] = 0;
			n_cigar[0] = n_cigar[1] = 0;
			if (popt->type != BWA_PET_STD && popt->type != BWA_PET_SOLID) continue; // other types of pairing is not considered
			for (k = 0; k < 2; ++k) { // p[1-k] is the reference read and p[k] is the read considered to be modified
				ubyte_t *seq;
				if (p[1-k]->type == BWA_TYPE_NO_MATCH) continue; // if p[1-k] is unmapped, skip
				if (popt->type == BWA_PET_STD) {
					if (p[1-k]->strand == 0) { // then the mate is on the reverse strand and has larger coordinate
						__set_rght_coor(beg[k], end[k], p[1-k], p[k]);
						seq = p[k]->rseq;
					} else { // then the mate is on forward stand and has smaller coordinate
						__set_left_coor(beg[k], end[k], p[1-k], p[k]);
						seq = p[k]->seq;
						seq_reverse(p[k]->len, seq, 0); // because ->seq is reversed; this will reversed back shortly
					}
				} else { // BWA_PET_SOLID
					if (p[1-k]->strand == 0) { // R3-F3 pairing
						if (k == 0) __set_left_coor(beg[k], end[k], p[1-k], p[k]); // p[k] is R3
						else __set_rght_coor(beg[k], end[k], p[1-k], p[k]); // p[k] is F3
						seq = p[k]->rseq;
						seq_reverse(p[k]->len, seq, 0); // because ->seq is reversed
					} else { // F3-R3 pairing
						if (k == 0) __set_rght_coor(beg[k], end[k], p[1-k], p[k]); // p[k] is R3
						else __set_left_coor(beg[k], end[k], p[1-k], p[k]); // p[k] is F3
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
	dbset_unload_pac(dbs);
}

void bwa_paired_sw(dbset_t *dbs, int n_seqs, bwa_seq_t *seqs[2], const pe_opt_t *popt, const isize_info_t *ii)
{
	int i;
	uint64_t n_tot[2] = {0,0};
	uint64_t n_mapped[2] = {0,0};
	bwa_paired_sw_data_t td;

	if (!popt->is_sw || ii->avg < 0.0) return;

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

//void bwa_sai2sam_pe_core(const char *prefix, char *const fn_sa[2], char *const fn_fa[2], pe_opt_t *popt)
void bwa_sai2sam_pe_core(pe_inputs_t* inputs, pe_opt_t *popt)
{
	extern bwa_seqio_t *bwa_open_reads(int mode, const char *fn_fa);
	int i, j, n_seqs, tot_seqs = 0;
	bwa_seq_t *seqs[2];
	bwa_seqio_t *ks[2];
	clock_t t;
	isize_info_t last_ii; // this is for the last batch of reads
	dbset_t *dbs = NULL;
	saiset_t *saiset = NULL;
	gap_opt_t *gopt = NULL;
	gap_opt_t *gopt0 = NULL;

	// initialization
	bwase_initialize(); // initialize g_log_n[] in bwase.c
	for (i = 1; i != 256; ++i) g_log_n[i] = (int)(4.343 * log(i) + 0.5);

	saiset = saiset_create(inputs->count, inputs->sai_pair.a);
	gopt0 = &saiset->opt[0];
	gopt = &saiset->opt[1];

	last_ii.avg = -1.0;

	ks[0] = bwa_open_reads(gopt0->mode, inputs->fq[0]);
	ks[1] = bwa_open_reads(gopt->mode, inputs->fq[1]);

	dbs = dbset_restore(inputs->count, inputs->prefixes.a, gopt->mode, popt->is_preload);
	srand48(dbs->bns[0]->bns->seed);

	// core loop
	dbset_print_sam_SQ(dbs);
	bwa_print_sam_PG();
	while ((seqs[0] = bwa_read_seq(ks[0], 0x40000, &n_seqs, gopt0->mode, gopt0->trim_qual)) != 0) {
		int cnt_chg;
		isize_info_t ii;

		seqs[1] = bwa_read_seq(ks[1], 0x40000, &n_seqs, gopt->mode, gopt->trim_qual);
		tot_seqs += n_seqs;
		t = clock();

		fprintf(stderr, "[bwa_sai2sam_pe_core] convert to sequence coordinate... \n");
		cnt_chg = bwa_cal_pac_pos_pe(dbs, n_seqs, seqs, saiset, &ii, popt, gopt, &last_ii);
		fprintf(stderr, "[bwa_sai2sam_pe_core] time elapses: %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC); t = clock();
		fprintf(stderr, "[bwa_sai2sam_pe_core] changing coordinates of %d alignments.\n", cnt_chg);

		fprintf(stderr, "[bwa_sai2sam_pe_core] align unmapped mate...\n");
		bwa_paired_sw(dbs, n_seqs, seqs, popt, &ii);
		fprintf(stderr, "[bwa_sai2sam_pe_core] time elapses: %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC); t = clock();

		fprintf(stderr, "[bwa_sai2sam_pe_core] refine gapped alignments... ");
		for (j = 0; j < 2; ++j)
			bwa_refine_gapped(dbs, n_seqs, seqs[j]);
		fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC); t = clock();

		fprintf(stderr, "[bwa_sai2sam_pe_core] print alignments... ");
		for (i = 0; i < n_seqs; ++i) {
			bwa_seq_t *p[2];
			p[0] = seqs[0] + i; p[1] = seqs[1] + i;
			if (p[0]->bc[0] || p[1]->bc[0]) {
				strcat(p[0]->bc, p[1]->bc);
				strcpy(p[1]->bc, p[0]->bc);
			}

			bwa_print_sam1(dbs, p[0], p[1], gopt->mode, gopt->max_top2);
			bwa_print_sam1(dbs, p[1], p[0], gopt->mode, gopt->max_top2);
		}
		fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC); t = clock();

		for (j = 0; j < 2; ++j)
			bwa_free_read_seq(n_seqs, seqs[j]);
		fprintf(stderr, "[bwa_sai2sam_pe_core] %d sequences have been processed.\n", tot_seqs);
		last_ii = ii;
	}

	// destroy
	dbset_destroy(dbs);
	saiset_destroy(saiset);

	for (i = 0; i < 2; ++i) {
		bwa_seq_close(ks[i]);
	}
}

static pe_inputs_t *pe_inputs_parse(int argc, char *argv[])
{
	pe_inputs_t *inputs = calloc(1, sizeof(pe_inputs_t));
	int i;
	inputs->count = 0;

	if (argc < 2) {
		fprintf(stderr, "not enough arguments!\n");
		exit(1);
	}

	inputs->fq[0] = argv[0];
	inputs->fq[1] = argv[1];

	i = 2;
	while (i < argc)
	{
		if (argc - i < 3) {
			fprintf(stderr, "[%s] insufficient arguments\n", __func__);
			exit(1);
		}

		kv_push(const char*, inputs->prefixes, argv[i++]);
		kv_push(const char**, inputs->sai_pair, (const char**)&argv[i]);
		i += 2;
		++inputs->count;
	}

	return inputs;
}

static void pe_inputs_destroy(pe_inputs_t *inputs)
{
	kv_destroy(inputs->prefixes);
	kv_destroy(inputs->sai_pair);
	free(inputs);
}

static void dump_pe_inputs(const pe_inputs_t *inputs)
{
	int i;
	fprintf(stderr, "[%s]: %d sets\n", __func__, inputs->count);
	fprintf(stderr, " - fastq files: %s, %s\n", inputs->fq[0], inputs->fq[1]);
	for (i = 0; i < inputs->count; ++i) {
		fprintf(stderr, " - ref: %s sai pair: <%s> <%s>\n",
			inputs->prefixes.a[i],
			inputs->sai_pair.a[i][0],
			inputs->sai_pair.a[i][1]);
	}
}

int bwa_sai2sam_pe(int argc, char *argv[])
{
	extern char *bwa_rg_line, *bwa_rg_id;
	extern int bwa_set_rg(const char *s);
	pe_inputs_t *inputs;
	int c;
	pe_opt_t *popt;
	popt = bwa_init_pe_opt();
	while ((c = getopt(argc, argv, "a:o:sPn:N:c:f:Ar:t:")) >= 0) {
		switch (c) {
		case 'r':
			if (bwa_set_rg(optarg) < 0) {
				fprintf(stderr, "[%s] malformated @RG line\n", __func__);
				return 1;
			}
			break;
		case 'a': popt->max_isize = atoi(optarg); break;
		case 'o': popt->max_occ = atoi(optarg); break;
		case 's': popt->is_sw = 0; break;
		case 'P': popt->is_preload = 1; break;
		case 'n': popt->n_multi = atoi(optarg); break;
		case 'N': popt->N_multi = atoi(optarg); break;
		case 't': popt->n_threads = atoi(optarg); break;
		case 'c': popt->ap_prior = atof(optarg); break;
		case 'f': xreopen(optarg, "w", stdout); break;
		case 'A': popt->force_isize = 1; break;
		default: return 1;
		}
	}

	if (optind + 5 > argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   bwa sampe [options] <in1.fq> <in2.fq> <prefix> <in1.sai> <in2.sai> [<prefix2> <in2,1.sai> <in2,2.sai> ...]\n\n");
		fprintf(stderr, "Options: -a INT   maximum insert size [%d]\n", popt->max_isize);
		fprintf(stderr, "         -o INT   maximum occurrences for one end [%d]\n", popt->max_occ);
		fprintf(stderr, "         -n INT   maximum hits to output for paired reads [%d]\n", popt->n_multi);
		fprintf(stderr, "         -N INT   maximum hits to output for discordant pairs [%d]\n", popt->N_multi);
		fprintf(stderr, "         -t INT   number of threads [%d]\n", popt->n_threads);
		fprintf(stderr, "         -c FLOAT prior of chimeric rate (lower bound) [%.1le]\n", popt->ap_prior);
		fprintf(stderr, "         -f FILE  sam file to output results to [stdout]\n");
		fprintf(stderr, "         -r STR   read group header line such as `@RG\\tID:foo\\tSM:bar' [null]\n");
		fprintf(stderr, "         -P       preload index into memory (for base-space reads only)\n");
		fprintf(stderr, "         -s       disable Smith-Waterman for the unmapped mate\n");
		fprintf(stderr, "         -A       disable insert size estimate (force -s)\n\n");
		fprintf(stderr, "Notes: 1. For SOLiD reads, <in1.fq> corresponds R3 reads and <in2.fq> to F3.\n");
		fprintf(stderr, "       2. For reads shorter than 30bp, applying a smaller -o is recommended to\n");
		fprintf(stderr, "          to get a sensible speed at the cost of pairing accuracy.\n");
		fprintf(stderr, "\n");
		return 1;
	}

	inputs = pe_inputs_parse(argc-optind, &argv[optind]);
	dump_pe_inputs(inputs);

	bwa_sai2sam_pe_core(inputs, popt);

	pe_inputs_destroy(inputs);
	free(bwa_rg_line); free(bwa_rg_id);
	free(popt);
	return 0;
}
