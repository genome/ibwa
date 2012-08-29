#include <unistd.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <string.h>
#include "bntseq.h"
#include "bwapair.h"
#include "bwaremap.h"
#include "bwasw.h"
#include "bwtaln.h"
#include "bwtcache.h"
#include "dbset.h"
#include "khash.h"
#include "kvec.h"
#include "saiset.h"
#include "stdaln.h"
#include "threadblock.h"
#include "utils.h"
#include "filter_alignments.h"

KHASH_MAP_INIT_INT64(64, bwtcache_itm_t)

typedef struct {
	int count;
	const char *fq[2];
	kvec_t(const char*) prefixes;
	kvec_t(const char**) sai_pair;
} pe_inputs_t;

typedef kvec_t(bwt_aln1_t) aln_buf_t;

typedef struct {
	dbset_t *dbs;
	const alngrp_t **buf[2];
	int n_seqs;
	bwa_seq_t *seqs[2];
	isize_info_t *ii;
	const pe_opt_t *opt;
	const gap_opt_t *gopt;
	int *cnt_chg;
} cal_pac_pos_params_t;

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
	po->remapping = 0;
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

static uint64_t __remap(const uint64_t pos, uint64_t len, int strand, uint32_t gap, const bwtdb_t *db, const bwtdb_t *target, int32_t *seqid, int *identical) {
	uint64_t x;
    const read_mapping_t *m;
    uint64_t relpos = pos;
    uint64_t shift = 0;

	if (!db->bns->remap) {/* not all sequences need remapping */
		*seqid = -1;
		return pos;
	}

	/* get the position relative to the particular sequence it is from */
	x = bwa_remap_position(db->bns, target->bns->bns, pos - db->offset+shift, seqid);
    m = &db->bns->mappings[*seqid]->map;
    relpos = pos - db->offset - db->bns->bns->anns[*seqid].offset;
    *identical = is_remapped_sequence_identical(m, relpos > gap ? relpos - gap : 0, len + gap);

	return x;
}

/* TODO: currently, the remapped dbidx is hard coded as 0, might want to change that in the future
 * to allow remappings to things other than the primary sequence */
#define remap(p, dbs, _dbidx, opt_remap) do { \
        uint64_t gap = (p)->n_gapo + (p)->n_gape; \
        uint64_t len = (p)->len; \
		const bwtdb_t *db = (dbs)->db[(_dbidx)]; \
		(p)->dbidx = (_dbidx); \
		(p)->remapped_dbidx = 0; \
		if ((opt_remap)) { \
			(p)->remapped_pos = __remap((p)->pos, len, (p)->strand, gap, db, (dbs)->db[0], &(p)->remapped_seqid, &(p)->remap_identical); \
		} else { \
			(p)->remapped_pos = (p)->pos; \
			(p)->remapped_seqid = -1; \
		} \
	} while (0);


static void bwa_cal_pac_pos_pe_thread(uint32_t idx, uint32_t size, void *data)
{
	cal_pac_pos_params_t const *tdata = (cal_pac_pos_params_t*)data;
	const dbset_t *dbs = tdata->dbs;
	const alngrp_t **buf[2] = {tdata->buf[0], tdata->buf[1]};
	int n_seqs = tdata->n_seqs;
	bwa_seq_t *seqs[2] = {tdata->seqs[0], tdata->seqs[1]};
	isize_info_t *ii = tdata->ii;
	const pe_opt_t *opt = tdata->opt;
	const gap_opt_t *gopt = tdata->gopt;
	alngrp_t aln[2] = {{0,0}, {0,0}};
	pos_arr_t arr = {0};
	int i,j;

	tdata->cnt_chg[idx] = 0;
	for (i = idx; i < n_seqs; i += size) {
		bwa_seq_t *p[2];
		for (j = 0; j < 2; ++j) {
			p[j] = seqs[j] + i;
			aln[j] = *buf[j][i];
		}

        compute_seq_coords_and_counts(dbs, opt->remapping, aln, &arr, p);
        for (j = 0; j < 2; ++j) {
            int max_diff = gopt->fnr > 0.0? bwa_cal_maxdiff(p[j]->len, BWA_AVG_ERR, gopt->fnr) : gopt->max_diff;
            if (p[j]->c1 || p[j]->c2)
                p[j]->seQ = p[j]->mapQ = bwa_approx_mapQ(p[j], max_diff);
        }

		if ((p[0]->type == BWA_TYPE_UNIQUE || p[0]->type == BWA_TYPE_REPEAT)
			&& (p[1]->type == BWA_TYPE_UNIQUE || p[1]->type == BWA_TYPE_REPEAT))
		{ // only when both ends mapped
/*
            int k, n_occ[2];
			for (j = 0; j < 2; ++j) {
				n_occ[j] = 0;
				for (k = 0; k < aln[j].n; ++k)
					n_occ[j] += aln[j].a[k].aln.l - aln[j].a[k].aln.k + 1;
			}
*/

			{
				pairing_param_t pairing_param = { p, &arr, aln, opt, gopt->s_mm, ii };
				tdata->cnt_chg[idx] += find_optimal_pair(&pairing_param);
			}
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
	alngrp_t **aln_buf[2];
	cal_pac_pos_params_t tp;

	aln_buf[0] = (alngrp_t**)calloc(n_seqs, sizeof(alngrp_t*));
	aln_buf[1] = (alngrp_t**)calloc(n_seqs, sizeof(alngrp_t*));

	// SE
	for (i = 0; i != n_seqs; ++i) {
		bwa_seq_t *p[2];
		for (j = 0; j < 2; ++j) {
			int main_idx = 0;
			p[j] = seqs[j] + i;
			p[j]->n_multi = 0;
			p[j]->extra_flag |= SAM_FPD | (j == 0? SAM_FR1 : SAM_FR2);
			aln_buf[j][i] = alngrp_create(dbs, saiset, j);

			// generate SE alignment and mapping quality
			select_sai(aln_buf[j][i], p[j], &main_idx);

			if (p[j]->type == BWA_TYPE_UNIQUE || p[j]->type == BWA_TYPE_REPEAT) {
				alignment_t *main_aln = &aln_buf[j][i]->a[main_idx];
				int max_diff = gopt->fnr > 0.0? bwa_cal_maxdiff(p[j]->len, BWA_AVG_ERR, gopt->fnr) : gopt->max_diff;
				p[j]->pos = bwtdb_sa2seq(dbs->db[main_aln->dbidx], p[j]->strand, p[j]->sa, p[j]->len);
				remap(p[j], dbs, main_aln->dbidx, opt->remapping);
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
		tp.buf[i] = (const alngrp_t**)aln_buf[i];
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
	for (i = 0; i < n_seqs; ++i) {
		alngrp_destroy(aln_buf[0][i]);
		alngrp_destroy(aln_buf[1][i]);
	}
	free(aln_buf[0]); free(aln_buf[1]);
	return cnt_chg;
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

	dbs = dbset_restore(inputs->count, inputs->prefixes.a, gopt->mode, popt->is_preload, popt->remapping);
	srand48(dbs->db[0]->bns->bns->seed);

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
		for (j = 0; j < 2; ++j) {
			bwa_refine_gapped(dbs, n_seqs, seqs[j]);
            /* refine_gapped changes pos, so we might need to update remapped_pos */
            for (i = 0; i < n_seqs; ++i) {
				remap(&seqs[j][i], dbs, seqs[j][i].dbidx, popt->remapping);
            }
        }

		fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC); t = clock();

		fprintf(stderr, "[bwa_sai2sam_pe_core] print alignments... ");
		for (i = 0; i < n_seqs; ++i) {
			bwa_seq_t *p[2];
			p[0] = seqs[0] + i; p[1] = seqs[1] + i;
			if (p[0]->bc[0] || p[1]->bc[0]) {
				strcat(p[0]->bc, p[1]->bc);
				strcpy(p[1]->bc, p[0]->bc);
			}

			// use remapped coords for printing
			if (popt->remapping) {
                uint64_t tmp;
				tmp = p[0]->pos; p[0]->pos = p[0]->remapped_pos; p[0]->remapped_pos = tmp;
				tmp = p[1]->pos; p[1]->pos = p[1]->remapped_pos; p[1]->remapped_pos = tmp;
			} else {
				p[0]->remapped_pos = p[0]->pos;
				p[1]->remapped_pos = p[1]->pos;
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
	int i = 0;

	if (argc < 2) {
		fprintf(stderr, "not enough arguments!\n");
		exit(1);
	}


	kv_push(const char*, inputs->prefixes, argv[i++]);
	kv_push(const char**, inputs->sai_pair, (const char**)&argv[i]);
	i += 2;

	inputs->fq[0] = argv[i++];
	inputs->fq[1] = argv[i++];

	inputs->count=1;
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
	fprintf(stderr, " - unaligned read files: %s, %s\n", inputs->fq[0], inputs->fq[1]);
	for (i = 0; i < inputs->count; ++i) {
		fprintf(stderr, " - ref: %s, sai pair: <%s> <%s>\n",
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
	while ((c = getopt(argc, argv, "a:o:sPn:N:c:f:ARr:t:")) >= 0) {
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
		case 'R': popt->remapping = 1; break;
		default: return 1;
		}
	}

	if (optind + 5 > argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   bwa sampe [options] <prefix> <in1.sai> <in2.sai> <in1.fq> <in2.fq> "
						"[<prefix2> <in2,1.sai> <in2,2.sai> <prefix3> ...]\n\n");
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
		fprintf(stderr, "         -R       enable compound sequence remapping\n");
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
