#include "bwapair.h"

#include "ksort.h"
#include <math.h>

typedef struct {
} pairing_internals_t;

extern int g_log_n[256]; // in bwase.c

static int position_lt(const position_t a, const position_t b) {
	return a.pos < b.pos;
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

	if ((a->pos == b->pos) && 
		(a->idx_and_end&1) == (b->idx_and_end&1) &&
		aln1->db != aln2->db
	) {
		return 1;
	}
	return 0;
}

static inline position_t *select_mapping(const alngrp_t aln[2], const pos_arr_t *arr, int begin, int end) {
	int i;
	position_t *best = &arr->a[begin];

	for (i = begin+1; i <= end; ++i) {
		position_t *p = &arr->a[i];
		if (__aln(*p, aln).aln.score > __aln(*best, aln).aln.score)
			best = p;
	}

	return best;
}

int find_optimal_pair(const pairing_param_t *param) {
    /* for convenience */
    bwa_seq_t **p = param->p;
    const pos_arr_t *arr = param->arr;
    const alngrp_t *aln = param->aln;
    const pe_opt_t *opt = param->opt;
    int s_mm = param->s_mm;
    const isize_info_t *ii = param->ii;

	int i, j, o_n, subo_n, cnt_chg = 0, low_bound = ii->low, max_len;
	position_t last_pos[2][2], o_pos[2]; /* TODO: make these into arrays of pointers */
	uint64_t subo_score, o_score;
	max_len = p[0]->full_len;
	if (max_len < p[1]->full_len) max_len = p[1]->full_len;
	if (low_bound < max_len) low_bound = max_len;

	// here v>=u. When ii is set, we check insert size with ii; otherwise with opt->max_isize

#define __pairing_aux(u,v) do {											\
		bwtint_t l = (v).pos + p[__aln_end(v)]->len - (u).pos;			\
		if ((u).pos != (uint64_t)-1 && (v).pos > (u).pos && l >= max_len\
			&& ((ii->high && l <= ii->high_bayesian) 					\
			|| (ii->high == 0 && l <= opt->max_isize))) 				\
		{																\
			uint64_t s = __aln(v, aln).aln.score + __aln(u, aln).aln.score;		\
			s *= 10;													\
			if (ii->high) s += (int)(-4.343 * log(.5 * erfc(M_SQRT1_2 * fabs(l - ii->avg) / ii->std)) + .499); \
			s = s<<32 | (uint32_t)hash_64((u).pos<<32 | (v).pos);		\
			if (s>>32 == o_score>>32) ++o_n;							\
			else if (s>>32 < o_score>>32) { subo_n += o_n; o_n = 1; }	\
			else ++subo_n;												\
			if (s < o_score) subo_score = o_score, o_score = s, o_pos[__aln_end(u)] = (u), o_pos[__aln_end(v)] = (v); \
			else if (s < subo_score) subo_score = s;					\
		}																\
	} while (0)

#define __pairing_aux2(q, w) do {										\
		const bwt_aln1_t *r = &__aln(w, aln).aln;						\
		(q)->extra_flag |= SAM_FPP;										\
		if ((q)->pos != (w).pos || (q)->strand != r->a) {				\
			(q)->n_mm = r->n_mm; (q)->n_gapo = r->n_gapo; (q)->n_gape = r->n_gape; (q)->strand = r->a; \
			(q)->score = r->score;										\
			(q)->pos = (w).pos;											\
			if ((q)->mapQ > 0) ++cnt_chg;								\
		}																\
	} while (0)

	o_score = subo_score = (uint64_t)-1;
	o_n = subo_n = 0;
	ks_introsort(position, arr->n, arr->a);
	for (j = 0; j < 2; ++j) last_pos[j][0].pos = last_pos[j][1].pos = (uint64_t)-1;
	if (opt->type == BWA_PET_STD) {
		i = 0;
		while (i < arr->n) {
			position_t *pos = &arr->a[i];
			int strand = aln[pos->idx_and_end&1].a[pos->idx_and_end>>1].aln.a;

			if (i < arr->n-1) {
				int k = i;
				while (mappings_overlap(pos, &arr->a[k+1], aln))
					k++;
				if (k > i) {
					pos = select_mapping(aln, arr, i, k);
					i = k;
				} else ++i;
			} else {
				++i;
			}

			if (strand == 1) { // reverse strand, then check
				int y = 1 - (pos->idx_and_end&1);
				__pairing_aux(last_pos[y][1], *pos);
				__pairing_aux(last_pos[y][0], *pos);
			} else { // forward strand, then push
				last_pos[pos->idx_and_end&1][0] = last_pos[pos->idx_and_end&1][1];
				last_pos[pos->idx_and_end&1][1] = *pos;
			}
		}
	} else if (opt->type == BWA_PET_SOLID) {
/*
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
*/
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
		rr[0] = __aln(o_pos[0], aln).aln.a;
		rr[1] = __aln(o_pos[1], aln).aln.a;
		if ((p[0]->pos == o_pos[0].pos && p[0]->strand == rr[0]) && (p[1]->pos == o_pos[1].pos && p[1]->strand == rr[1])) { // both ends not moved
			if (p[0]->mapQ > 0 && p[1]->mapQ > 0) {
				int mapQ = p[0]->mapQ + p[1]->mapQ;
				if (mapQ > 60) mapQ = 60;
				p[0]->mapQ = p[1]->mapQ = mapQ;
			} else {
				if (p[0]->mapQ == 0) p[0]->mapQ = (mapQ_p + 7 < p[1]->mapQ)? mapQ_p + 7 : p[1]->mapQ;
				if (p[1]->mapQ == 0) p[1]->mapQ = (mapQ_p + 7 < p[0]->mapQ)? mapQ_p + 7 : p[0]->mapQ;
			}
		} else if (p[0]->pos == o_pos[0].pos && p[0]->strand == rr[0]) { // [1] moved
			p[1]->seQ = 0; p[1]->mapQ = p[0]->mapQ;
			if (p[1]->mapQ > mapQ_p) p[1]->mapQ = mapQ_p;
		} else if (p[1]->pos == o_pos[1].pos && p[1]->strand == rr[1]) { // [0] moved
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
