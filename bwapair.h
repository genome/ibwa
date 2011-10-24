#ifndef BWAPAIR_H
#define BWAPAIR_H

#include "saiset.h"
#include "dbset.h"
#include "kvec.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
	double avg, std, ap_prior;
	bwtint_t low, high, high_bayesian;
} isize_info_t;

typedef struct {
	uint64_t pos;
    uint64_t remapped_pos;
	uint32_t idx_and_end;
    uint32_t dbidx;
    uint32_t remapped_dbidx;
	int32_t remapped_seqid;
} position_t;
/* Find the alignment object given a position_t and array of alignments */
#define __aln_end(x)    ((x).idx_and_end&1)
#define __aln_idx(x)    ((x).idx_and_end>>1)
#define __aln(x, aln) (aln[__aln_end(x)].a[__aln_idx(x)])

typedef kvec_t(position_t) pos_arr_t;

typedef struct {
    bwa_seq_t **p;
    const pos_arr_t *arr;
    const alngrp_t *aln;
    const pe_opt_t *opt;
    int s_mm;
    const isize_info_t *ii;
} pairing_param_t;

int find_optimal_pair(const pairing_param_t *p);

#ifdef __cplusplus
}
#endif

#endif /* BWAPAIR_H */
