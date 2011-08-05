#ifndef SAISET_H
#define SAISET_H

#include <stdio.h>
#include "bntseq.h"
#include "bwtaln.h"
#include "kvec.h"
#include "dbset.h"

typedef struct {
    bwt_aln1_t aln;
    int remapped;
    uint32_t dbidx;
    bwtdb_t *db;
} alignment_t;

typedef kvec_t(alignment_t) alngrp_t;

typedef struct {
    int count;
    FILE **fp[2];
    gap_opt_t opt[2];
} saiset_t;

#ifdef __cplusplus
extern "C" {
#endif

    saiset_t *saiset_create(int n, const char **files[]);
    void saiset_destroy(saiset_t *saiset);

    alngrp_t *alngrp_create(const dbset_t *dbs, const saiset_t *saiset, int which);
    void alngrp_destroy(alngrp_t *ag);

    void select_sai(const alngrp_t *aln, bwa_seq_t *s, int *main_idx);
    void select_sai_multi(const alngrp_t *ag, bwa_seq_t *s, int n_multi);

#ifdef __cplusplus
}
#endif

#endif /* SAISET_H */
