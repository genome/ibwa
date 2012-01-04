#ifndef __filter_alignments_h__
#define  __filter_alignments_h__

#include "dbset.h"
#include "bwtaln.h"
#include "saiset.h"
#include "bwapair.h"

#ifdef __cplusplus
extern "C" {
#endif

    void compute_seq_coords_and_counts(
        const dbset_t* dbs,
        int do_remap,
        const alngrp_t aln[2],
        pos_arr_t* out_arr,
        bwa_seq_t** p
        );

#ifdef __cplusplus
}
#endif

#endif /*  __filter_alignments_h__ */
