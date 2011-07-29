#ifndef BWASW_H
#define BWASW_H

#include "dbset.h"
#include "bwtaln.h"
#include "bwapair.h"

#ifdef __cplusplus
extern "C" {
#endif

    void bwa_paired_sw(dbset_t *dbs, int n_seqs, bwa_seq_t *seqs[2], const pe_opt_t *popt, const isize_info_t *ii);

#ifdef __cplusplus
}
#endif

#endif /* BWASW_H */
