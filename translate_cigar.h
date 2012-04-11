#pragma once

#include "bwtaln.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

bwa_cigar_t* translate_cigar(
    const char* cigar,
    uint32_t start,
    bwa_cigar_t* read_cigar,
    int n_cigar,
    int read_len,
    int* n_cigar_out);

#ifdef __cplusplus
}
#endif /* __cplusplus */
