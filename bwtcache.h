#ifndef BWTCACHE_H
#define BWTCACHE_H

#include "bwt.h"
#include "bwtaln.h"
#include "kvec.h"

#include <stdint.h>

typedef struct _bwtcache_t bwtcache_t;

typedef enum {
    eUNINITIALIZED,
    eLOADING,
    eINITIALIZED
} cache_item_state_t;

typedef struct {
    uint64_t n;
    uint64_t *a;
} poslist_t;

typedef struct {
    cache_item_state_t state;
    poslist_t pos;
} bwtcache_itm_t;

#ifdef __cplusplus
extern "C" {
#endif

    poslist_t bwt_cached_sa(uint64_t offset, bwtcache_t *c, const bwt_t *const bwt[2], const bwt_aln1_t *a, uint32_t seqlen);

    bwtcache_t *bwtcache_create();
    void bwtcache_destroy(bwtcache_t *c);

#ifdef __cplusplus
}
#endif


#endif /* BWTCACHE_H */
