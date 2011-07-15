#ifndef DBSET_H
#define DBSET_H

#include "bntseq.h"
#include "bwt.h"
#include "bwtaln.h"
#include "bwtcache.h"

#include <stdint.h>

typedef struct {
    bntseq_t *bns;
    ubyte_t *data;
    int remap; /* set if the sequence is logically remapped onto another*/
} seq_t;

typedef struct {
    const char *prefix;
    bwt_t *bwt[2]; 
    bwtcache_t *bwtcache;
    uint64_t offset;
    seq_t *bns;
    seq_t *ntbns;
} bwtdb_t;

typedef struct {
    int count; 
    int color_space;
    int preload;
    bwtdb_t **db;
    uint32_t **coord_map;
    seq_t **bns;
    seq_t **ntbns;
    uint64_t l_pac;
    uint64_t total_bwt_seq_len[2];
} dbset_t;


#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

    int coord2idx(const dbset_t *dbs, int64_t pos);
    dbset_t *dbset_restore(int count, const char **prefixes, int mode, int preload);
    void dbset_destroy(dbset_t *dbs);

    void dbset_load_sa(dbset_t *dbs, int which);
    void dbset_unload_sa(dbset_t *dbs, int which);

    void dbset_load_pac(dbset_t *dbs);
    void dbset_unload_pac(dbset_t *dbs);

    void dbset_load_ntpac(dbset_t *dbs);
    void dbset_unload_ntpac(dbset_t *dbs);

    uint32_t dbset_extract_sequence(const dbset_t *dbs, seq_t **seqs, ubyte_t* ref_seq, uint64_t beg, uint32_t len);
    int dbset_coor_pac2real(const dbset_t *dbs, int64_t pac_coor, int len, int32_t *real_seq, 
                            const bntseq_t **bns, uint64_t *offset);

    uint64_t bwtdb_sa2seq(const bwtdb_t *db, int strand, uint32_t sa, uint32_t seq_len);
    poslist_t bwtdb_cached_sa2seq(const bwtdb_t *db, const bwt_aln1_t* aln, uint32_t seq_len);

    void dbset_print_sam_SQ(const dbset_t *dbs);

#ifdef __cplusplus
}
#endif /* __cplusplus */


#endif /* DBSET_H */
