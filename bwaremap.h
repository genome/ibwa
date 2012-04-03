#ifndef __coordmap_h__
#define __coordmap_h__

#include "bntseq.h"
#include "bwtaln.h"

#include <stdint.h>

typedef struct
{
    char *seqname;
    int exact;
    uint32_t start;
    uint32_t stop;
    char *cigar;
} read_mapping_t;

typedef struct {
    read_mapping_t map;
    uint64_t target_tid_offset;
} bnsremap_t;

typedef struct {
    bntseq_t *bns;
    ubyte_t *data;
    int remap; /* set if the sequence is logically remapped onto another*/
    bnsremap_t **mappings;
} seq_t;

#ifdef __cplusplus
extern "C" {
#endif

    int load_remappings(seq_t* seq, const char* path);


    int read_mapping_extract(const char *str, read_mapping_t *m);
    void read_mapping_destroy(read_mapping_t *m);
    int is_remapped_sequence_identical(const read_mapping_t *m, uint32_t start, uint32_t len);
    // int remap_read_coordinates(const read_mapping_t *m, uint32_t *remapped, uint32_t len);
    int remap_cigar(const char *cigar, uint32_t *result, uint32_t pos, uint32_t seqlen);

    uint64_t bwa_remap_position(const seq_t* bns, const bntseq_t* target, uint64_t pac_coor, int32_t *seqid);
    uint64_t bwa_remap_position_with_seqid(const seq_t* bns, const bntseq_t* target, uint64_t pac_coor, int32_t seqid);
    int update_cigar(bwa_seq_t* seq, uint32_t pos, const seq_t* bns);

#ifdef __cplusplus
}
#endif

#endif /* __coordmap_h__ */
