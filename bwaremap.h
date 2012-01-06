#ifndef __coordmap_h__
#define __coordmap_h__

#include "bntseq.h"

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

    int can_remap(const char* str);
    int read_mapping_extract(const char *str, read_mapping_t *m);
    void read_mapping_destroy(read_mapping_t *m);
    int is_remapped_sequence_identical(const read_mapping_t *m, uint32_t start, uint32_t len);
    int remap_read_coordinates(const read_mapping_t *m, uint32_t *remapped, uint32_t len);
    void bwa_remap_dump(uint32_t len, uint32_t *map, const char *outfile);
    void bwa_remap_load(uint32_t *len, uint32_t **map, const char *infile);
    uint64_t bwa_remap_position(const seq_t* bns, const bntseq_t* target, uint64_t pac_coor, int32_t *seqid);
    uint64_t bwa_remap_position_with_seqid(const seq_t* bns, const bntseq_t* target, uint64_t pac_coor, int32_t seqid);

#ifdef __cplusplus
}
#endif

#endif /* __coordmap_h__ */
