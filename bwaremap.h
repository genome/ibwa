#ifndef __coordmap_h__
#define __coordmap_h__

#include "bntseq.h"

#include <stdint.h>

typedef struct
{
    char *seqname;
    uint32_t start;
    uint32_t stop;
    char *cigar;
} read_mapping_t;

int can_remap(const char* str);
int read_mapping_extract(const char *str, read_mapping_t *m);
void read_mapping_destroy(read_mapping_t *m);
int remap_read_coordinates(const read_mapping_t *m, uint32_t *remapped, uint32_t len);
uint32_t *bwa_remap_bns(const bntseq_t* bns);
void bwa_remap_dump(uint32_t len, uint32_t *map, const char *outfile);
void bwa_remap_load(uint32_t *len, uint32_t **map, const char *infile);
uint32_t bwa_remap_position(const bntseq_t* bns, uint64_t pac_coor);

#endif /* __coordmap_h__ */
