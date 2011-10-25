#include "bwaremap.h"
#include "byteorder.h"
#include "utils.h"

#include <limits.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

int can_remap(const char* str) {
    return strncmp("REMAP-", str, 6) == 0;
}

int read_mapping_extract(const char *str, read_mapping_t *m) {
    const char *beg;
    char *end;

    m->exact = 0;
    if (!can_remap(str))
        return 0;

    beg = str + 6;
    end = strchr(beg, '|');
    if (beg == end || end == 0)
        return 0;

    m->seqname = malloc(end-beg+1);
    strncpy(m->seqname, beg, end-beg);
    m->seqname[end-beg]=0;
    beg = end+1;

    if (strncmp("exact", beg, 5) == 0) {
        m->exact = 1;
        m->start = 0;
        m->stop = 0;
        m->cigar = 0;
        m->var_start = 0;
        m->var_stop = 0;
        return 1;
    }


    m->start = strtoul(beg, &end, 10);
    if (end == beg || *end != '|')
        return 0;

    beg = end+1;
    m->stop = strtoul(beg, &end, 10) + 1; /* stop = one past the last base */
    if (end == beg || *end != '|')
        return 0;

    beg = end+1;
    end = strchr(beg, '|');
    m->cigar = malloc(end-beg+1);
    strncpy(m->cigar, beg, end-beg);
    m->cigar[end-beg]=0;

    beg = end+1;
    end = strchr(beg, '|');
    m->var_start = strtoul(beg, &end, 10);
    if (end == beg || *end != '|')
        return 0;

    beg = end+1;
    end = strchr(beg, '|');
    m->var_stop = strtoul(beg, &end, 10);
    if (end == beg || *end != '|')
        return 0;

    return 1;
}

void read_mapping_destroy(read_mapping_t *m) {
    free(m->seqname);
    if (m->cigar)
        free(m->cigar);
}

int remap_read_coordinates(const read_mapping_t *m, uint32_t *remapped, uint32_t len) {
    uint32_t i;
    uint32_t pos = m->start;
    const char* cigar = m->cigar;
    uint32_t *cur = remapped;

    if (m->exact) {
        for (i = 0; i < len; ++i)
            *cur++ = pos++;
        return 1;
    }

    while (*cigar) {
        char* end;
        uint32_t n = strtoul(cigar, &end, 10);

        if (end == cigar) {
            fprintf(stderr, "[remap_coordinates] expected number in cigar string '%s' at pos %ld\n",
                m->cigar, cigar-m->cigar);
            return 0;
        }

        if (cur-remapped+n > len) {
            fprintf(stderr, "[remap_coordinates] cigar '%s' string implies length > read mapping (%ld vs %d)\n",
                m->cigar, cur-remapped+n, len);
            return 0;
        }
 
        cigar = end;

        switch (*cigar) {
        case 'M':
        case 'X':
        case '=':
            for (i = 0; i < n; ++i) *cur++ = pos++;
            break;
        case 'I':
            for (i = 0; i < n; ++i) *cur++ = pos;
            break;
        case 'D':
            pos += n;
            break;
        default:
            fprintf(stderr, "[remap_coordinates] unknown cigar op '%c' in cigar string '%s'\n",
                *cigar, m->cigar);
            return 0;
            break;
        } 
        cigar++;
    }

    return 1;
}

uint64_t bwa_remap_position(const seq_t* bns, const bntseq_t* tgtbns, uint64_t pac_coor, int32_t* seqid) {
    *seqid = bns_seq_for_pos(bns->bns, pac_coor);
    return bwa_remap_position_with_seqid(bns, tgtbns, pac_coor, *seqid); 
}

/* looking up the seqid is expensive. bwa_remap_position allows the value to be stored and then
 * ..._with_seqid can be used later. */
uint64_t bwa_remap_position_with_seqid(const seq_t* bns, const bntseq_t* tgtbns, uint64_t pac_coor, int32_t seqid) {
    uint32_t rv = 0;
    int32_t target_idx = 0;
    const read_mapping_t *m = &bns->mappings[seqid]->map;
    uint32_t *map = calloc(bns->bns->anns[seqid].len, sizeof(uint32_t));
    if (!m)
        err_fatal(__func__, "No read mapping for sequence id %d\n", seqid);

    if (!remap_read_coordinates(m, map, bns->bns->anns[seqid].len))
        err_fatal(__func__, "Failed to remap coordinates");

    target_idx = bns_seq_by_name(tgtbns, m->seqname);
    if (target_idx < 0)
        err_fatal(__func__, "Failed to locate remapping target: %s\n", m->seqname);

    rv = map[pac_coor - bns->bns->anns[seqid].offset];
    free(map);
    if (rv < m->start || rv > m->stop)
        err_fatal(__func__, "remapped position out of range (%u should be in [%u, %u])\n", rv, m->start, m->stop);

    return rv + tgtbns->anns[target_idx].offset;
}
