#include "bwaremap.h"
#include "byteorder.h"
#include "utils.h"

#include <errno.h>
#include <limits.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>


static int can_remap(const char* str) {
    int ndash = 0;
    int npipe = 0;
    while (*str) {
        if (*str == '-') ++ndash;
        if (*str == '|') ++npipe;
        ++str;
    }
    return ndash == 1 && npipe == 2;
}

static int cigar_gap_opens(const char* cigar) {
    int rv = 0;
    while (*cigar) {
        if (*cigar == 'I' || *cigar == 'D' || *cigar == 'N')
            ++rv;
        ++cigar;
    }
    return rv;
}

/* return value:
 *  -1 - error processing remapping file (fatal error)
 *   0 - no remapping file (this is not an error)
 *   1 - remappings loaded
 */
int load_remappings(seq_t* seq, const char* path) {
    const int bufsz = 65536;
    char buf[bufsz];
    FILE* fp = fopen(path, "r");
    long line = 0;
    int i = 0;
    char* rv = NULL;
    if (fp == NULL) {
        fprintf(stderr, "No remapping file %s: (%s)\n",
            path, strerror(errno));
        return 0;
    }

    seq->mappings = calloc(seq->bns->n_seqs, sizeof(bnsremap_t*));
    rv = fgets(buf, bufsz, fp);
    while (!feof(fp) && rv != NULL) {
        size_t cigarBufSize = bufsz;
        size_t cigarLen = 0;
        char* nlpos = strchr(buf, '\n');
        if (nlpos == NULL ) {
            fprintf(stderr, "Line %ld too long in %s\n", line, path);
            return -1;
        }
        *nlpos = 0;

        ++line;
        if (buf[0] != '>') {
            fprintf(stderr, 
                "Unexpected character '%c' at start of line %ld in "
                "file '%s'. Expected '>'.\n", buf[0], line, path);
            return -1;
        }

        /* TODO: check for OOM */
        seq->mappings[i] = calloc(1, sizeof(bnsremap_t));
        seq->mappings[i]->map.cigar = calloc(cigarBufSize, sizeof(char));
        if (!read_mapping_extract(buf+1, &seq->mappings[i]->map)) {
            fprintf(stderr, "Failed to extract read mapping from string '%s' "
                "in file %s, line %ld\n", buf+1, path, line);
            return -1;
        }

        while (!feof(fp) && (rv = fgets(buf, bufsz, fp)) != NULL && buf[0] != '>') {
            size_t len = strlen(buf);
            nlpos = strchr(buf, '\n');
            if (nlpos == NULL ) {
                fprintf(stderr, "Line %ld too long in %s\n", line, path);
                return -1;
            }
            *nlpos = 0;

            ++line;
            if (cigarLen + len > cigarBufSize) {
                while (cigarLen + len > cigarBufSize)
                    cigarBufSize *= 2;

                seq->mappings[i]->map.cigar = realloc(seq->mappings[i]->map.cigar, cigarBufSize);
                if (seq->mappings[i]->map.cigar == NULL) {
                    fprintf(stderr, "load_remappings: failed to allocate memory for cigar string\n");
                    return -1;
                }
            }

            strcat(seq->mappings[i]->map.cigar, buf);
        }
        ++line;
        seq->mappings[i]->map.n_gapo = cigar_gap_opens(seq->mappings[i]->map.cigar);
        ++i;
    }

    fclose(fp);
    return 1;
}

int read_mapping_extract(const char *str, read_mapping_t *m) {
    const char *beg;
    char *end;

    m->exact = 0;
    if (!can_remap(str))
        return 0;

    beg = strchr(str, '-'); 
    if (beg == 0)
        return 0;
    ++beg;
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
        return 1;
    }


    m->start = strtoul(beg, &end, 10);
    if (end == beg || *end != '|')
        return 0;
    --m->start;

    beg = end+1;
    m->stop = strtoul(beg, &end, 10) + 1; /* stop = one past the last base */
    if (end == beg || *end != '\0')
        return 0;

    return 1;
}

void read_mapping_destroy(read_mapping_t *m) {
    free(m->seqname);
    if (m->cigar)
        free(m->cigar);
}

int is_remapped_sequence_identical(const read_mapping_t *m, uint32_t start, uint32_t len) {
    uint32_t pos = 0;
    const char* cigar = m->cigar;
    char last_op = 0;
    uint32_t last_len = 0;
    if (m->exact)
        return 1;

    while (pos <= start && *cigar) {
        char* end;
        last_len = strtoul(cigar, &end, 10);

        if (end == cigar) {
            fprintf(stderr, "[remap_coordinates] expected number in cigar string '%s' at pos %ld\n",
                m->cigar, cigar-m->cigar);
            return 0;
        }

        cigar = end;
        last_op = *cigar;
        switch (last_op) {
        case 'M':
        case 'X':
        case '=':
        case 'N':
        case 'D':
            pos += last_len;
            break;
        case 'I':
            break;
        default:
            fprintf(stderr, "invalid cigar character '%c'\n", last_op);
            return 0;
        }
        cigar++;
    }

    if (pos > start) {
        return (last_op == 'M' || last_op == '=')
            && last_len - start > len;
    } else if (pos == last_len) {
        fprintf(stderr, "failed to parse cigar string '%s'\n", m->cigar);
        return 0;
    }

    return 0;
}

int remap_cigar(const char *cigar, uint32_t *result, uint32_t pos, uint32_t seqlen) {
    const char *p_cigar = cigar;
    uint32_t altpos = 0;
    uint32_t refpos = 0;
    char last_op = 0;
    uint32_t last_len = 0;

    if (pos >= seqlen) {
        fprintf(stderr, "[remap_coordinates] requested pos %u > sequence length %u\n", pos, seqlen);
        return 0;
    }

    while (altpos <= pos && *p_cigar) {
        char* end;
        last_len = strtoul(p_cigar, &end, 10);

        if (end == p_cigar) {
            fprintf(stderr, "[remap_coordinates] expected number in cigar string '%s' at pos %ld\n",
                cigar, p_cigar - cigar);
            return 0;
        }

        p_cigar = end;
        last_op = *p_cigar;
        switch (last_op) {
            case 'M':
            case 'X':
            case '=':
                refpos += last_len;
                altpos += last_len;
                break;

            case 'N':
            case 'D':
                refpos += last_len;
                break;

            case 'I':
                altpos += last_len;
                break;

            default:
                fprintf(stderr, "invalid cigar character '%c'\n", last_op);
                return 0;
        }
        p_cigar++;
    }

    if (altpos > seqlen) {
        fprintf(stderr, "[remap_coordinates] cigar '%s' string implies length > read mapping (%u vs %u)\n",
            cigar, altpos, seqlen);
        return 0;
    }

    if (altpos == pos) {
        *result = refpos;
        return 1;
    } else if (altpos > pos) {
        switch (last_op) {
            case 'M':
            case 'X':
            case '=':
                *result = refpos - (altpos-pos);
                break;

            case 'I':
                *result = refpos;
                break;

            default:
                fprintf(stderr, "Error remapping cigar string %s:, pos=%u\n", cigar, pos);
                return 0;
                break;
        }
    } else {
        fprintf(stderr, "failed to parse cigar string '%s'\n", cigar);
        return 0;
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
    if (!m)
        err_fatal(__func__, "No read mapping for sequence id %d\n", seqid);

    target_idx = bns_seq_by_name(tgtbns, m->seqname);
    if (target_idx < 0)
        err_fatal(__func__, "Failed to locate remapping target: %s\n", m->seqname);

    /* TODO: get rid of remap_read_coordinates using an array! */
    if (!m->exact) {
        uint32_t offset = 0;
        uint32_t altpos = pac_coor - bns->bns->anns[seqid].offset;
        if (!remap_cigar(m->cigar, &offset, altpos, bns->bns->anns[seqid].len))
            err_fatal(__func__, "Failed to remap coordinates");

        rv = m->start + offset;
    } else {
        rv = pac_coor - bns->bns->anns[seqid].offset;
    }

    if (!m->exact && (rv < m->start || rv > m->stop)) {
        err_fatal(__func__, "remapped position out of range (%u should be in [%u, %u])\n", rv, m->start, m->stop);
    }

    return rv + tgtbns->anns[target_idx].offset;
}
