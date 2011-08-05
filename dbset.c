#include "dbset.h"
#include "bwaremap.h"
#include "kstring.h"
#include "utils.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

extern char *bwa_rg_line;
extern char *bwa_rg_id;
bntseq_t *bwa_open_nt(const char *prefix);

/* binary search through sequences to find the one containing pos */
int coord2idx(const dbset_t *dbs, int64_t pos) {
    int left = 0, right = dbs->count;
    int mid = 0;

    if (pos > dbs->l_pac)
        return -1;

    while (left < right) {
        mid = (left + right) >> 1;
        if (pos > dbs->db[mid]->offset) {
            if (mid == dbs->count-1) break;
            if (pos < dbs->db[mid+1]->offset) break;
            left = mid + 1;
        } else if (pos < dbs->db[mid]->offset) {
            right = mid;
        } else {
            break;
        }
    }

    assert(mid < dbs->count);
    return mid;
}

static bwtdb_t *bwtdb_load(const char *prefix) {
    bwtdb_t *db = calloc(1, sizeof(bwtdb_t));

    db->prefix = prefix;
    db->bwtcache = bwtcache_create();

    return db;
}

static void bwtdb_load_sa(bwtdb_t *db, int which) {
    char path[PATH_MAX];
    const char *bwt_suffix = ".bwt";
    const char *sa_suffix = ".sa";

    if (db->bwt[which] != NULL)
        return;

    if (which == 1) {
        bwt_suffix = ".rbwt";
        sa_suffix = ".rsa";
    }

    strcpy(path, db->prefix); strcat(path, bwt_suffix);
    db->bwt[which] = bwt_restore_bwt(path);
    strcpy(path, db->prefix); strcat(path, sa_suffix);
    bwt_restore_sa(path, db->bwt[which]);
}

static void bwtdb_unload_sa(bwtdb_t *db, int which) {
    if (db->bwt[which])
        bwt_destroy(db->bwt[which]);
    db->bwt[which] = NULL;
}

static void bwtdb_destroy(bwtdb_t *db) {
    bwtdb_unload_sa(db, 0);
    bwtdb_unload_sa(db, 1);
    bwtcache_destroy(db->bwtcache);
    free(db);
}

static seq_t *seq_restore(const char *prefix, const char *extension) {
    char path[PATH_MAX];
    strcat(strcpy(path, prefix), extension);
    seq_t *s = calloc(1, sizeof(seq_t));
    s->bns = bns_restore(path);
    s->remap = s->bns->n_seqs > 0 && can_remap(s->bns->anns[0].name);
    if (s->remap)
        fprintf(stderr, " - Remapping enabled for sequence %s\n", prefix);
    return s;
}

static void seq_load_pac(seq_t *s) {
    assert(s->data == NULL);
    s->data = (ubyte_t*)calloc(s->bns->l_pac/4+1, 1);
    rewind(s->bns->fp_pac);
    fread(s->data, 1, s->bns->l_pac/4+1, s->bns->fp_pac);
}

static void seq_unload_pac(seq_t *s) {
    if (s->data) free(s->data);
    s->data = NULL;
}

static void seq_destroy(seq_t *s) {
    if (!s) return;
    seq_unload_pac(s);
    if (s->bns) bns_destroy(s->bns);
    free(s);
}

dbset_t *dbset_restore(int count, const char **prefixes, int mode, int preload) {
    int i;
    dbset_t *dbs = calloc(1, sizeof(dbset_t));
    dbs->count = count;
    dbs->db = (bwtdb_t**)calloc(count, sizeof(bwtdb_t*));
    dbs->coord_map = calloc(count, sizeof(uint32_t*));
    dbs->bns = calloc(count, sizeof(seq_t*));
    dbs->ntbns = calloc(count, sizeof(seq_t*));
    dbs->preload = preload;

    dbs->color_space = !(mode & BWA_MODE_COMPREAD);

    for (i = 0; i < count; ++i) {
        dbs->db[i] = bwtdb_load(prefixes[i]);
        dbs->db[i]->offset = dbs->l_pac;
        dbs->bns[i] = seq_restore(prefixes[i], "");
        dbs->db[i]->bns = dbs->bns[i];
        dbs->l_pac += dbs->bns[i]->bns->l_pac;


/* TODO: get rid of these, we just need the length of the bwt */
        bwtdb_load_sa(dbs->db[i], 0);
        bwtdb_load_sa(dbs->db[i], 1);

        dbs->total_bwt_seq_len[0] += dbs->db[i]->bwt[0]->seq_len;
        dbs->total_bwt_seq_len[1] += dbs->db[i]->bwt[1]->seq_len;

        if (dbs->color_space) {
            dbs->ntbns[i] = seq_restore(prefixes[i], ".nt");
            dbs->db[i]->ntbns = dbs->ntbns[i];
            dbs->color_space = 1;
        } else if (preload) {
            seq_load_pac(dbs->bns[i]);
            bwtdb_load_sa(dbs->db[i], 0);
            bwtdb_load_sa(dbs->db[i], 1);
        }
    }

    return dbs;
}

void dbset_destroy(dbset_t *dbs) {
    int i;
    for (i = 0; i < dbs->count; ++i) {
        bwtdb_destroy(dbs->db[i]);
        seq_destroy(dbs->bns[i]);
        seq_destroy(dbs->ntbns[i]);
        if (dbs->coord_map[i])
            free(dbs->coord_map[i]);
    }
    free(dbs->db);
    free(dbs->coord_map);
    free(dbs->bns);
    free(dbs->ntbns);
    free(dbs);
}

void dbset_load_sa(dbset_t *dbs, int which) {
    int i;
    assert(which == 0 || which == 1);
    if (dbs->preload) return;
    for (i = 0; i < dbs->count; ++i) {
        if (dbs->db[i]->bwt[which] == NULL)
            bwtdb_load_sa(dbs->db[which], which);
    }
}

void dbset_unload_sa(dbset_t *dbs, int which) {
    int i;
    if (dbs->preload) return;
    for (i = 0; i < dbs->count; ++i)
        bwtdb_unload_sa(dbs->db[i], which);
}

void dbset_load_pac(dbset_t *dbs) {
    int i;
    if (dbs->preload)
        return;

    for (i = 0; i < dbs->count; ++i)
        if (dbs->bns[i]->data == NULL)
            seq_load_pac(dbs->bns[i]);
}

void dbset_unload_pac(dbset_t *dbs) {
    int i;
    if (dbs->preload)
        return;

    for (i = 0; i < dbs->count; ++i)
        seq_unload_pac(dbs->bns[i]);
}

void dbset_load_ntpac(dbset_t *dbs) {
    int i;
    for (i = 0; i < dbs->count; ++i)
        if (dbs->ntbns[i]->data == NULL)
            seq_load_pac(dbs->ntbns[i]);
}

void dbset_unload_ntpac(dbset_t *dbs) {
    int i;
    if (dbs->preload)
        return;

    for (i = 0; i < dbs->count; ++i)
        seq_unload_pac(dbs->ntbns[i]);
}

uint64_t bwtdb_sa2seq(const bwtdb_t *db, int strand, uint32_t sa, uint32_t seq_len) {
    if (strand) {
        return db->offset + bwt_sa(db->bwt[0], sa);
    } else {
        return db->offset + db->bwt[1]->seq_len - (bwt_sa(db->bwt[1], sa) + seq_len);
    }
}

int dbset_coor_pac2real(const dbset_t *dbs, int64_t pac_coor, int len, int32_t *real_seq, 
                        const bntseq_t **bns, uint64_t *offset)
{
    int idx = coord2idx(dbs, pac_coor);
    *bns = dbs->bns[idx]->bns;
    *offset = dbs->db[idx]->offset;
    return bns_coor_pac2real(*bns, pac_coor - *offset, len, real_seq);
}

poslist_t bwtdb_cached_sa2seq(const bwtdb_t *db, const bwt_aln1_t* aln, uint32_t seq_len) {
    return bwt_cached_sa(db->offset, db->bwtcache, (const bwt_t **const)db->bwt, aln, seq_len);
}

uint32_t dbset_extract_remapped(const dbset_t *dbs, seq_t **seqs, uint32_t dbidx, int32_t seqid, ubyte_t* ref_seq, uint64_t beg, uint32_t len) {
    bwtdb_t *db = dbs->db[dbidx];
    uint64_t seq_begin;
    uint32_t total = 0;
    bntann1_t *ann;

    if (seqid < 0 || !db->bns->remap) {
        return dbset_extract_sequence(dbs, seqs, ref_seq, beg, len);
    }

    ann = db->bns->bns->anns + seqid;
    seq_begin = db->offset + ann->offset;

    if (beg < seq_begin) {
        uint64_t remapped_begin = bwa_remap_position_with_seqid(db->bns->bns, ann->offset, seqid);
        uint64_t sublen = seq_begin - beg;
        uint64_t offset = remapped_begin - sublen;
        if (sublen > remapped_begin)
            err_fatal(__func__, "request too far ahead of remapped region");
        total += dbset_extract_sequence(dbs, seqs, &ref_seq[total], offset, sublen);
    }

    if (total < len) {
        uint32_t sublen = len - total;
        if (sublen > ann->len)
            sublen = ann->len;
        total += dbset_extract_sequence(dbs, seqs, &ref_seq[total], seq_begin, sublen);
    }

    if (total < len) {
        uint64_t remapped_end = bwa_remap_position_with_seqid(db->bns->bns, ann->offset + ann->len-1, seqid)+1;
        total += dbset_extract_sequence(dbs, seqs, &ref_seq[total], remapped_end, len-total);
    }

    if (total != len) {
        err_fatal(__func__, "logic error: got %lu bases instead of %lu\n", total, len);
    }

    return total;
}

uint32_t dbset_extract_sequence(const dbset_t *dbs, seq_t **seqs, ubyte_t* ref_seq, uint64_t beg, uint32_t len) {
    uint32_t total = 0;
    while (total < len) {
        int64_t idx, pos;
        seq_t *s; 
        bwtdb_t *db;
        if (beg >= dbs->l_pac) break;
        idx = coord2idx(dbs, beg);
        s = seqs[idx];
        db = dbs->db[idx];
        pos = beg - db->offset;
        /* TODO: this is wrong */
        while (pos < s->bns->l_pac && total < len) {
            ref_seq[total++] = bns_pac(s->data, pos);
            ++pos; /* bns_pac is a macro referencing idx more than once */
        }
        beg = pos + db->offset;
    }
    return total;
}

void dbset_print_sam_SQ(const dbset_t *dbs) {
    int i, j;
    for (i = 0; i < dbs->count; ++i) {
        seq_t *s = dbs->bns[i];
        for (j = 0; j < s->bns->n_seqs; ++j)
            printf("@SQ\tSN:%s\tLN:%d\n", 
                s->bns->anns[j].name, s->bns->anns[j].len);
    }
    if (bwa_rg_line) printf("%s\n", bwa_rg_line);
}
