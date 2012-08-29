#include "filter_alignments.h"
#include "bwaremap.h"
#include "bwase.h"

#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>

#define MIN_HASH_WIDTH 1000

using namespace std;

namespace {
    static uint64_t __remap(const uint64_t pos, uint64_t len, uint32_t gap, const bwtdb_t *db, const bwtdb_t *target, int32_t *seqid, int *identical) {
        uint64_t x;
        const read_mapping_t *m;
        uint64_t relpos = pos;

        if (!db->bns->remap) {/* not all sequences need remapping */
            *seqid = -1;
            return pos;
        }

        /* get the position relative to the particular sequence it is from */
        x = bwa_remap_position(db->bns, target->bns->bns, pos - db->offset, seqid);
        m = &db->bns->mappings[*seqid]->map;
        relpos = pos - db->offset - db->bns->bns->anns[*seqid].offset;
        *identical = is_remapped_sequence_identical(m, relpos > gap ? relpos - gap : 0, len + gap);

        return x;
    }
}

/* TODO: currently, the remapped dbidx is hard coded as 0, might want to change that in the future
 * to allow remappings to things other than the primary sequence */
#define remap(p, dbs, _dbidx, opt_remap) do { \
        uint64_t gap = (p)->n_gapo + (p)->n_gape; \
        uint64_t len = (p)->len; \
		const bwtdb_t *db = (dbs)->db[(_dbidx)]; \
		(p)->dbidx = (_dbidx); \
		(p)->remapped_dbidx = 0; \
		if ((opt_remap)) { \
			(p)->remapped_pos = __remap((p)->pos, len, gap, db, (dbs)->db[0], &(p)->remapped_seqid, &(p)->remap_identical); \
		} else { \
			(p)->remapped_pos = (p)->pos; \
			(p)->remapped_seqid = -1; \
		} \
	} while (0);


void compute_seq_coords_and_counts(
    const dbset_t* dbs,
    int do_remap,
    const alngrp_t aln[2],
    pos_arr_t* out_arr,
    bwa_seq_t** p
    )
{
    struct db_and_score {
        int dbidx;
        int score;
    };

    out_arr->n = 0;
    for (int j = 0; j < 2; ++j) {
        map<uint64_t, alignment_t*> pos2score;
        int min_score = numeric_limits<int>::max();
        for (unsigned k = 0; k < aln[j].n; ++k) {
            alignment_t *ar = &aln[j].a[k];
            min_score = std::min(min_score, ar->aln.score);

            bwtint_t l;
            if (ar->aln.l - ar->aln.k + 1 >= MIN_HASH_WIDTH) { // then check hash table
                bwtdb_t* db = dbs->db[ar->dbidx];
                /* TODO: cache remappings */
                poslist_t pos = bwtdb_cached_sa2seq(dbs->db[ar->dbidx], &ar->aln, p[j]->len);
                for (l = 0; l < pos.n; ++l) {
                    position_t alnpos = {0};
                    alnpos.pos = pos.a[l];
                    if (alnpos.pos < db->offset || alnpos.pos >= db->offset + db->bns->bns->l_pac)
                        continue;
                    alnpos.len = p[j]->len;
                    alnpos.n_gape = ar->aln.n_gape;
                    alnpos.n_gapo = ar->aln.n_gapo;
                    alnpos.score = ar->aln.score;
                    remap(&alnpos, dbs, ar->dbidx, do_remap);
                    alnpos.idx_and_end = k<<1 | j;
                    kv_push(position_t, *out_arr, alnpos);

                    auto inserted = pos2score.insert(make_pair(alnpos.remapped_pos, ar));
                    if (!inserted.second) {
                        if (ar->aln.score < inserted.first->second->aln.score)
                            inserted.first->second = ar;
                    }
                }
            } else { // then calculate on the fly
                bwtdb_t* db = dbs->db[ar->dbidx];
                for (l = ar->aln.k; l <= ar->aln.l; ++l) {
                    position_t alnpos = {0};
                    alnpos.pos = bwtdb_sa2seq(db, ar->aln.a, l, p[j]->len);
                    if (alnpos.pos < db->offset || alnpos.pos >= db->offset + db->bns->bns->l_pac)
                        continue;
                    alnpos.len = p[j]->len;
                    alnpos.n_gape = ar->aln.n_gape;
                    alnpos.n_gapo = ar->aln.n_gapo;
                    alnpos.score = ar->aln.score;
                    remap(&alnpos, dbs, ar->dbidx, do_remap);
                    alnpos.idx_and_end = k<<1 | j;
                    kv_push(position_t, *out_arr, alnpos);
                    auto inserted = pos2score.insert(make_pair(alnpos.remapped_pos, ar));
                    if (!inserted.second) {
                        if (ar->aln.score < inserted.first->second->aln.score)
                            inserted.first->second = ar;
                    }
                }
            }
        }

        size_t totalReadCounts[2] = {0};
        size_t primaryReadCounts[2] = {0};
        for (auto i = pos2score.begin(); i != pos2score.end(); ++i) {
            int idx = i->second->aln.score == min_score ? 0 : 1;
            ++totalReadCounts[idx];
            if (i->second->dbidx == 0)
                ++primaryReadCounts[idx];
        }

        if ((p[j]->c1 = primaryReadCounts[0]) == 0)
            p[j]->c1 = totalReadCounts[0];

        p[j]->c2 = primaryReadCounts[1];
        if (p[j]->c1 != 0)
            p[j]->type = p[j]->c1 > 1 ? BWA_TYPE_REPEAT : BWA_TYPE_UNIQUE;
    }
}
