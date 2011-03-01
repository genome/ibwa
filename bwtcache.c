#include "bwtcache.h"

#include "utils.h"
#include "khash.h"

#include <pthread.h>

#define psafe(expr, msg) xassert((expr)==0, msg)

KHASH_MAP_INIT_INT64(64, bwtcache_itm_t)

struct _bwtcache_t {
    kh_64_t *hash;

#ifdef HAVE_PTHREAD
    pthread_mutex_t mtx;
    pthread_mutex_t cond_mtx;
    pthread_cond_t cond;
#endif /* HAVE_PTHREAD */

    uint64_t cache_waits;
};

bwtcache_t *bwtcache_create() {
    bwtcache_t *c = calloc(1, sizeof(bwtcache_t));
    c->hash = kh_init(64);
 
#ifdef HAVE_PTHREAD
    psafe(pthread_mutex_init(&c->mtx, NULL), "failed to initialize cache mutex");
    psafe(pthread_mutex_init(&c->cond_mtx, NULL), "failed to initialize cache cond mutex");
    psafe(pthread_cond_init(&c->cond, NULL), "failed to initialize condition variable");
#endif /* HAVE_PTHREAD */

    return c;
}

bwtcache_itm_t bwt_cached_sa(bwtcache_t *c, const bwt_t *bwt[2], const bwt_aln1_t *a, uint32_t seqlen) {
    bwtint_t l;
    bwtcache_itm_t itm;
    uint64_t key = (uint64_t)a->k<<32 | a->l;
    itm = bwtcache_get(c, key);
    if (itm.state == eUNINITIALIZED) {
        itm.n = a->l - a->k + 1;
        itm.a = (bwtint_t*)malloc(sizeof(bwtint_t) * itm.n);
        for (l = a->k; l <= a->l; ++l)
            itm.a[l - a->k] = a->a? bwt_sa(bwt[0], l) : bwt[1]->seq_len - (bwt_sa(bwt[1], l) + seqlen);
        bwtcache_put(c, key, &itm);
    } else if (itm.state == eLOADING) {
        return bwtcache_wait(c, key);
    }

    return itm; 
} 

void bwtcache_destroy(bwtcache_t *c) {
	khint_t iter;

#ifdef HAVE_PTHREAD
    psafe(pthread_mutex_destroy(&c->mtx), "failed to destroy mutex");
    psafe(pthread_mutex_destroy(&c->cond_mtx), "failed to destroy mutex");
    psafe(pthread_cond_destroy(&c->cond), "failed to destroy condition variable");
#endif /* HAVE_PTHREAD */

    fprintf(stderr, "[%s] %lu cache waits encountered\n", __func__, c->cache_waits);
	for (iter = kh_begin(c->hash); iter != kh_end(c->hash); ++iter)
		if (kh_exist(c->hash, iter)) free(kh_val(c->hash, iter).a);
	kh_destroy(64, c->hash);
    free(c);
}

bwtcache_itm_t bwtcache_get(bwtcache_t* c, uint64_t key) {
    khint_t iter;
    int ret;
    bwtcache_itm_t *item; 
    bwtcache_itm_t rv; 

#ifdef HAVE_PTHREAD
    psafe(pthread_mutex_lock(&c->mtx), "failed to lock mutex");
#endif /* HAVE_PTHREAD */

    iter = kh_put(64, c->hash, key, &ret);
    item = &kh_val(c->hash, iter);
    rv = *item;
    if (ret) {
        item->state = eLOADING;
        rv.state = eUNINITIALIZED;
    } else {
        rv = kh_val(c->hash, iter);
    }
#ifdef HAVE_PTHREAD
    psafe(pthread_mutex_unlock(&c->mtx), "failed to unlock mutex");    
#endif /* HAVE_PTHREAD */
    return rv;
}

void bwtcache_put(bwtcache_t *c, uint64_t key, bwtcache_itm_t *value) {
    khint_t iter;
    int ret;
    bwtcache_itm_t *item; 
#ifdef HAVE_PTHREAD
    psafe(pthread_mutex_lock(&c->mtx), "failed to lock mutex");
#endif /* HAVE_PTHREAD */
    iter = kh_put(64, c->hash, key, &ret);
    item = &kh_val(c->hash, iter);
    if (item->state == eINITIALIZED) {
        free(value->a);
    } else {
        psafe(pthread_mutex_lock(&c->cond_mtx), "failed to lock mutex");
        *item = *value;
        item->state = eINITIALIZED;
        pthread_cond_broadcast(&c->cond);
        psafe(pthread_mutex_unlock(&c->cond_mtx), "failed to unlock mutex");    
    }
#ifdef HAVE_PTHREAD
    psafe(pthread_mutex_unlock(&c->mtx), "failed to unlock mutex");    
#endif /* HAVE_PTHREAD */
}

bwtcache_itm_t bwtcache_wait(bwtcache_t *c, uint64_t key) {
    khint_t iter;
    bwtcache_itm_t *item;
    int ret;
    ++c->cache_waits;
#ifdef HAVE_PTHREAD
    psafe(pthread_mutex_lock(&c->cond_mtx), "failed to lock mutex");
#endif /* HAVE_PTHREAD */
    iter = kh_put(64, c->hash, key, &ret);
    item = &kh_val(c->hash, iter);
    while (item->state != eINITIALIZED) {
#ifdef HAVE_PTHREAD
        psafe(pthread_cond_wait(&c->cond, &c->cond_mtx), "failed to wait on condition variable");
#endif /* HAVE_PTHREAD */
        iter = kh_put(64, c->hash, key, &ret);
        item = &kh_val(c->hash, iter);
    }
#ifdef HAVE_PTHREAD
    psafe(pthread_mutex_unlock(&c->cond_mtx), "failed to unlock mutex");    
#endif /* HAVE_PTHREAD */
    return *item; 
}
