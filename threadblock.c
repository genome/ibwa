#include "threadblock.h"

#include "utils.h"

#include <stdlib.h>
#ifdef HAVE_PTHREAD
#include <pthread.h>
#endif /* HAVE_PTHREAD */

typedef struct {
    pthread_t tid;
    uint32_t size;
    uint32_t idx;
    void *data;
    threadblockfn_t func;
} __threadblock_t;

static void *__threadblock_dispatch(void *data) {
    __threadblock_t *tb = (__threadblock_t*)data;
    tb->func(tb->idx, tb->size, tb->data);
    return NULL;
}

void threadblock_exec(uint32_t size, threadblockfn_t func, void *data)
{
    uint32_t i;

    __threadblock_t *tb = calloc(size, sizeof(__threadblock_t));
    for (i = 0; i < size; ++i) {
        tb[i].size = size;
        tb[i].idx = i;
        tb[i].data = data;
        tb[i].func = func;

#ifdef HAVE_PTHREAD
        if (pthread_create(&tb[i].tid, NULL, __threadblock_dispatch, tb+i) != 0)
            err_fatal_simple("thread creation failed.");
#else /* HAVE_PTHREAD */
        __threadblock_dispatch(tb+i); /* dispatch inline */
#endif /* HAVE_PTHREAD */
    }

#ifdef HAVE_PTHREAD
    for (i = 0; i < size; ++i)
        pthread_join(tb[i].tid, NULL);
#endif /* HAVE_PTHREAD */
    free(tb);
}
