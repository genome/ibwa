#ifndef THREADBLOCK_H
#define THREADBLOCK_H

#include <stdint.h>

/* thread block functions look like:
 *
 * void *myfunc(uint32_t index, uint32_t nThreads, void* data)
 *
 * and will be called by 'nThreads' threads, and for the ith thread, index=i
 *
 * It should be noted that the return value of thread block functions is
 * ignored.
 */

typedef void (*threadblockfn_t)(uint32_t, uint32_t, void*);
void threadblock_exec(uint32_t size, threadblockfn_t func, void *data);

#endif /* THREADBLOCK_H */
