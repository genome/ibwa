#ifndef __byteorder_h__
#define __byteorder_h__

#include <stdint.h>

int is_big_endian();

uint16_t swap_endian_2(uint16_t v);
void *swap_endian_2p(void *x);
uint32_t swap_endian_4(uint32_t v);
void *swap_endian_4p(void *x);
uint64_t swap_endian_8(uint64_t v);
void *swap_endian_8p(void *x);

#endif /* __byteorder_h__ */
