/** Inline or macro bit twiddling utilities.
 * @file bithacks.h
 */
#pragma once

#include <stdint.h>
#include <unistd.h>
#include <strings.h>

/** Find the lowest power of 2 greater than @a x.
 * Per section 3-2 of Warren, Henry S. Hacker's Delight. Addison-Wesley, 2002. */
static inline size_t ceil_pow2(size_t x)
{
     --x;
     x |= x >> 1;
     x |= x >> 2;
     x |= x >> 4;
     x |= x >> 8;
     x |= x >>16;
     return ++x;
}

/** Reverse the nibbles in @a v.
 * Per section 7-1 of Warren, Henry S. Hacker's Delight. Addison-Wesley, 2002. */
static inline uint64_t reverse_nibbles_64(uint64_t v)
{
     v = (v & 0xf0f0f0f0f0f0f0f0LL)>>4  | (v & 0x0f0f0f0f0f0f0f0fLL)<<4;
     v = (v & 0xff00ff00ff00ff00LL)>>8  | (v & 0x00ff00ff00ff00ffLL)<<8;
     v = (v & 0xffff0000ffff0000LL)>>16 | (v & 0x0000ffff0000ffffLL)<<16;
     v = (v & 0xffffffff00000000LL)>>32 | (v & 0x00000000ffffffffLL)<<32;
     return v;
}

/** Find the next bit set in a bitmap, toggling it off in the bitmap
 * before returning its index.
 * @param[inout] p bitmap
 */
static inline uint32_t bitmap_next(uint32_t *p) {
     uint32_t t = *p & -*p;
     *p ^= t;
     return ffs(t)-1;
}
