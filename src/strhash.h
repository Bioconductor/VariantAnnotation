#ifndef _STRHASH_H_
#define _STRHASH_H_

#include "samtools/khash.h"

KHASH_SET_INIT_STR(strhash)

khash_t(strhash) *_strhash_new();
void _strhash_free(khash_t(strhash) *str);
const char *_strhash_put(khash_t(strhash) *str, const char *value);

#endif
