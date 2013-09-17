#include <Rdefines.h>
#include "strhash.h"
#include "utilities.h"

khash_t(strhash) *_strhash_new()
{
    return kh_init(strhash);
}

void _strhash_free(khash_t(strhash) *str)
{
    khiter_t key;
    for (key = kh_begin(str); key != kh_end(str); ++key) {
        if (kh_exist(str, key))
            Free(kh_key(str, key));
    }
    kh_destroy(strhash, str);
}

const char *_strhash_put(khash_t(strhash) *str, const char *value)
{
    int ret;
    khiter_t key = kh_get(strhash, str, value);
    if (key == kh_end(str))
        key = kh_put(strhash, str, Strdup(value), &ret);
    return kh_key(str, key);
}

