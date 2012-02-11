#ifndef _RLE_H
#define _RLE_H

#include <Rdefines.h>

struct rle_t {
    int len, size, *length;
    char **value;
};

struct rle_t *rle_new(const int size);
void rle_free(struct rle_t *rle);
void rle_grow(struct rle_t *rle, int size);
void rle_append(struct rle_t *rle, const char *value);
SEXP rle_as_Rle(struct rle_t *rle);

#endif
