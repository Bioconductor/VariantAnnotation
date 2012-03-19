#include "rle.h"
#include "utilities.h"

/* rle representation */

static const double RLE_GROW = 1.6;

struct rle_t *rle_new(const int size)
{
    struct rle_t *rle = Calloc(1, struct rle_t);
    rle->size = size;
    rle->len = 0;
    rle->value = Calloc(size, char *);
    rle->length = Calloc(size, int);
    return rle;
}

void rle_free(struct rle_t *rle)
{
    for (int i = 0; i < rle->len; ++i)
        Free(rle->value[i]);
    Free(rle->value);
    Free(rle->length);
    Free(rle);
}

void rle_grow(struct rle_t *rle, int size)
{
    rle->value = Realloc(rle->value, size, char *);
    rle->length = Realloc(rle->length, size, int);
    rle->size = size;
}

void rle_append(struct rle_t *rle, const char *value)
{
    if (0 == rle->len || 0 != strcmp(value, rle->value[rle->len - 1]))
    {
        if (rle->len == rle->size)
            rle_grow(rle, RLE_GROW * rle->size);
        rle->value[rle->len] = Strdup(value);
        rle->length[rle->len] = 1;
        rle->len += 1;
    } else
        rle->length[rle->len - 1] += 1;
}

SEXP rle_as_Rle(struct rle_t *rle)
{
    SEXP value, length, nmspc, fun, expr, result;

    PROTECT(value = Rf_allocVector(STRSXP, rle->len));
    PROTECT(length = Rf_allocVector(INTSXP, rle->len));
    for (int i = 0; i < rle->len; ++i) {
        SET_STRING_ELT(value, i, mkChar(rle->value[i]));
        INTEGER(length)[i] = rle->length[i];
    }

    PROTECT(nmspc = get_namespace("IRanges"));
    PROTECT(fun = findFun(install("Rle"), nmspc));
    PROTECT(expr = lang3(fun, value, length));
    result = eval(expr, R_GlobalEnv);
    UNPROTECT(5);
    return result;
}
