#include "rle.h"
#include "utilities.h"

/* rle representation */

static const double RLE_GROW = 1.6;

struct rle_t *rle_new(const int size)
{
    struct rle_t *rle = Calloc(1, struct rle_t);
    rle->size = size;
    rle->len = 0;
    rle->rle = Rf_allocVector(VECSXP, 2);
    R_PreserveObject(rle->rle);
    SET_VECTOR_ELT(rle->rle, 0, Rf_allocVector(STRSXP, size));
    SET_VECTOR_ELT(rle->rle, 1, Rf_allocVector(INTSXP, size));
    rle->value = STRING_PTR(VECTOR_ELT(rle->rle, 0));
    rle->length = INTEGER(VECTOR_ELT(rle->rle, 1));
    return rle;
}

void rle_free(struct rle_t *rle)
{
    R_ReleaseObject(rle->rle);
    Free(rle);
}

void rle_grow(struct rle_t *rle, int size)
{
    SEXP updt;

    rle->size = size;

    updt = Rf_lengthgets(VECTOR_ELT(rle->rle, 0), size);
    SET_VECTOR_ELT(rle->rle, 0, updt) ;
    rle->value = STRING_PTR(updt);

    updt = Rf_lengthgets(VECTOR_ELT(rle->rle, 1), size);
    SET_VECTOR_ELT(rle->rle, 1, updt) ;
    rle->length = INTEGER(updt);
}

void rle_append(struct rle_t *rle, const char *value)
{
    SEXP val = PROTECT(mkChar(value));
    if (rle->len == 0 || val != rle->value[rle->len - 1])
    {
        if (rle->len == rle->size)
            rle_grow(rle, RLE_GROW * rle->size);
        rle->value[rle->len] = val;
        rle->length[rle->len] = 1;
        rle->len += 1;
    } else
        rle->length[rle->len - 1] += 1;
    UNPROTECT(1);
}

SEXP rle_as_Rle(struct rle_t *rle)
{
    SEXP sxp, nmspc, fun, expr;
    if (rle->len != Rf_length(VECTOR_ELT(rle->rle, 0))) {
        SEXP elt;
        PROTECT(sxp = Rf_allocVector(VECSXP, 2));

        elt = Rf_lengthgets(VECTOR_ELT(rle->rle, 0), rle->len);
        SET_VECTOR_ELT(sxp, 0, elt);

        elt = Rf_lengthgets(VECTOR_ELT(rle->rle, 1), rle->len);
        SET_VECTOR_ELT(sxp, 1, elt);
    } else
        PROTECT(sxp = rle->rle); /* not necessary; simpler UNPROTECT... */

    PROTECT(nmspc = get_namespace("IRanges"));
    PROTECT(fun = findFun(install("Rle"), nmspc));
    PROTECT(expr = lang3(fun, VECTOR_ELT(sxp, 0), VECTOR_ELT(sxp, 1)));
    sxp = eval(expr, R_GlobalEnv);
    UNPROTECT(4);
    return sxp;
}
