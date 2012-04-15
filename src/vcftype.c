#include "vcftype.h"

struct vcftype_t *_vcftype_new(SEXPTYPE type, int nrow, int ncol,
                               const Rboolean isArray)
{
    struct vcftype_t *vcftype = Calloc(1, struct vcftype_t);
    vcftype->isArray = isArray;
    vcftype->ncol = ncol;
    vcftype->type = type;

    return _vcftype_grow(vcftype, nrow);
}

void _vcftype_free(struct vcftype_t *vcftype)
{
    if (NULL == vcftype)
        return;
    int sz = vcftype->nrow * (0 == vcftype->ncol ? 1 : vcftype->ncol);
    switch (vcftype->type) {
    case NILSXP:
        break;
    case LGLSXP:
        Free(vcftype->u.logical);
        break;
    case INTSXP:
        Free(vcftype->u.integer);
        break;
    case REALSXP:
        Free(vcftype->u.numeric);
        break;
    case STRSXP:
        for (int i = 0; i < sz; ++i)
            Free(vcftype->u.character[i]);
        Free(vcftype->u.character);
        break;
    case VECSXP:
        for (int i = 0; i < sz; ++i)
            _vcftype_free(vcftype->u.list[i]);
        Free(vcftype->u.list);
        break;
    default:
        Rf_error("(internal) unhandled type '%s'",
                 type2char(vcftype->type));
    }
    Free(vcftype);
}

void *vcf_Realloc(void * p, size_t n)
{
    /* Realloc(p, 0, *) fails inappropriately */
    if (n == 0) {
        Free(p);
        p = NULL;
    } else {
        p = R_chk_realloc(p, n);
    }
    return p;
}

struct vcftype_t *_vcftype_grow(struct vcftype_t * vcftype, int nrow)
{
    if (NULL == vcftype)
        return vcftype;
    int ncol = (0 == vcftype->ncol) ? 1 : vcftype->ncol;
    int osz = vcftype->nrow * ncol, sz = nrow * ncol;
    switch (vcftype->type) {
    case NILSXP:
        break;
    case LGLSXP:
        vcftype->u.logical = (int *)
            vcf_Realloc(vcftype->u.logical, sz * sizeof(int));
        for (int i = osz; i < sz; ++i)
            vcftype->u.logical[i] = FALSE;
        break;
    case INTSXP:
        vcftype->u.integer = (int *)
            vcf_Realloc(vcftype->u.integer, sz * sizeof(int));
        for (int i = osz; i < sz; ++i)
            vcftype->u.integer[i] = R_NaInt;
        break;
    case REALSXP:
        vcftype->u.numeric = (double *)
            vcf_Realloc(vcftype->u.numeric, sz * sizeof(double));
        for (int i = osz; i < sz; ++i)
            vcftype->u.numeric[i] = R_NaReal;
        break;
    case STRSXP:
        vcftype->u.character = (char **)
            vcf_Realloc(vcftype->u.character, sz * sizeof(char *));
        for (int i = osz; i < sz; ++i)
            vcftype->u.character[i] = NULL;
        break;
    case VECSXP:
        vcftype->u.list = (struct vcftype_t **)
            vcf_Realloc(vcftype->u.list, sz * sizeof(struct vcftype_t *));
        for (int i = osz; i < sz; ++i)
            vcftype->u.list[i] = NULL;
        break;
    default:
        Rf_error("(internal) unhandled type '%s'",
                 type2char(vcftype->type));
    }
    vcftype->nrow = nrow;

    return vcftype;
}

#define TPOSE(to, from, nrow, ncol)                                \
    for (int j = 0; j < (ncol); ++j)                               \
        for (int i = 0; i < (nrow); ++i)                           \
            *(to)++ = (from)[i * (ncol) + j]

SEXP _vcftype_as_SEXP(struct vcftype_t *vcftype)
{
    if (NULL == vcftype || NILSXP == vcftype->type)
        return R_NilValue;

    const int ncol = vcftype->isArray ? vcftype->ncol : 1,
        nrow = vcftype->nrow;
    SEXP ans = PROTECT(Rf_allocVector(vcftype->type, nrow * ncol));

    /* FIXME: transpose matricies */
    switch (vcftype->type) {
    case LGLSXP: {
        int *val = LOGICAL(ans);
        TPOSE(val, vcftype->u.logical, nrow, ncol);
    }
        break;
    case INTSXP: {
        int *val = INTEGER(ans);
        TPOSE(val, vcftype->u.integer, nrow, ncol);
    }
        break;
    case REALSXP: {
        double *val = REAL(ans);
        TPOSE(val, vcftype->u.numeric, nrow, ncol);
    }
        break;
    case STRSXP: {
        SEXP elt;
        int idx = 0;
        for (int j = 0; j < ncol; ++j)
            for (int i = 0; i < nrow; ++i) {
                const char *s = vcftype->u.character[i * ncol + j];
                elt = (NULL == s) ? R_NaString : mkChar(s);
                SET_STRING_ELT(ans, idx++, elt);
                Free(vcftype->u.character[i * ncol + j]);
                vcftype->u.character[i * ncol + j] = NULL;
            }
    }
        break;
    case VECSXP: {
        SEXP elt;
        int idx = 0;
        for (int j = 0; j < ncol; ++j)
            for (int i = 0; i < nrow; ++i) {
                struct vcftype_t *t = vcftype->u.list[i * ncol + j];
                elt = (NULL == t) ? R_NilValue : _vcftype_as_SEXP(t);
                SET_VECTOR_ELT(ans, idx++, elt);
                vcftype->u.list[i * ncol + j] = NULL;
        }
    }
        break;
    default:
        Rf_error("(internal) unhandled type '%s'",
                 type2char(vcftype->type));
    }
    if (TRUE == vcftype->isArray) {
        SEXP dim = PROTECT(Rf_allocVector(INTSXP, 2));
        INTEGER(dim)[0] = vcftype->nrow;
        INTEGER(dim)[1] = vcftype->ncol;
        Rf_setAttrib(ans, R_DimSymbol, dim);
        UNPROTECT(1);
    }

    _vcftype_free(vcftype);
    UNPROTECT(1);
    return ans;
}
