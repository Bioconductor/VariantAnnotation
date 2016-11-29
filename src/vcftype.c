#include "vcftype.h"

struct vcftype_t *_vcftype_new(SEXPTYPE type, SEXPTYPE listtype,
                               char number, const char *charDotAs,
                               int nrow, int ncol, int ndim, int arrayDim)
{
    struct vcftype_t *vcftype = Calloc(1, struct vcftype_t);
    vcftype->type = type;
    vcftype->listtype = listtype; /* VECSXP: ragged array */
    vcftype->number = number;     /* 'A' or '.' for ragged array only */
    vcftype->charDotAs = charDotAs;
    vcftype->ncol = ncol;
    vcftype->ndim = ndim;
    vcftype->arrayDim = arrayDim;

    return _vcftype_grow(vcftype, nrow);
}

void _vcftype_free(struct vcftype_t *vcftype)
{
    if (NULL == vcftype)
        return;
    int sz = vcftype->nrow * vcftype->ncol * vcftype->ndim;
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
        if (NULL != vcftype->u.character)
            Free(vcftype->u.character);
        break;
    case VECSXP:
        if (NULL != vcftype->u.list) {
            for (int i = 0; i < sz; ++i)
                if (NULL != vcftype->u.list[i])
                    _vcftype_free(vcftype->u.list[i]);
            Free(vcftype->u.list);
        }
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
    int ncol = vcftype->ncol, ndim = vcftype->ndim,
        o_nrow = vcftype->nrow, osz = o_nrow * ncol * ndim,
        sz = nrow * ncol * ndim;

    if (nrow < 0)
        Rf_error("(internal) _vcftype_grow 'nrow' < 0");
    if (sz < 0)
        Rf_error("(internal) _vcftype_grow 'sz' < 0; cannot allocate memory?");

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
        vcftype->u.character = (const char **)
            vcf_Realloc(vcftype->u.character, sz * sizeof(const char *));
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

#define TPOSE(to, from, nrow, ncol, ndim)                          \
    for (int k = 0; k < (ndim); ++k)                               \
        for (int j = 0; j < (ncol); ++j)                           \
            for (int i = 0; i < (nrow); ++i)                       \
                *(to)++ = (from)[i * (ncol) * (ndim) + j * (ndim) + k]

SEXP _vcftype_as_SEXP(struct vcftype_t *vcftype)
{
    if (NULL == vcftype || NILSXP == vcftype->type)
        return R_NilValue;

    const int ncol = vcftype->ncol, ndim = vcftype->ndim,
        nrow = vcftype->nrow, sz = nrow * ncol * ndim;
    SEXP ans = PROTECT(Rf_allocVector(vcftype->type, sz));
    int *ival, idx;
    double *dval;

    switch (vcftype->type) {
    case LGLSXP:
        ival = LOGICAL(ans);
        TPOSE(ival, (const int *) vcftype->u.logical, nrow, ncol, ndim);
        Free(vcftype->u.logical);
        break;
    case INTSXP:
        ival = INTEGER(ans);
        TPOSE(ival, (const int *) vcftype->u.integer, nrow, ncol, ndim);
        Free(vcftype->u.integer);
        break;
    case REALSXP:
        dval = REAL(ans);
        TPOSE(dval, (const double *) vcftype->u.numeric, nrow, ncol, ndim);
        Free(vcftype->u.numeric);
        break;
    case STRSXP:
        idx = 0;
        for (int k = 0; k < ndim; ++k)
            for (int j = 0; j < ncol; ++j)
                for (int i = 0; i < nrow; ++i) {
                    const int idx0 = i * ncol * ndim + j * ndim + k;
                    const char * const s = vcftype->u.character[idx0];
                    const SEXP elt = (NULL == s) ? R_NaString : mkChar(s);
                    SET_STRING_ELT(ans, idx++, elt);
                }
        Free(vcftype->u.character);
        break;
    case VECSXP:
        idx = 0;
        for (int k = 0; k < ndim; ++k)
            for (int j = 0; j < ncol; ++j)
                for (int i = 0; i < nrow; ++i) {
                    const int idx0 = i * ncol * ndim + j * ndim + k;
                    struct vcftype_t *t = vcftype->u.list[idx0];
                    const SEXP elt = (NULL == t) ?
                        Rf_allocVector(vcftype->listtype, 0) :
                        _vcftype_as_SEXP(t);
                    SET_VECTOR_ELT(ans, idx++, elt);
                }
        Free(vcftype->u.list);
        break;
    default:
        Rf_error("(internal) unhandled type '%s'",
                 type2char(vcftype->type));
    }

    if (vcftype->arrayDim > 1) {
        SEXP dim = PROTECT(Rf_allocVector(INTSXP, vcftype->arrayDim));
        INTEGER(dim)[0] = nrow;
        if (vcftype->arrayDim == 2) {
            INTEGER(dim)[1] = ncol * ndim;
        } else {
            INTEGER(dim)[1] = ncol;
            INTEGER(dim)[2] = ndim;
        }
        Rf_setAttrib(ans, R_DimSymbol, dim);
        UNPROTECT(1);
    }

    _vcftype_free(vcftype);
    UNPROTECT(1);
    return ans;
}

void _vcftype_set(struct vcftype_t *vcftype,
                  const int idx, const char *field)
{
    if (NULL == vcftype)
        return;
    switch (vcftype->type) {
    case NILSXP:
        break;
    case LGLSXP:
        vcftype->u.logical[idx] = TRUE;
        break;
    case INTSXP:
        vcftype->u.integer[idx] =
            ('.' == *field) ? R_NaInt : atoi(field);
        break;
    case REALSXP:
        vcftype->u.numeric[idx] =
            ('.' == *field) ? R_NaReal : atof(field);
        break;
    case STRSXP:
        vcftype->u.character[idx] =
            ('.' == *field) ? vcftype->charDotAs : field;
        break;
    default:
        Rf_error("(internal) unhandled field type '%s'",
                 type2char(vcftype->type));
    }
}

void _vcftype_padarray(struct vcftype_t *vcftype,
                       const int irow, const int icol,
                       khash_t(strhash) *str, const int ragged_n)
{
    if (NULL == vcftype)
        return;
    const int offset = irow * vcftype->ncol + icol;
    if (vcftype->u.list[offset] != NULL)
        return;
    _vcftype_setarray(vcftype, irow, icol, "", ragged_n, str);
}

void _vcftype_setarray(struct vcftype_t *vcftype,
                       const int irow, const int icol, char *field,
                       int ragged_n, khash_t(strhash) *str)
{
    struct it_t it;
    char *ifld;

    if (NULL == vcftype)
        return;

    if (VECSXP == vcftype->type) { /* ragged array */
        /* 'G': one value per genotype */
        if (vcftype->number == 'G')
            ragged_n = ((ragged_n + 1) * (ragged_n + 2))/2;
        /* 'R': one value per alternate allele + ref*/
        else if (vcftype->number == 'R')
            ragged_n = ragged_n + 1;
        /* 'A': one value per alternate allele */
        else if (vcftype->number != 'A')
            ragged_n = _vcftype_ragged_n(field);
 
        /* allocate and fill */
        const int offset = irow * vcftype->ncol + icol;
        vcftype->u.list[offset] =
            _vcftype_new(vcftype->listtype, NILSXP, '\0',
                         vcftype->charDotAs, ragged_n, 1, 1, 0);
        ifld = it_init(&it, field, ',');
        for (int k = 0; k < ragged_n; ++k) {
            if ('\0' == *ifld)
                ifld = ".";
            _vcftype_set(vcftype->u.list[offset], k, _strhash_put(str, ifld));
            ifld = it_next(&it);
        }
    } else {                       /* array */
        const int offset = (irow * vcftype->ncol + icol) * vcftype->ndim;
        ifld = it_init(&it, field, ',');
        for (int k = 0; k < vcftype->ndim; ++k) {
            _vcftype_set(vcftype, offset + k, _strhash_put(str, ifld));
            ifld = it_next(&it);
        }
    }
}
