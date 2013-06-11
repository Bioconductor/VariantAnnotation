#ifndef _VCFTYPE_H_
#define _VCFTYPE_H_

#include <stdlib.h>
#include <string.h>
#include <Rdefines.h>
#include "utilities.h"

struct vcftype_t {
    SEXPTYPE type, listtype;    /* listtype for ragged arrays */
    char number,                /* 'A' or '.' for ragged array only */
        *charDotAs;              /* '.' enocding when in character vectors */
    int nrow, ncol, ndim, arrayDim;
    union {
        int *logical;
        int *integer;
        double *numeric;
        char **character;
        struct vcftype_t **list;
    } u;
};

struct vcftype_t *_vcftype_new(SEXPTYPE type, SEXPTYPE listtype,
                               char number, char *charDotAs,
                               int nrow, int ncol, int ndim,
                               int arrayDim);
void _vcftype_free(struct vcftype_t *vcftype);
struct vcftype_t *_vcftype_grow(struct vcftype_t *vcftype, int nrow);
SEXP _vcftype_as_SEXP(struct vcftype_t *vcftype);
void _vcftype_set(struct vcftype_t *vcftype,
                  const int idx, const char *field);
void _vcftype_setarray(struct vcftype_t *vcftype,
                       const int irow, const int icol, char *field,
                       int ragged_n);
void _vcftype_padarray(struct vcftype_t *vcftype,
                       const int irow, const int icol, const int ragged_n);

static inline int _vcftype_ragged_n(const char *a)
{
    int n = (*a == '\0') ? 0 : 1;
    while (*a != '\0')
        if (*a++ == ',') ++n;
    return n;
}

#endif
