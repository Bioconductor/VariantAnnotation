#ifndef _VCFTYPE_H_
#define _VCFTYPE_H_

#include <stdlib.h>
#include <string.h>
#include <Rdefines.h>
#include "utilities.h"

struct vcftype_t {
    SEXPTYPE type;
    Rboolean isArray;
    int nrow, ncol;
    union {
        int *logical;
        int *integer;
        double *numeric;
        char **character;
        struct vcftype_t **list;
    } u;
};

struct vcftype_t *_vcftype_new(SEXPTYPE type, int nrow, int ncol,
                               const Rboolean isArray);
void _vcftype_free(struct vcftype_t *vcftype);
struct vcftype_t *_vcftype_grow(struct vcftype_t *vcftype, int nrow);
SEXP _vcftype_as_SEXP(struct vcftype_t *vcftype);

static inline void _vcftype_set(struct vcftype_t *vcftype,
                                const int idx, const char *field)
{
    switch (vcftype->type) {
    case NILSXP:
        break;
    case INTSXP:
        vcftype->u.integer[idx] = atoi(field);
        break;
    case REALSXP:
        vcftype->u.numeric[idx] =
            ('.' == *field) ? R_NaReal : atof(field);
        break;
    case STRSXP:
        vcftype->u.character[idx] = Strdup(field);
        break;
    default:
        Rf_error("(internal) unhandled field type '%s'",
                 type2char(vcftype->type));
    }
}

#endif
