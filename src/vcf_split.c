#include <stdlib.h>
#include <string.h>

#include "vcf_split.h"

/* iterator to return delimited fields as null-terminated */

struct it {
    char *str;
    char delim;
};

char *
_it_next(struct it *it)
{
    char *curr= it->str;
    while ('\0' != *it->str && it->delim != *it->str)
        it->str++;
    if ('\0' != *it->str)
        *it->str++ = '\0';
    return curr;
}

char *
_it_init(struct it *it, char *str, char delim)
{
    it->str = str;
    it->delim = delim;
    return _it_next(it);
}

/* split 'vcf' into list of matricies based on 'map' */

SEXP
vcf_split(SEXP vcf, SEXP sample, SEXP map)
{
    int i, j;
    const int
        vcf_n = Rf_length(vcf),
        samp_n = Rf_length(sample),
        map_n = Rf_length(map);

    SEXP nms = GET_NAMES(map);
    int fmtidx, sampleidx,
        *mapidx = (int *) R_alloc(sizeof(int), map_n);

    /* allocate result, set names on list, dimnames on matricies */
    SEXP result, dimnames;
    PROTECT(result = Rf_allocVector(VECSXP, map_n));
    result = Rf_namesgets(result, nms);
    PROTECT(dimnames = Rf_allocVector(VECSXP, 2));
    SET_VECTOR_ELT(dimnames, 0, R_NilValue);
    SET_VECTOR_ELT(dimnames, 1, sample);
    for (j = 0; j < map_n; ++j) {
        SEXPTYPE type = TYPEOF(VECTOR_ELT(map, j));
        if (NILSXP == type) {
            SET_VECTOR_ELT(result, j, R_NilValue);
            continue;
        }
        SEXP elt = Rf_allocMatrix(type, vcf_n, samp_n);
        SET_VECTOR_ELT(result, j, elt);
        switch (type) {
        case INTSXP:
            for (i = 0; i < vcf_n * samp_n; ++i)
                INTEGER(elt)[i] = R_NaInt;
            break;
        case REALSXP:
            for (i = 0; i < vcf_n * samp_n; ++i)
                REAL(elt)[i] = R_NaReal;
            break;
        case STRSXP:
            for (i = 0; i < vcf_n * samp_n; ++i)
                SET_STRING_ELT(elt, i, R_NaString);
            break;
        default:
            Rf_error("(internal) unhandled type '%s'", type2char(type));
        }
        elt = Rf_dimnamesgets(elt, dimnames);
    }
    UNPROTECT(1);

    /* parse each line of vcf */
    for (i = 0; i < vcf_n; i++) {
        struct it it0, it1;
        char *record, *sample, *fmt, *field;

        record = strdup(CHAR(STRING_ELT(vcf, i)));

        /* first tab-delimited string is 'FORMAT' */
        fmt = _it_init(&it0, record, '\t');
        for (field = _it_init(&it1, fmt, ':'), fmtidx = 0;
             '\0' != *field;
             field = _it_next(&it1), fmtidx++)
        {
            for (j = 0; j < map_n; ++j) {
                if (0L == strcmp(field, CHAR(STRING_ELT(nms, j))))
                    break;
            }
            if (map_n == j)
                Rf_error("record %d field %d FORMAT '%s' not found",
                         i + 1, fmtidx + 1, field);
            mapidx[fmtidx] = j;
        }

        /* process samples */
        for (sample = _it_next(&it0), sampleidx = 0;
             '\0' != *sample;
             sample = _it_next(&it0), sampleidx++)
        {
            for (field = _it_init(&it1, sample, ':'), fmtidx = 0;
                 '\0' != *field;
                 field = _it_next(&it1), fmtidx++)
            {
                SEXP matrix = VECTOR_ELT(result, mapidx[fmtidx]);
                int midx = sampleidx * vcf_n + i;
                switch (TYPEOF(matrix)) {
                case NILSXP:
                    break;
                case INTSXP:
                    INTEGER(matrix)[midx] = atoi(field);
                    break;
                case REALSXP:
                    REAL(matrix)[midx] = atof(field);
                    break;
                case STRSXP:
                    SET_STRING_ELT(matrix, midx, mkChar(field));
                    break;
                default:
                    Rf_error("(internal) unhandled type '%s'",
                             type2char(TYPEOF(matrix)));
                }
            }
        }
        free(record);
    }

    /* remove NULL elements of result */
    for (i = 0, j = 0; i < Rf_length(result); ++i)
        if (R_NilValue != VECTOR_ELT(result, i))
            SET_VECTOR_ELT(result, j++, VECTOR_ELT(result, i));
    result = Rf_lengthgets(result, j);

    UNPROTECT(1);
    return result;
}
