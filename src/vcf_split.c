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

const struct fld_fmt {
    const char *name;
    SEXPTYPE type;
} FLD_FMT[] = {
    {"CHROM", STRSXP}, {"POS", INTSXP}, {"ID", STRSXP},
    {"REF", STRSXP}, {"ALT", STRSXP}, {"QUAL", REALSXP},
    {"FILTER", STRSXP}, {"INFO", STRSXP}, {"GENO", VECSXP}
};

const int N_FLDS = 9;

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

    /* allocate result and fixed fields */
    SEXP result, geno, eltnms;
    PROTECT(result = Rf_allocVector(VECSXP, N_FLDS));
    for (i = 0; i < N_FLDS - 1; ++i)
        SET_VECTOR_ELT(result, i,
                       Rf_allocVector(FLD_FMT[i].type, vcf_n));
    geno = Rf_allocVector(VECSXP, map_n);
    SET_VECTOR_ELT(result, N_FLDS - 1, geno);

    PROTECT(eltnms = Rf_allocVector(STRSXP, N_FLDS));
    for (i = 0; i < N_FLDS; ++i)
        SET_STRING_ELT(eltnms, i, mkChar(FLD_FMT[i].name));
    result = Rf_namesgets(result, eltnms);
    UNPROTECT(1);

    /* allocate geno, including matricies */
    geno = Rf_namesgets(geno, nms);
    PROTECT(eltnms = Rf_allocVector(VECSXP, 2));
    SET_VECTOR_ELT(eltnms, 0, R_NilValue);
    SET_VECTOR_ELT(eltnms, 1, sample);
    for (j = 0; j < map_n; ++j) {
        SEXPTYPE type = TYPEOF(VECTOR_ELT(map, j));
        if (NILSXP == type) {
            SET_VECTOR_ELT(geno, j, R_NilValue);
            continue;
        }
        SEXP elt = Rf_allocMatrix(type, vcf_n, samp_n);
        SET_VECTOR_ELT(geno, j, elt);
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
        elt = Rf_dimnamesgets(elt, eltnms);
    }
    UNPROTECT(1);

    /* parse each line */
    for (i = 0; i < vcf_n; i++) {
        struct it it0, it1;
        char *record, *sample, *fmt, *field;

        record = strdup(CHAR(STRING_ELT(vcf, i)));

        /* 'fixed' fields */
        for (field = _it_init(&it0, record, '\t'), j = 0;
             j < N_FLDS - 1;
             field = _it_next(&it0), ++j)
        {
            SEXP elt = VECTOR_ELT(result, j);
            switch(TYPEOF(elt)) {
            case INTSXP:
                INTEGER(elt)[i] = atoi(field);
                break;
            case REALSXP:
                REAL(elt)[i] = atof(field);
                break;
            case STRSXP:
                SET_STRING_ELT(elt, i, mkChar(field));
                break;
            default:
                Rf_error("(internal) unhandled fixed field type '%s'",
                         type2char(TYPEOF(elt)));
            }
        }

        /* 'FORMAT' field */
        fmt = field;
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

        /* 'samples' field(s) */
        for (sample = _it_next(&it0), sampleidx = 0;
             '\0' != *sample;
             sample = _it_next(&it0), sampleidx++)
        {
            for (field = _it_init(&it1, sample, ':'), fmtidx = 0;
                 '\0' != *field;
                 field = _it_next(&it1), fmtidx++)
            {
                SEXP matrix = VECTOR_ELT(geno, mapidx[fmtidx]);
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

    /* remove NULL elements of geno */
    for (i = 0, j = 0; i < Rf_length(geno); ++i)
        if (R_NilValue != VECTOR_ELT(geno, i))
            SET_VECTOR_ELT(geno, j++, VECTOR_ELT(geno, i));
    geno = Rf_lengthgets(geno, j);

    UNPROTECT(1);
    return result;
}
