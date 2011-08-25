#include <stdlib.h>
#include <string.h>

#include "vcf_split.h"

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
        SEXP elt = Rf_allocMatrix(TYPEOF(VECTOR_ELT(map, j)),
                                  vcf_n, samp_n);
        SET_VECTOR_ELT(result, j, elt);
        elt = Rf_dimnamesgets(elt, dimnames);
    }
    UNPROTECT(1);

    /* parse each line of vcf */
    for (i = 0; i < vcf_n; i++) {
        const char *s0 = CHAR(STRING_ELT(vcf, i));
        char *record = strdup(s0);
        char *sample, *savesample, *fmt;
        char *field, *savefield;

        fmt = strtok_r(record, "\t", &savesample);

        /* first tab-delimited string is 'FORMAT' */
        field = strtok_r(fmt, ":", &savefield);
        fmtidx = 0;
        while (NULL != field) {
            for (j = 0; j < map_n; ++j) {
                const char *id = CHAR(STRING_ELT(nms, j));
                if (0L == strcmp(field, id))
                    break;
            }
            if (map_n == j) {
                UNPROTECT(1);
                Rf_error("record %d field %d FORMAT '%s' not found",
                         i + 1, fmtidx + 1, field);
            }
            mapidx[fmtidx++] = j;
            field = strtok_r(NULL, ":", &savefield);
        }

        /* process samples */
        sample = strtok_r(NULL, "\t", &savesample);
        sampleidx = 0;
        while (NULL != sample) {
            field = strtok_r(sample, ":", &savefield);
            fmtidx = 0;
            while (NULL != field) {
                SEXP matrix = VECTOR_ELT(result, mapidx[fmtidx]);
                int midx = sampleidx * vcf_n + i;
                switch (TYPEOF(matrix)) {
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
                    Rf_error("(internal) unhandled SEXPTYPE");
                }
                fmtidx++;
                field = strtok_r(NULL, ":", &savefield);
            }
            sampleidx++;
            sample = strtok_r(NULL, "\t", &savesample);
        }
        free(record);
    }

    UNPROTECT(1);
    return result;
}
