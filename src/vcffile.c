#include <stdlib.h>
#include <string.h>
#include "vcffile.h"

/* iterator to return null-terminated delimited fields */

struct it {
    char *str;
    char delim;
};

char *_it_next(struct it *it)
{
    char *curr = it->str;
    while ('\0' != *it->str && it->delim != *it->str)
        it->str++;
    if ('\0' != *it->str)
        *it->str++ = '\0';
    return curr;
}

char *_it_init(struct it *it, char *str, char delim)
{
    it->str = str;
    it->delim = delim;
    return _it_next(it);
}

/* parse a single 'vcf' line into vectors and matricies 'map' */

const struct fld_fmt {
    const char *name;
    SEXPTYPE type;
} FLD_FMT[] = {
    {
    "CHROM", STRSXP}, {
    "POS", INTSXP}, {
    "ID", STRSXP}, {
    "REF", STRSXP}, {
    "ALT", STRSXP}, {
    "QUAL", REALSXP}, {
    "FILTER", STRSXP}, {
    "INFO", STRSXP}, {
    "GENO", VECSXP}
};

const int N_FLDS = sizeof(FLD_FMT) / sizeof(FLD_FMT[0]);

SEXP
_alloc_types_list(int vcf_n, int col_n, SEXP map, SEXP eltnms)
{
    int i, j, map_n = Rf_length(map);
    SEXP types;

    /* no INFO or GENO information in header */
    if (map_n == 0) {
        PROTECT(types = Rf_allocVector(VECSXP, 1));
        SEXP elt = Rf_allocMatrix(STRSXP, vcf_n, 1);
        SET_VECTOR_ELT(types, 0, elt);
        for (i = 0; i < vcf_n; ++i)
            SET_STRING_ELT(elt, i, R_NaString);
    }
    else {
        PROTECT(types = Rf_allocVector(VECSXP, map_n));

        for (j = 0; j < map_n; ++j) {
            SEXPTYPE type = TYPEOF(VECTOR_ELT(map, j));
            if (NILSXP == type) {
                SET_VECTOR_ELT(types, j, R_NilValue);
                continue;
            }
            SEXP elt = Rf_allocMatrix(type, vcf_n, col_n);
            SET_VECTOR_ELT(types, j, elt);
            switch (type) {
            case LGLSXP:
                for (i = 0; i < vcf_n * col_n; ++i)
                    LOGICAL(elt)[i] = FALSE;
                break;
            case INTSXP:
                for (i = 0; i < vcf_n * col_n; ++i)
                    INTEGER(elt)[i] = R_NaInt;
                break;
            case REALSXP:
                for (i = 0; i < vcf_n * col_n; ++i)
                    REAL(elt)[i] = R_NaReal;
                break;
            case STRSXP:
                for (i = 0; i < vcf_n * col_n; ++i)
                    SET_STRING_ELT(elt, i, R_NaString);
                break;
            default:
                Rf_error("(internal) unhandled type '%s'",
                         type2char(type));
            }
            if (R_NilValue != eltnms)
                elt = Rf_dimnamesgets(elt, eltnms);
        }
    }

    UNPROTECT(1);
    return types;
}

SEXP
_trim_null(SEXP data, SEXP nms)
{
   int i, j = 0;
   for (i = 0; i < Rf_length(data); ++i) {
       if (R_NilValue != VECTOR_ELT(data, i)) {
           SET_VECTOR_ELT(data, j, VECTOR_ELT(data, i));
           SET_STRING_ELT(nms, j, STRING_ELT(nms, i));
           j++;
       }
   }
   PROTECT(nms = Rf_lengthgets(nms, j));
   PROTECT(data = Rf_lengthgets(data, j));
   data = Rf_namesgets(data, nms);
   UNPROTECT(2);

   return data;
}

SEXP _split_vcf(SEXP vcf, SEXP sample, SEXP imap, SEXP gmap)
{
    int i, j;
    const int
        vcf_n = Rf_length(vcf),
        samp_n = Rf_length(sample),
        imap_n = Rf_length(imap),
        gmap_n = Rf_length(gmap);

    SEXP gnms = GET_NAMES(gmap);
    SEXP inms = GET_NAMES(imap);
    int fmtidx, sampleidx;
    int *gmapidx = (int *) R_alloc(sizeof(int), gmap_n), imapidx;

    /* allocate result and first 7 fixed fields */
    SEXP result, info, geno, eltnms;
    PROTECT(result = Rf_allocVector(VECSXP, N_FLDS));
    for (i = 0; i < N_FLDS - 2; ++i)
        SET_VECTOR_ELT(result, i, Rf_allocVector(FLD_FMT[i].type, vcf_n));

    PROTECT(eltnms = Rf_allocVector(STRSXP, N_FLDS));
    for (i = 0; i < N_FLDS; ++i)
        SET_STRING_ELT(eltnms, i, mkChar(FLD_FMT[i].name));
    result = Rf_namesgets(result, eltnms);
    UNPROTECT(1);

    /* allocate info */
    PROTECT(info = _alloc_types_list(vcf_n, 1, imap, R_NilValue));

    /* allocate geno, including matricies */
    PROTECT(eltnms = Rf_allocVector(VECSXP, 2));
    SET_VECTOR_ELT(eltnms, 0, R_NilValue);
    SET_VECTOR_ELT(eltnms, 1, sample);
    PROTECT(geno = _alloc_types_list(vcf_n, samp_n, gmap, eltnms));

    /* parse each line */
    for (i = 0; i < vcf_n; i++) {
        struct it it0, it1, it2;
        char *record, *sample, *field, *ifld, *ikey, *fmt;

        record = strdup(CHAR(STRING_ELT(vcf, i)));

        /* first 7 'fixed' fields */
        for (field = _it_init(&it0, record, '\t'), j = 0;
             j < N_FLDS - 2; field = _it_next(&it0), ++j) {
            SEXP elt = VECTOR_ELT(result, j);
            switch (TYPEOF(elt)) {
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

        /* 'INFO' field */
        int midx = i;
        /* INFO field not parsed if no header information present */
        if (imap_n == 0) {
            SEXP matrix = VECTOR_ELT(info, 0);
            SET_STRING_ELT(matrix, midx, mkChar(field)); 
        } else { 
            for (ifld = _it_init(&it1, field, ';'); '\0' != *ifld;
                ifld = _it_next(&it1)) {
                ikey = _it_init(&it2, ifld, '=');
                for (imapidx = 0; imapidx < imap_n; ++imapidx) {
                    if (0L == strcmp(ikey, CHAR(STRING_ELT(inms, imapidx))))
                        break;
                }
                if (imap_n == imapidx)
                    Rf_error("record %d INFO '%s' not found", i + 1, ikey);

                SEXP matrix = VECTOR_ELT(info, imapidx);
                if (LGLSXP == TYPEOF(matrix)) {
                    LOGICAL(matrix)[midx] = TRUE;
                } else {
                    field = _it_next(&it2);
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
        }

        /* 'FORMAT' field */
        field = _it_next(&it0);
        fmt = field;
        for (field = _it_init(&it2, fmt, ':'), fmtidx = 0;
             '\0' != *field; field = _it_next(&it2), fmtidx++) {
            for (j = 0; j < gmap_n; ++j) {
                if (0L == strcmp(field, CHAR(STRING_ELT(gnms, j))))
                    break;
            }
            if (gmap_n == j)
                Rf_error("record %d field %d FORMAT '%s' not found",
                         i + 1, fmtidx + 1, field);
            gmapidx[fmtidx] = j;
        }

        /* 'samples' field(s) */
        for (sample = _it_next(&it0), sampleidx = 0;
             '\0' != *sample; sample = _it_next(&it0), sampleidx++) {
            for (field = _it_init(&it2, sample, ':'), fmtidx = 0;
                 '\0' != *field; field = _it_next(&it2), fmtidx++) {
                SEXP matrix = VECTOR_ELT(geno, gmapidx[fmtidx]);
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
    if (NILSXP == TYPEOF(inms)) {
        PROTECT(inms = Rf_allocVector(STRSXP, 1));
        SET_STRING_ELT(inms, 0, mkChar("INFO"));
        UNPROTECT(1);
    }
    /* remove NULL elements of info and geno  */
    SET_VECTOR_ELT(result, N_FLDS - 2, _trim_null(info, inms));
    SET_VECTOR_ELT(result, N_FLDS - 1, _trim_null(geno, gnms));

    UNPROTECT(4);
    return result;
}

SEXP scan_vcf(SEXP tbx, SEXP sample, SEXP imap, SEXP gmap)
{
    SEXP result = PROTECT(Rf_allocVector(VECSXP, Rf_length(tbx)));

    for (int i = 0; i < Rf_length(tbx); ++i) {
        SEXP record = _split_vcf(VECTOR_ELT(tbx, i), sample, imap, gmap);
        SET_VECTOR_ELT(result, i, record);
    }

    UNPROTECT(1);
    return result;
}

SEXP scan_vcf_connection(SEXP txt, SEXP sample, SEXP imap, SEXP gmap)
{
    SEXP result = PROTECT(Rf_allocVector(VECSXP, 1));

    SET_VECTOR_ELT(result, 0, _split_vcf(txt, sample, imap, gmap));

    UNPROTECT(1);
    return result;
}
