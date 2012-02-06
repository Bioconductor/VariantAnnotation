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
    {"CHROM", STRSXP},
    {"POS", INTSXP},
    {"ID", STRSXP},
    {"REF", STRSXP},
    {"ALT", STRSXP},
    {"QUAL", REALSXP},
    {"FILTER", STRSXP},
    {"INFO", STRSXP},
    {"GENO", VECSXP}
};

const int N_FLDS = sizeof(FLD_FMT) / sizeof(FLD_FMT[0]);

SEXP _types_list_alloc(int vcf_n, int col_n, SEXP map, SEXP eltnms)
{
    int i, j, map_n = Rf_length(map);
    SEXP types;

    /* case of no INFO or GENO information in header */
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

void _types_list_grow(SEXP types, int vcf_n, int col_n)
{
    SEXP elt, elt_n, dim_n;
    SEXPTYPE type;
    int curr_n, new_n = vcf_n * col_n, i;

    for (int j = 0; j < Rf_length(types); ++j) {
        elt = VECTOR_ELT(types, j);
        if (R_NilValue == elt)
            continue;

        curr_n = Rf_length(elt);
        PROTECT(elt_n = Rf_lengthgets(elt, new_n));
        PROTECT(dim_n = Rf_allocVector(INTSXP, 2));
        INTEGER(dim_n)[0] = vcf_n;
        INTEGER(dim_n)[1] = col_n;
        SET_DIM(elt_n, dim_n);
        SET_DIMNAMES(elt_n, GET_DIMNAMES(elt));
        SET_VECTOR_ELT(types, j, elt_n);
        UNPROTECT(2);
        type = TYPEOF(elt_n);
        switch (type) {
            case LGLSXP:
                for (i = curr_n + 1; i < new_n; ++i)
                    LOGICAL(elt_n)[i] = FALSE;
                break;
            case INTSXP:
                for (i = curr_n + 1; i < new_n; ++i)
                    INTEGER(elt_n)[i] = R_NaInt;
                break;
            case REALSXP:
                for (i = curr_n + 1; i < new_n; ++i)
                    REAL(elt_n)[i] = R_NaReal;
                break;
            case STRSXP:
                for (i = curr_n + 1; i < new_n; ++i)
                    SET_STRING_ELT(elt_n, i, R_NaString);
                break;
            default:
                Rf_error("(internal) unhandled type '%s'",
                         type2char(type));
        }
    }
}

void _types_list_transpose(SEXP types)
{
    SEXP elt, dims, orig_dimnames, dimnames, tpose;
    int nrow, ncol, len, l_1, i, j;
    for (int t = 0; t < Rf_length(types); ++t) {
        elt = VECTOR_ELT(types, t);
        if (R_NilValue == elt)
            continue;

        dims = GET_DIM(elt);
        nrow = INTEGER(dims)[0];
        ncol = INTEGER(dims)[1];
        len = Rf_length(elt);
        l_1 = len - 1;

        PROTECT(tpose = Rf_allocMatrix(TYPEOF(elt), ncol, nrow));
        switch (TYPEOF(elt)) {
        case LGLSXP:
            for (i = 0, j = 0; i < len; i++, j += nrow) {
                if (j > l_1) j -= l_1;
                LOGICAL(tpose)[i] = LOGICAL(elt)[j];
            }
            break;
        case INTSXP:
            for (i = 0, j = 0; i < len; i++, j += nrow) {
                if (j > l_1) j -= l_1;
                INTEGER(tpose)[i] = INTEGER(elt)[j];
            }
            break;
        case REALSXP:
            for (i = 0, j = 0; i < len; i++, j += nrow) {
                if (j > l_1) j -= l_1;
                REAL(tpose)[i] = REAL(elt)[j];
            }
            break;
        case STRSXP:
            for (i = 0, j = 0; i < len; i++, j += nrow) {
                if (j > l_1) j -= l_1;
                SET_STRING_ELT(tpose, i, STRING_ELT(elt, j));
            }
            break;
        default:
            Rf_error("(internal) 'transpose' unhandled type '%s'",
                     type2char(TYPEOF(elt)));
        }

        orig_dimnames = GET_DIMNAMES(elt);
        PROTECT(dimnames = Rf_allocVector(VECSXP, 2));
        SET_VECTOR_ELT(dimnames, 0, VECTOR_ELT(orig_dimnames, 1));
        SET_VECTOR_ELT(dimnames, 1, VECTOR_ELT(orig_dimnames, 0));
        SET_DIMNAMES(tpose, dimnames);
        UNPROTECT(1);

        SET_VECTOR_ELT(types, t, tpose);
        UNPROTECT(1);
    }
}

SEXP _trim_null(SEXP data, SEXP nms)
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

SEXP _split_vcf(SEXP vcf, SEXP sample, SEXP fmap, SEXP imap, SEXP gmap)
{
    int i, j;
    const int
        vcf_n = Rf_length(vcf),
        samp_n = Rf_length(sample),
        imap_n = Rf_length(imap),
        gmap_n = Rf_length(gmap);

    SEXP inms = GET_NAMES(imap);
    SEXP gnms = GET_NAMES(gmap);
    int fmtidx, sampleidx;
    int *gmapidx = (int *) R_alloc(sizeof(int), gmap_n), imapidx;

    if (Rf_length(fmap) != N_FLDS - 2)
        Rf_error("[internal] 'fixed' field length %d does not equal %d",
                 Rf_length(fmap), N_FLDS - 2);

    /* allocate result and first 7 fixed fields */
    SEXP result, info, geno, eltnms;
    PROTECT(result = Rf_allocVector(VECSXP, N_FLDS));
    for (i = 0; i < N_FLDS - 2; ++i) {
        SEXP elt = R_NilValue;
        if (R_NilValue != VECTOR_ELT(fmap, i))
            elt = Rf_allocVector(FLD_FMT[i].type, vcf_n);
        SET_VECTOR_ELT(result, i, elt);
    }

    PROTECT(eltnms = Rf_allocVector(STRSXP, N_FLDS));
    for (i = 0; i < N_FLDS; ++i)
        SET_STRING_ELT(eltnms, i, mkChar(FLD_FMT[i].name));
    result = Rf_namesgets(result, eltnms);
    UNPROTECT(1);

    /* allocate info */
    info = _types_list_alloc(vcf_n, 1, imap, R_NilValue);
    SET_VECTOR_ELT(result, N_FLDS - 2, info);

    /* allocate geno, including matricies */
    PROTECT(eltnms = Rf_allocVector(VECSXP, 2));
    SET_VECTOR_ELT(eltnms, 0, R_NilValue);
    SET_VECTOR_ELT(eltnms, 1, sample);
    geno = _types_list_alloc(vcf_n, samp_n, gmap, eltnms);
    SET_VECTOR_ELT(result, N_FLDS - 1, geno);
    UNPROTECT(1);

    /* parse each line */
    for (i = 0; i < vcf_n; i++) {
        struct it it0, it1, it2;
        char *record, *sample, *field, *ifld, *ikey, *fmt;
        char *dot = ".";

        record = strdup(CHAR(STRING_ELT(vcf, i)));

        /* first 7 'fixed' fields */
        for (field = _it_init(&it0, record, '\t'), j = 0;
             j < N_FLDS - 2; field = _it_next(&it0), ++j) {
            SEXP elt = VECTOR_ELT(result, j);
            switch (TYPEOF(elt)) {
            case NILSXP:
                break;
            case INTSXP:
                INTEGER(elt)[i] = atoi(field);
                break;
            case REALSXP:
                if (strcmp(field, dot) == 0)
                    REAL(elt)[i] = R_NaReal;
                else
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

    UNPROTECT(1);
    return result;
}

SEXP scan_vcf_connection(SEXP txt, SEXP sample, SEXP fmap, SEXP imap,
                         SEXP gmap)
{
    SEXP result = PROTECT(Rf_allocVector(VECSXP, 1));

    SET_VECTOR_ELT(result, 0, _split_vcf(txt, sample, fmap, imap, gmap));

    UNPROTECT(1);
    return result;
}

SEXP _vcf_allocate(const int vcf_n, const int samp_n,
                   SEXP sample, SEXP fmap, SEXP imap, SEXP gmap)
{
    /* allocate result and first 7 fixed fields */
    SEXP vcf, info, geno, elt, eltnms;

    if (Rf_length(fmap) != N_FLDS - 2)
        Rf_error("[internal] 'fixed' field length %d does not equal %d",
                 Rf_length(fmap), N_FLDS - 2);

    PROTECT(vcf = Rf_allocVector(VECSXP, N_FLDS));
    PROTECT(eltnms = Rf_allocVector(STRSXP, N_FLDS));
    for (int i = 0; i < N_FLDS; ++i)
        SET_STRING_ELT(eltnms, i, mkChar(FLD_FMT[i].name));
    vcf = Rf_namesgets(vcf, eltnms);
    UNPROTECT(1);

    /* allocate fixed fields */
    for (int i = 0; i < N_FLDS - 2; ++i) {
        elt = R_NilValue;
        if (R_NilValue != VECTOR_ELT(fmap, i))
            elt = Rf_allocVector(FLD_FMT[i].type, vcf_n);
        SET_VECTOR_ELT(vcf, i, elt);
    }

    /* allocate info */
    info = _types_list_alloc(vcf_n, 1, imap, R_NilValue);
    SET_VECTOR_ELT(vcf, N_FLDS - 2, info);

    /* allocate _transposed_ GENO */
    PROTECT(eltnms = Rf_allocVector(VECSXP, 2));
    SET_VECTOR_ELT(eltnms, 0, sample);
    SET_VECTOR_ELT(eltnms, 1, R_NilValue);
    geno = _types_list_alloc(samp_n, vcf_n, gmap, eltnms);
    SET_VECTOR_ELT(vcf, N_FLDS - 1, geno);
    UNPROTECT(1);

    UNPROTECT(1);
    return vcf;
}

void _vcf_grow(SEXP vcf, int vcf_n, int samp_n)
{
    SEXP elt;

    for (int i = 0; i < N_FLDS - 2; ++i) {
        elt = VECTOR_ELT(vcf, i);
        if (R_NilValue != elt)
            SET_VECTOR_ELT(vcf, i, Rf_lengthgets(elt, vcf_n));
    }

    elt = VECTOR_ELT(vcf, N_FLDS - 2);
    _types_list_grow(elt, vcf_n, 1);

    /* _transposed_ matrix */
    elt = VECTOR_ELT(vcf, N_FLDS - 1);
    _types_list_grow(elt, samp_n, vcf_n);
}

SEXP tabix_as_vcf(tabix_t *tabix, ti_iter_t iter, const int size,
                  const Rboolean grow, SEXP state)
{
    SEXP sample = VECTOR_ELT(state, 0), fmap = VECTOR_ELT(state, 1),
        imap = VECTOR_ELT(state, 2), gmap = VECTOR_ELT(state, 3);
    const int col_n = Rf_length(sample), imap_n = Rf_length(imap),
        gmap_n = Rf_length(gmap);
    int vcf_n = size;

    SEXP inms = GET_NAMES(imap);
    SEXP gnms = GET_NAMES(gmap);
    int fmtidx, sampleidx, imapidx,
        *gmapidx = (int *) R_alloc(sizeof(int), gmap_n);

    const double SCALE = 1.6;
    int buflen = 4096;
    char *buf = Calloc(buflen, char);

    int linelen;
    const char *line;

    int irec = 0, j;
    SEXP vcf = _vcf_allocate(vcf_n, col_n, sample, fmap, imap, gmap),
        info, geno, elt;
    PROTECT(vcf);
    info = VECTOR_ELT(vcf, N_FLDS - 2);
    geno = VECTOR_ELT(vcf, N_FLDS - 1);

    while (NULL != (line = ti_read(tabix, iter, &linelen))) {

        struct it it0, it1, it2;
        char *sample, *field, *ifld, *ikey, *fmt;
        char *dot = ".";

        if (irec == vcf_n) {
            if (!grow)
                break;
            vcf_n *= SCALE;
            _vcf_grow(vcf, vcf_n, col_n);
            info = VECTOR_ELT(vcf, N_FLDS - 2);
            geno = VECTOR_ELT(vcf, N_FLDS - 1);
        }

        if (linelen + 1 > buflen) {
            Free(buf);
            buflen = 2 * linelen;
            buf = Calloc(buflen, char);
        }
        memcpy(buf, line, linelen);
        buf[linelen] = '\0';

        /* first 7 'fixed' fields */
        for (field = _it_init(&it0, buf, '\t'), j = 0;
             j < N_FLDS - 2; field = _it_next(&it0), ++j) {
            elt = VECTOR_ELT(vcf, j);
            switch (TYPEOF(elt)) {
            case NILSXP:
                break;
            case INTSXP:
                INTEGER(elt)[irec] = atoi(field);
                break;
            case REALSXP:
                if (strcmp(field, dot) == 0)
                    REAL(elt)[irec] = R_NaReal;
                else
                    REAL(elt)[irec] = atof(field);
                break;
            case STRSXP:
                SET_STRING_ELT(elt, irec, mkChar(field));
                break;
            default:
                Rf_error("(internal) unhandled fixed field type '%s'",
                         type2char(TYPEOF(elt)));
            }
        }

        /* 'INFO' field */
        int midx = irec;
        /* INFO field not parsed if no header information present */
        if (imap_n == 0) {
            elt = VECTOR_ELT(info, 0);
            SET_STRING_ELT(elt, midx, mkChar(field));
        } else {
            for (ifld = _it_init(&it1, field, ';'); '\0' != *ifld;
                 ifld = _it_next(&it1)) {
                ikey = _it_init(&it2, ifld, '=');
                for (imapidx = 0; imapidx < imap_n; ++imapidx) {
                    if (0L == strcmp(ikey, CHAR(STRING_ELT(inms, imapidx))))
                        break;
                }
                if (imap_n == imapidx)
                    Rf_error("record %d INFO '%s' not found",
                             irec + 1, ikey);

                elt = VECTOR_ELT(info, imapidx);
                if (LGLSXP == TYPEOF(elt)) {
                    LOGICAL(elt)[midx] = TRUE;
                } else {
                    field = _it_next(&it2);
                    switch (TYPEOF(elt)) {
                    case NILSXP:
                        break;
                    case INTSXP:
                        INTEGER(elt)[midx] = atoi(field);
                        break;
                    case REALSXP:
                        REAL(elt)[midx] = atof(field);
                        break;
                    case STRSXP:
                        SET_STRING_ELT(elt, midx, mkChar(field));
                        break;
                    default:
                        Rf_error("(internal) unhandled type '%s'",
                                 type2char(TYPEOF(elt)));
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
                         irec + 1, fmtidx + 1, field);
            gmapidx[fmtidx] = j;
        }

        /* 'samples' field(s) */
        for (sample = _it_next(&it0), sampleidx = 0;
             '\0' != *sample; sample = _it_next(&it0), sampleidx++) {
            for (field = _it_init(&it2, sample, ':'), fmtidx = 0;
                 '\0' != *field; field = _it_next(&it2), fmtidx++) {
                elt = VECTOR_ELT(geno, gmapidx[fmtidx]);
                /* GENO matrix entries are _transposed_ for easy grow */
                int midx = irec * col_n + sampleidx;
                switch (TYPEOF(elt)) {
                case NILSXP:
                    break;
                case INTSXP:
                    INTEGER(elt)[midx] = atoi(field);
                    break;
                case REALSXP:
                    REAL(elt)[midx] = atof(field);
                    break;
                case STRSXP:
                    SET_STRING_ELT(elt, midx, mkChar(field));
                    break;
                default:
                    Rf_error("(internal) unhandled type '%s'",
                             type2char(TYPEOF(elt)));
                }
            }
        }

        irec += 1;
    }

    _vcf_grow(vcf, irec, col_n);
    if (NILSXP == TYPEOF(inms))
        inms = mkString("INFO");
    SET_VECTOR_ELT(vcf, N_FLDS - 2, _trim_null(info, inms));
    SET_VECTOR_ELT(vcf, N_FLDS - 1, _trim_null(geno, gnms));
    _types_list_transpose(VECTOR_ELT(vcf, N_FLDS - 1));

    Free(buf);
    UNPROTECT(1);
    return vcf;
}
