#include <stdlib.h>
#include <string.h>
#include "vcffile.h"

/* iterator to return null-terminated delimited fields */

struct it {
    char *str;
    char delim;
};

static inline char *_it_next(struct it *it)
{
    char *curr = it->str;
    while ('\0' != *it->str && it->delim != *it->str)
        it->str++;
    if ('\0' != *it->str)
        *it->str++ = '\0';
    return curr;
}

static char *_it_init(struct it *it, char *str, char delim)
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
    {"CHROM", STRSXP}, {"POS", INTSXP}, {"ID", STRSXP},
    {"REF", STRSXP}, {"ALT", STRSXP}, {"QUAL", REALSXP},
    {"FILTER", STRSXP}, {"INFO", STRSXP}, {"GENO", VECSXP}
};

const int N_FLDS = sizeof(FLD_FMT) / sizeof(FLD_FMT[0]);

struct vcf_parse_t {
    SEXP vcf, info, geno;
    const char **inms, **gnms;
    int imap_n, gmap_n, samp_n, *gmapidx;
};

static SEXP _types_alloc(const int vcf_n, const int col_n,
                         SEXP map, SEXP eltnms)
{
    int i, j, map_n = Rf_length(map);
    SEXP types, elt;

    /* case of no INFO or GENO information in header */
    if (map_n == 0) {
        PROTECT(types = Rf_allocVector(VECSXP, 1));
        elt = Rf_allocMatrix(STRSXP, vcf_n, 1);
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
            elt = Rf_allocMatrix(type, vcf_n, col_n);
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

static void _types_grow(SEXP types, const int vcf_n, const int col_n)
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

static void _types_transpose(SEXP types)
{
    SEXP elt, tpose, fun;
    for (int t = 0; t < Rf_length(types); ++t) {
        elt = VECTOR_ELT(types, t);
        if (R_NilValue == elt)
            continue;
        fun = PROTECT(findFun(install("t"), R_BaseEnv));
        tpose = eval(lang2(fun, elt), R_BaseEnv);
        SET_VECTOR_ELT(types, t, tpose);
        UNPROTECT(1);
    }
}

static SEXP _trim_null(SEXP data, const char **cnms)
{
    SEXP nms = PROTECT(NEW_CHARACTER(Rf_length(data)));
    int i, j = 0;
    for (i = 0; i < Rf_length(data); ++i) {
        if (R_NilValue != VECTOR_ELT(data, i)) {
            SET_VECTOR_ELT(data, j, VECTOR_ELT(data, i));
            SET_STRING_ELT(nms, j, mkChar(cnms[i]));
            j++;
        }
    }
    PROTECT(nms = Rf_lengthgets(nms, j));
    PROTECT(data = Rf_lengthgets(data, j));
    data = Rf_namesgets(data, nms);
    UNPROTECT(3);

    return data;
}

static SEXP _vcf_allocate(const int vcf_n, SEXP sample,
                          SEXP fmap, SEXP imap, SEXP gmap)
{
    /* allocate result and first 7 fixed fields */
    SEXP vcf, info, geno, elt, eltnms;
    const int samp_n = Rf_length(sample);

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
    info = _types_alloc(vcf_n, 1, imap, R_NilValue);
    SET_VECTOR_ELT(vcf, N_FLDS - 2, info);

    /* allocate _transposed_ GENO */
    PROTECT(eltnms = Rf_allocVector(VECSXP, 2));
    SET_VECTOR_ELT(eltnms, 0, sample);
    SET_VECTOR_ELT(eltnms, 1, R_NilValue);
    geno = _types_alloc(samp_n, vcf_n, gmap, eltnms);
    SET_VECTOR_ELT(vcf, N_FLDS - 1, geno);
    UNPROTECT(1);

    UNPROTECT(1);
    return vcf;
}

static void _vcf_grow(SEXP vcf, const int vcf_n, const int samp_n)
{
    SEXP elt;

    for (int i = 0; i < N_FLDS - 2; ++i) {
        elt = VECTOR_ELT(vcf, i);
        if (R_NilValue != elt)
            SET_VECTOR_ELT(vcf, i, Rf_lengthgets(elt, vcf_n));
    }

    elt = VECTOR_ELT(vcf, N_FLDS - 2);
    _types_grow(elt, vcf_n, 1);

    /* _transposed_ matrix */
    elt = VECTOR_ELT(vcf, N_FLDS - 1);
    _types_grow(elt, samp_n, vcf_n);
}

static void _vcf_parse(char *line, const int irec,
                       const struct vcf_parse_t *param)
{
    SEXP vcf = param->vcf, info = param->info, geno=param->geno, elt;
    const int samp_n = param->samp_n, imap_n = param->imap_n,
        gmap_n = param->gmap_n;
    const char **inms = param->inms, **gnms = param->gnms;
    int fmtidx, sampleidx, imapidx, *gmapidx = param->gmapidx;

    int idx = irec, j;
    struct it it0, it1, it2;
    char *sample, *field, *ifld, *ikey, *fmt;
    char *dot = ".";

    /* first 7 'fixed' fields */
    for (field = _it_init(&it0, line, '\t'), j = 0;
         j < N_FLDS - 2; field = _it_next(&it0), ++j) {
        elt = VECTOR_ELT(vcf, j);
        switch (TYPEOF(elt)) {
        case NILSXP:
            break;
        case INTSXP:
            INTEGER(elt)[idx] = atoi(field);
            break;
        case REALSXP:
            if (strcmp(field, dot) == 0)
                REAL(elt)[idx] = R_NaReal;
            else
                REAL(elt)[idx] = atof(field);
            break;
        case STRSXP:
            SET_STRING_ELT(elt, idx, mkChar(field));
            break;
        default:
            Rf_error("(internal) unhandled fixed field type '%s'",
                     type2char(TYPEOF(elt)));
        }
    }

    /* 'INFO' field */
    if (imap_n == 0) {          /* no header; parse as char */
        elt = VECTOR_ELT(info, 0);
        SET_STRING_ELT(elt, idx, mkChar(field));
    } else {
        for (ifld = _it_init(&it1, field, ';'); '\0' != *ifld;
             ifld = _it_next(&it1)) {
            ikey = _it_init(&it2, ifld, '=');
            for (imapidx = 0; imapidx < imap_n; ++imapidx) {
                if (0L == strcmp(ikey, inms[imapidx]))
                    break;
            }
            if (imap_n == imapidx)
                Rf_error("record %d INFO '%s' not found", idx + 1,
                         ikey);

            elt = VECTOR_ELT(info, imapidx);
            if (LGLSXP == TYPEOF(elt)) {
                LOGICAL(elt)[idx] = TRUE;
            } else {
                field = _it_next(&it2);
                switch (TYPEOF(elt)) {
                case NILSXP:
                    break;
                case INTSXP:
                    INTEGER(elt)[idx] = atoi(field);
                    break;
                case REALSXP:
                    REAL(elt)[idx] = atof(field);
                    break;
                case STRSXP:
                    SET_STRING_ELT(elt, idx, mkChar(field));
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
            if (0L == strcmp(field, gnms[j]))
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
            idx = irec * samp_n + sampleidx;
            switch (TYPEOF(elt)) {
            case NILSXP:
                break;
            case INTSXP:
                INTEGER(elt)[idx] = atoi(field);
                break;
            case REALSXP:
                REAL(elt)[idx] = atof(field);
                break;
            case STRSXP:
                SET_STRING_ELT(elt, idx, mkChar(field));
                break;
            default:
                Rf_error("(internal) unhandled type '%s'",
                         type2char(TYPEOF(elt)));
            }
        }
    }
}

SEXP scan_vcf_connection(SEXP txt, SEXP sample, SEXP fmap, SEXP imap,
                         SEXP gmap)
{
    struct vcf_parse_t param;
    const int vcf_n = Rf_length(txt);
    SEXP result = PROTECT(Rf_allocVector(VECSXP, 1));

    param.vcf = _vcf_allocate(vcf_n, sample, fmap, imap, gmap);
    SET_VECTOR_ELT(result, 0, param.vcf);

    param.info = VECTOR_ELT(param.vcf, N_FLDS - 2);
    param.imap_n = Rf_length(imap);
    param.inms = (const char **) R_alloc(sizeof(char *), param.imap_n);
    for (int j = 0; j < param.imap_n; ++j)
        param.inms[j] = CHAR(STRING_ELT(GET_NAMES(imap), j));

    param.samp_n = Rf_length(sample);
    param.geno = VECTOR_ELT(param.vcf, N_FLDS - 1);
    param.gmap_n = Rf_length(gmap);
    param.gnms = (const char **) R_alloc(sizeof(char *), param.gmap_n);
    for (int j = 0; j < param.gmap_n; ++j)
        param.gnms[j] = CHAR(STRING_ELT(GET_NAMES(gmap), j));
    param.gmapidx = (int *) R_alloc(sizeof(int), param.gmap_n);

    /* parse each line */
    for (int irec = 0; irec < vcf_n; irec++) {
        char *line = strdup(CHAR(STRING_ELT(txt, irec)));
        _vcf_parse(line, irec, &param);
        free(line);
    }

    if (NULL == param.inms) {
        param.inms = (const char **) R_alloc(sizeof(const char *), 1);
        param.inms[0] = "INFO";
    }
    SET_VECTOR_ELT(param.vcf, N_FLDS - 2,
                   _trim_null(param.info, param.inms));
    SET_VECTOR_ELT(param.vcf, N_FLDS - 1,
                   _trim_null(param.geno, param.gnms));
    _types_transpose(VECTOR_ELT(param.vcf, N_FLDS - 1));

    UNPROTECT(1);
    return result;
}

SEXP tabix_as_vcf(tabix_t *tabix, ti_iter_t iter,
                  const int *keep, const int size,
                  const Rboolean grow, SEXP state)
{
    SEXP sample = VECTOR_ELT(state, 0), fmap = VECTOR_ELT(state, 1),
        imap = VECTOR_ELT(state, 2), gmap = VECTOR_ELT(state, 3);
    int vcf_n = size;
    const int samp_n = Rf_length(sample);
    struct vcf_parse_t param;

    param.vcf = _vcf_allocate(vcf_n, sample, fmap, imap, gmap);
    PROTECT(param.vcf);

    param.info = VECTOR_ELT(param.vcf, N_FLDS - 2);
    param.imap_n = Rf_length(imap);
    param.inms =
        (const char **) R_alloc(sizeof(const char *), param.imap_n);
    for (int j = 0; j < param.imap_n; ++j)
        param.inms[j] = CHAR(STRING_ELT(GET_NAMES(imap), j));

    param.samp_n = Rf_length(sample);
    param.geno = VECTOR_ELT(param.vcf, N_FLDS - 1);
    param.gmap_n = Rf_length(gmap);
    param.gnms =
        (const char **) R_alloc(sizeof(const char *), param.gmap_n);
    for (int j = 0; j < param.gmap_n; ++j)
        param.gnms[j] = CHAR(STRING_ELT(GET_NAMES(gmap), j));
    param.gmapidx = (int *) R_alloc(sizeof(int), param.gmap_n);

    const double SCALE = 1.6;
    int buflen = 4096;
    char *buf = Calloc(buflen, char);

    int linelen;
    const char *line;

    int irec = 0, trec=1;
    while (NULL != (line = ti_read(tabix, iter, &linelen))) {

        if (irec == vcf_n) {
            if (!grow)
                break;
            vcf_n *= SCALE;
            _vcf_grow(param.vcf, vcf_n, samp_n);
            param.info = VECTOR_ELT(param.vcf, N_FLDS - 2);
            param.geno = VECTOR_ELT(param.vcf, N_FLDS - 1);
        }

        if (NULL != keep) {
            if (trec++ != *keep) continue;
            keep++;
        }

        if (linelen + 1 > buflen) {
            Free(buf);
            buflen = 2 * linelen;
            buf = Calloc(buflen, char);
        }
        memcpy(buf, line, linelen);
        buf[linelen] = '\0';

        _vcf_parse(buf, irec, &param);
        irec += 1;
    }

    _vcf_grow(param.vcf, irec, samp_n);
    if (NULL == param.inms) {
        param.inms = (const char **) R_alloc(sizeof(const char *), 1);
        param.inms[0] = "INFO";
    }
    SET_VECTOR_ELT(param.vcf, N_FLDS - 2,
                   _trim_null(param.info, param.inms));
    SET_VECTOR_ELT(param.vcf, N_FLDS - 1,
                   _trim_null(param.geno, param.gnms));
    _types_transpose(VECTOR_ELT(param.vcf, N_FLDS - 1));

    Free(buf);
    UNPROTECT(1);
    return param.vcf;
}
