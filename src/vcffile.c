#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include "vcffile.h"
#include "rle.h"
#include "dna_hash.h"
#include "utilities.h"
#include "IRanges_interface.h"

const struct fld_fmt {
    const char *name;
    SEXPTYPE type;
} FLD_FMT[] = {
    /* {"CHROM", VECSXP}, {"POS", INTSXP}, {"ID", STRSXP}, */
    /* {"REF", STRSXP}, */
    {"rowData", S4SXP}, {"REF", S4SXP}, {"ALT", STRSXP},
    {"QUAL", REALSXP}, {"FILTER", STRSXP}, {"INFO", STRSXP},
    {"GENO", VECSXP}
};

enum { ROWDATA_IDX = 0, REF_IDX, ALT_IDX, QUAL_IDX, FILTER_IDX,
       INFO_IDX, GENO_IDX };
enum { POS_IDX = 0, ID_IDX };

const int N_FLDS = sizeof(FLD_FMT) / sizeof(FLD_FMT[0]);

struct vcf_parse_t {
    SEXP vcf, info, geno;
    const char **inms, **gnms;
    int vcf_n, imap_n, gmap_n, samp_n, *gmapidx;
    struct rle_t *chrom;
    struct dna_hash_t *ref;
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

static SEXP _types_transpose(SEXP types)
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
    return types;
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
    /* allocate result and fixed fields */
    SEXP vcf, rowData, info, geno, elt, eltnms;
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

    /* fixed fields */
    SET_VECTOR_ELT(vcf, ROWDATA_IDX, Rf_allocVector(VECSXP, 2));
    rowData = VECTOR_ELT(vcf, ROWDATA_IDX);
    SET_VECTOR_ELT(rowData, POS_IDX, Rf_allocVector(INTSXP, vcf_n));
    SET_VECTOR_ELT(rowData, ID_IDX, Rf_allocVector(STRSXP, vcf_n));

    SET_VECTOR_ELT(vcf, REF_IDX, R_NilValue);

    for (int i = ALT_IDX; i <= FILTER_IDX; ++i) {
        elt = R_NilValue;
        if (R_NilValue != VECTOR_ELT(fmap, i))
            elt = Rf_allocVector(FLD_FMT[i].type, vcf_n);
        SET_VECTOR_ELT(vcf, i, elt);
    }

    /* info */
    info = _types_alloc(vcf_n, 1, imap, R_NilValue);
    SET_VECTOR_ELT(vcf, INFO_IDX, info);

    /* _transposed_ GENO */
    PROTECT(eltnms = Rf_allocVector(VECSXP, 2));
    SET_VECTOR_ELT(eltnms, 0, sample);
    SET_VECTOR_ELT(eltnms, 1, R_NilValue);
    geno = _types_alloc(samp_n, vcf_n, gmap, eltnms);
    SET_VECTOR_ELT(vcf, GENO_IDX, geno);
    UNPROTECT(1);

    UNPROTECT(1);
    return vcf;
}

static void _vcf_grow(SEXP vcf, const int vcf_n, const int samp_n)
{
    SEXP rowData, elt;

    rowData = VECTOR_ELT(vcf, ROWDATA_IDX);
    for (int i = POS_IDX; i <= ID_IDX; ++i) {
        elt = VECTOR_ELT(rowData, i);
        SET_VECTOR_ELT(rowData, i, Rf_lengthgets(elt, vcf_n));
    }

    for (int i = ALT_IDX; i <= FILTER_IDX; ++i) {
        elt = VECTOR_ELT(vcf, i);
        if (R_NilValue != elt)
            SET_VECTOR_ELT(vcf, i, Rf_lengthgets(elt, vcf_n));
    }

    elt = VECTOR_ELT(vcf, INFO_IDX);
    _types_grow(elt, vcf_n, 1);

    /* _transposed_ matrix */
    elt = VECTOR_ELT(vcf, GENO_IDX);
    _types_grow(elt, samp_n, vcf_n);
}

static void _vcf_parse(char *line, const int irec,
                       const struct vcf_parse_t *param)
{
    SEXP vcf = param->vcf, info = param->info, geno=param->geno,
        rowData, elt;
    const int samp_n = param->samp_n, imap_n = param->imap_n,
        gmap_n = param->gmap_n;
    const char **inms = param->inms, **gnms = param->gnms;
    int fmtidx, sampleidx, imapidx, *gmapidx = param->gmapidx;

    int idx = irec, j;
    struct it_t it0, it1, it2;
    char *sample, *field, *ifld, *ikey, *fmt;

    /* fixed fields */
    char *chrom = it_init(&it0, line, '\t'); /* CHROM */
    rle_append(param->chrom, chrom);
    rowData = VECTOR_ELT(vcf, ROWDATA_IDX);
    field = it_next(&it0);      /* POS */
    INTEGER(VECTOR_ELT(rowData, POS_IDX))[irec] = atoi(field);
    field = it_next(&it0);      /* ID */
    if ('.' == *field && '\0' == *(field + 1)) {
        /* construct ID if missing */
        const int len = strlen(chrom);
        field = chrom;
        *(field + len) = ':';   /* chrom\0pos\0 ==> chrom:pos\0 */
    }
    SET_STRING_ELT(VECTOR_ELT(rowData, ID_IDX), irec, mkChar(field));
    field = it_next(&it0);      /* REF */
    dna_hash_append(param->ref, field);
    for (field = it_next(&it0), j = ALT_IDX; j <= FILTER_IDX;
         field = it_next(&it0), ++j)
    {
        elt = VECTOR_ELT(vcf, j);
        switch (TYPEOF(elt)) {
        case NILSXP:
            break;
        case INTSXP:
            INTEGER(elt)[idx] = atoi(field);
            break;
        case REALSXP:
            REAL(elt)[idx] = ('.' == *field) ? R_NaReal : atof(field);
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
        for (ifld = it_init(&it1, field, ';'); '\0' != *ifld;
             ifld = it_next(&it1)) {
            ikey = it_init(&it2, ifld, '=');
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
                field = it_next(&it2);
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
    field = it_next(&it0);
    fmt = field;
    for (field = it_init(&it2, fmt, ':'), fmtidx = 0;
         '\0' != *field; field = it_next(&it2), fmtidx++) {
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
    for (sample = it_next(&it0), sampleidx = 0;
         '\0' != *sample; sample = it_next(&it0), sampleidx++) {
        for (field = it_init(&it2, sample, ':'), fmtidx = 0;
             '\0' != *field; field = it_next(&it2), fmtidx++) {
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

static void _vcf_as_S4(struct vcf_parse_t *param)
{
    /* ref: DNAStringSet */
    SEXP dna = dna_hash_as_DNAStringSet(param->ref);
    SET_VECTOR_ELT(param->vcf, REF_IDX, dna);

    /* rowData: GRanges */
    SEXP rowData, seqnames, start, width, names;
    PROTECT(seqnames = rle_as_Rle(param->chrom));
    rowData = VECTOR_ELT(param->vcf, ROWDATA_IDX);
    start = VECTOR_ELT(rowData, POS_IDX);
    names = VECTOR_ELT(rowData, ID_IDX);
    width = get_XVectorList_width(dna);

    SEXP ranges, nmspc, fun, expr;
    PROTECT(ranges =  new_IRanges("IRanges", start, width, names));
    PROTECT(nmspc = get_namespace("GenomicRanges"));
    PROTECT(fun = findFun(install("GRanges"), nmspc));
    PROTECT(expr = lang3(fun, seqnames, ranges));

    SET_VECTOR_ELT(param->vcf, ROWDATA_IDX, eval(expr, R_GlobalEnv));
    UNPROTECT(5);
}

static void _vcf_types_tidy(struct vcf_parse_t *param)
{
    if (NULL == param->inms) {
        param->inms = (const char **) R_alloc(sizeof(const char *), 1);
        param->inms[0] = "INFO";
    }
    param->info = _trim_null(param->info, param->inms);
    SET_VECTOR_ELT(param->vcf, INFO_IDX, param->info);

    param->geno = _trim_null(param->geno, param->gnms);
    SET_VECTOR_ELT(param->vcf, GENO_IDX, param->geno);

    param->geno = _types_transpose(param->geno);
    SET_VECTOR_ELT(param->vcf, GENO_IDX, param->geno);
}

static struct vcf_parse_t *_vcf_parse_new(int vcf_n, SEXP sample,
                                          SEXP fmap, SEXP imap,
                                          SEXP gmap)
{
    struct vcf_parse_t *param = Calloc(1, struct vcf_parse_t);
    SEXP elt;

    param->vcf_n = vcf_n;

    elt = _vcf_allocate(param->vcf_n, sample, fmap, imap, gmap);
    PROTECT(param->vcf = elt);
    param->info = VECTOR_ELT(param->vcf, INFO_IDX);
    param->geno = VECTOR_ELT(param->vcf, GENO_IDX);

    param->imap_n = Rf_length(imap);
    param->inms =
        (const char **) R_alloc(sizeof(const char *), param->imap_n);
    for (int j = 0; j < param->imap_n; ++j)
        param->inms[j] = CHAR(STRING_ELT(GET_NAMES(imap), j));

    param->samp_n = Rf_length(sample);
    param->gmap_n = Rf_length(gmap);
    param->gnms =
        (const char **) R_alloc(sizeof(const char *), param->gmap_n);
    for (int j = 0; j < param->gmap_n; ++j)
        param->gnms[j] = CHAR(STRING_ELT(GET_NAMES(gmap), j));
    param->gmapidx = (int *) R_alloc(sizeof(int), param->gmap_n);
    param->chrom = rle_new(param->vcf_n);
    param->ref = dna_hash_new(param->vcf_n);

    UNPROTECT(1);
    return param;
}

static void _vcf_parse_free(struct vcf_parse_t *param)
{
    rle_free(param->chrom);
    dna_hash_free(param->ref);
    Free(param);
}

static void _vcf_parse_grow(struct vcf_parse_t *param, int size)
{
    static const double SCALE = 1.6;
    if (0 == size)
        size = param->vcf_n < 2 ? 2 : param->vcf_n * SCALE;
    _vcf_grow(param->vcf, size, param->samp_n);
    param->info = VECTOR_ELT(param->vcf, N_FLDS - 2);
    param->geno = VECTOR_ELT(param->vcf, N_FLDS - 1);
    param->vcf_n = size;
}

SEXP scan_vcf_connection(SEXP txt, SEXP sample, SEXP fmap, SEXP imap,
                         SEXP gmap)
{
    struct vcf_parse_t *param;
    SEXP result;

    PROTECT(result = Rf_allocVector(VECSXP, 1));
    param = _vcf_parse_new(Rf_length(txt), sample, fmap, imap, gmap);
    SET_VECTOR_ELT(result, 0, param->vcf);

    /* parse each line */
    for (int irec = 0; irec < param->vcf_n; irec++) {
        char *line = strdup(CHAR(STRING_ELT(txt, irec)));
        _vcf_parse(line, irec, param);
        free(line);
    }

    _vcf_as_S4(param);
    _vcf_types_tidy(param);
    _vcf_parse_free(param);

    UNPROTECT(1);
    return result;
}

SEXP scan_vcf_character(SEXP file, SEXP yield,
                        /* SEXP grow, */
                        SEXP sample, SEXP fmap, SEXP imap, SEXP gmap)
{
    struct vcf_parse_t *param;
    SEXP result;

    if (!IS_INTEGER(yield) || 1L != Rf_length(yield))
        Rf_error("'yield' must be integer(1)");
    if (!IS_CHARACTER(file) || 1L != Rf_length(file))
        Rf_error("'file' must be character(1) or as on ?scanVcf");

    PROTECT(result = Rf_allocVector(VECSXP, 1));
    param = _vcf_parse_new(INTEGER(yield)[0], sample, fmap, imap,
                           gmap);
    SET_VECTOR_ELT(result, 0, param->vcf);

    const int BUFLEN = 4096;
    char *buf0 = Calloc(BUFLEN, char);
    char *buf = buf0, *end = buf0 + BUFLEN;

    gzFile gz = gzopen(CHAR(STRING_ELT(file, 0)), "rb");
    int curr, prev, irec = 0;

    if (Z_NULL == gz) {
        Free(param);
        Rf_error("failed to open file");
    }

    prev = gztell(gz);
    while (Z_NULL != gzgets(gz, buf, end - buf)) {
        curr = gztell(gz);
        if (curr - prev == end - buf0 - 1 && *(end - 1) == '\0') {
            const int len0 = end - buf0, len1 = len0 * 1.6;
            buf0 = Realloc(buf0, len1, char);
            buf = buf0 + len0 - 1;
            end = buf0 + len1;
            continue;
        }
        if ('#' == *buf || '\0' == *buf)
            continue;

        if (irec == param->vcf_n) {
            /* if (!grow_b) */
            /*     break; */
            _vcf_parse_grow(param, 0);
        }

        _vcf_parse(buf, irec, param);
        irec += 1;

        prev = curr;
        buf = buf0;
    }

    gzclose(gz);
    free(buf0);

    _vcf_grow(param->vcf, irec, param->samp_n);
    _vcf_as_S4(param);
    _vcf_types_tidy(param);
    _vcf_parse_free(param);

    UNPROTECT(1);
    return result;
}

SEXP tabix_as_vcf(tabix_t *tabix, ti_iter_t iter,
                  const int *keep, const int size,
                  const Rboolean grow, SEXP state)
{
    struct vcf_parse_t *param;
    SEXP result;

    param = _vcf_parse_new(size, VECTOR_ELT(state, 0),
                           VECTOR_ELT(state, 1), VECTOR_ELT(state, 2),
                           VECTOR_ELT(state, 3));
    PROTECT(result = param->vcf);

    int BUFLEN = 4096;
    char *buf = Calloc(BUFLEN, char);

    int linelen;
    const char *line;

    int irec = 0, trec=1;
    while (NULL != (line = ti_read(tabix, iter, &linelen))) {

        if (irec == param->vcf_n) {
            if (!grow)
                break;
            _vcf_parse_grow(param, 0);
        }

        if (NULL != keep) {
            if (trec++ != *keep) continue;
            keep++;
        }

        if (linelen + 1 > BUFLEN) {
            Free(buf);
            BUFLEN = 2 * linelen;
            buf = Calloc(BUFLEN, char);
        }
        memcpy(buf, line, linelen);
        buf[linelen] = '\0';

        _vcf_parse(buf, irec, param);
        irec += 1;
    }
    Free(buf);

    _vcf_grow(param->vcf, irec, param->samp_n);
    _vcf_as_S4(param);
    _vcf_types_tidy(param);
    _vcf_parse_free(param);

    UNPROTECT(1);
    return result;
}
