#include <zlib.h>
#include "vcffile.h"
#include "vcftype.h"
#include "rle.h"
#include "dna_hash.h"
#include "utilities.h"
#include "IRanges_interface.h"

enum { ROWDATA_IDX = 0, REF_IDX, ALT_IDX, QUAL_IDX, FILTER_IDX,
       INFO_IDX, GENO_IDX };
enum { POS_IDX = 0, ID_IDX };

const int N_FLDS = 7;

struct parse_t {
    struct vcftype_t *vcf;
    struct rle_t *chrom;
    struct dna_hash_t *ref;
    int vcf_n, imap_n, gmap_n, samp_n, *gmapidx;
    const char **inms, **gnms;
};

static struct vcftype_t *_types_alloc(const int vcf_n, const int col_n,
                                      SEXP map, const Rboolean isArray)
{
    struct vcftype_t *types;
    const int map_n = Rf_length(map);

    if (map_n == 0) {           /* no INFO or GENO in header */
        types = _vcftype_new(VECSXP, 0, 0, FALSE);
    } else {
        types = _vcftype_new(VECSXP, map_n, 0, FALSE);
        for (int j = 0; j < map_n; ++j) {
            SEXPTYPE type = TYPEOF(VECTOR_ELT(map, j));
            types->u.list[j] = _vcftype_new(type, vcf_n, col_n, isArray);
        }
    }

    return types;
}

static void _types_grow(struct vcftype_t *types, int vcf_n)
{
    for (int j = 0; j < types->nrow; ++j)
        types->u.list[j] = _vcftype_grow(types->u.list[j], vcf_n);
}

static SEXP _trim_null(SEXP data, const char **cnms)
{
    SEXP nms = PROTECT(NEW_CHARACTER(Rf_length(data)));
    int j = 0;
    for (int i = 0; i < Rf_length(data); ++i) {
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

static struct vcftype_t *_vcf_alloc(const int vcf_n, SEXP sample,
                                    SEXP fmap, SEXP imap, SEXP gmap)
{
    struct vcftype_t *vcf, *rowData;
    const int samp_n = Rf_length(sample);

    if (Rf_length(fmap) != N_FLDS - 2)
        Rf_error("[internal] 'fixed' field length %d does not equal %d",
                 Rf_length(fmap), N_FLDS - 2);
    vcf = _vcftype_new(VECSXP, N_FLDS, 0, FALSE);

    /* fixed fields */
    rowData = _vcftype_new(VECSXP, 2, 0, FALSE);
    rowData->u.list[POS_IDX] = _vcftype_new(INTSXP, vcf_n, 0, FALSE);
    rowData->u.list[ID_IDX] = _vcftype_new(STRSXP, vcf_n, 0, FALSE);
    vcf->u.list[ROWDATA_IDX] = rowData;

    for (int i = ALT_IDX; i <= FILTER_IDX; ++i) {
        SEXPTYPE type = TYPEOF(VECTOR_ELT(fmap, i));
        vcf->u.list[i] = _vcftype_new(type, vcf_n, 0, FALSE);
    }

    /* info, geno */
    vcf->u.list[INFO_IDX] = _types_alloc(vcf_n, 1, imap, TRUE);
    vcf->u.list[GENO_IDX] = _types_alloc(vcf_n, samp_n, gmap, TRUE);

    return vcf;
}

static void _vcf_grow(struct vcftype_t * vcf, int vcf_n)
{
    struct vcftype_t *elt;

    elt = vcf->u.list[ROWDATA_IDX];
    for (int i = POS_IDX; i <= ID_IDX; ++i)
        elt->u.list[i] = _vcftype_grow(elt->u.list[i], vcf_n);

    for (int i = ALT_IDX; i <= FILTER_IDX; ++i)
        vcf->u.list[i] = _vcftype_grow(vcf->u.list[i], vcf_n);

    _types_grow(vcf->u.list[INFO_IDX], vcf_n);
    _types_grow(vcf->u.list[GENO_IDX], vcf_n);
}

static void _parse(char *line, const int irec,
                   const struct parse_t *parse)
{
    struct vcftype_t *vcf = parse->vcf, *rowData, *elt;
    const int samp_n = parse->samp_n, imap_n = parse->imap_n,
        gmap_n = parse->gmap_n;
    const char **inms = parse->inms, **gnms = parse->gnms;
    int fmtidx, sampleidx, imapidx, *gmapidx = parse->gmapidx;

    int idx = irec, j;
    struct it_t it0, it1, it2;
    char *sample, *field, *ifld, *ikey, *fmt;

    /* fixed fields */
    char *chrom = it_init(&it0, line, '\t'); /* CHROM */
    rle_append(parse->chrom, chrom);
    rowData = vcf->u.list[ROWDATA_IDX];
    field = it_next(&it0);      /* POS */
    rowData->u.list[POS_IDX]->u.integer[irec] = atoi(field);
    field = it_next(&it0);      /* ID */
    if ('.' == *field && '\0' == *(field + 1)) {
        /* construct ID if missing: chrom\0pos\0 ==> chrom:pos\0 */
        field = chrom;
        *(field + strlen(chrom)) = ':';
    }
    rowData->u.list[ID_IDX]->u.character[irec] = Strdup(field);

    field = it_next(&it0);      /* REF */
    dna_hash_append(parse->ref, field);
    for (field = it_next(&it0), j = ALT_IDX; j <= FILTER_IDX;
         field = it_next(&it0), ++j)
    {
        elt = vcf->u.list[j];
        _vcftype_set(elt, idx, field);
    }

    /* INFO */
    struct vcftype_t *info = vcf->u.list[INFO_IDX];
    if (1 == imap_n && NULL == inms) { /* no header; parse as char */
        elt = info->u.list[0];
        elt->u.character[idx] = Strdup(field);
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

            elt = info->u.list[imapidx];
            if (LGLSXP == elt->type) {
                elt->u.logical[idx] = TRUE;
            } else {
                field = it_next(&it2);
                _vcftype_set(elt, idx, field);
            }
        }
    }

    /* FORMAT */
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

    /* sample(s) */
    struct vcftype_t *geno = vcf->u.list[GENO_IDX];
    for (sample = it_next(&it0), sampleidx = 0;
         '\0' != *sample; sample = it_next(&it0), sampleidx++) {
        for (field = it_init(&it2, sample, ':'), fmtidx = 0;
             '\0' != *field; field = it_next(&it2), fmtidx++) {
            elt = geno->u.list[ gmapidx[fmtidx] ];
            idx = irec * samp_n + sampleidx;
            _vcftype_set(elt, idx, field);
        }
    }
}

static SEXP _vcf_as_SEXP(struct parse_t *parse, SEXP fmap, SEXP sample)
{
    SEXP result = PROTECT(_vcftype_as_SEXP(parse->vcf));

    /* ref: DNAStringSet */
    SEXP dna = dna_hash_as_DNAStringSet(parse->ref);
    SET_VECTOR_ELT(result, REF_IDX, dna);

    /* rowData: GRanges */
    SEXP rowData, seqnames, start, width, names;
    PROTECT(seqnames = rle_as_Rle(parse->chrom));
    rowData = VECTOR_ELT(result, ROWDATA_IDX);
    start = VECTOR_ELT(rowData, POS_IDX);
    names = VECTOR_ELT(rowData, ID_IDX);
    width = get_XVectorList_width(dna);

    SEXP ranges, nmspc, fun, expr;
    PROTECT(ranges =  new_IRanges("IRanges", start, width, names));
    PROTECT(nmspc = get_namespace("GenomicRanges"));
    PROTECT(fun = findFun(install("GRanges"), nmspc));
    PROTECT(expr = lang3(fun, seqnames, ranges));

    SET_VECTOR_ELT(result, ROWDATA_IDX, eval(expr, R_GlobalEnv));
    UNPROTECT(5);

    /* names */
    SEXP nms, sxp = Rf_getAttrib(fmap, R_NamesSymbol), elt;
    PROTECT(nms = Rf_allocVector(STRSXP, Rf_length(result)));
    for (int i = ROWDATA_IDX; i <= FILTER_IDX; ++i)
        SET_STRING_ELT(nms, i, STRING_ELT(sxp, i));
    SET_STRING_ELT(nms, INFO_IDX, mkChar("INFO"));
    SET_STRING_ELT(nms, GENO_IDX, mkChar("GENO"));
    Rf_namesgets(result, nms);
    UNPROTECT(1);

    PROTECT(nms = Rf_allocVector(STRSXP, parse->imap_n));
    if (1 == parse->imap_n && NULL == parse->inms)
        SET_STRING_ELT(nms, 0, R_NaString);
    else
        for (int i = 0; i < parse->imap_n; ++i)
            SET_STRING_ELT(nms, i, mkChar(parse->inms[i]));
    Rf_namesgets(VECTOR_ELT(result, INFO_IDX), nms);
    UNPROTECT(1);

    PROTECT(nms = Rf_allocVector(STRSXP, parse->gmap_n));
    for (int i = 0; i < parse->gmap_n; ++i)
        SET_STRING_ELT(nms, i, mkChar(parse->gnms[i]));
    Rf_namesgets(VECTOR_ELT(result, GENO_IDX), nms);
    UNPROTECT(1);

    PROTECT(nms = Rf_allocVector(VECSXP, 2));
    SET_VECTOR_ELT(nms, 0, R_NilValue);
    SET_VECTOR_ELT(nms, 1, sample);
    sxp = VECTOR_ELT(result, GENO_IDX);
    for (int i = 0; i < Rf_length(sxp); ++i) {
        elt = VECTOR_ELT(sxp, i);
        if (R_NilValue != elt)
            Rf_dimnamesgets(elt, nms);
    }
    UNPROTECT(1);

    UNPROTECT(1);
    return result;
}

static void _vcf_types_tidy(struct parse_t *parse, SEXP result)
{
    SEXP elt;
    if (NULL == parse->inms) {
        parse->inms = (const char **) R_alloc(sizeof(const char *), 1);
        parse->inms[0] = "INFO";
    }
    elt = VECTOR_ELT(result, INFO_IDX);
    SET_VECTOR_ELT(result, INFO_IDX, _trim_null(elt, parse->inms));

    elt = VECTOR_ELT(result, GENO_IDX);
    SET_VECTOR_ELT(result, GENO_IDX, _trim_null(elt, parse->gnms));
}

static struct parse_t *_parse_new(int vcf_n, SEXP sample, SEXP fmap,
                                  SEXP imap, SEXP gmap)
{
    struct parse_t *parse = Calloc(1, struct parse_t);

    parse->vcf_n = vcf_n;
    parse->vcf = _vcf_alloc(parse->vcf_n, sample, fmap, imap, gmap);

    parse->imap_n = Rf_length(imap);
    if (1 == parse->imap_n && R_NilValue == GET_NAMES(imap))
        parse->inms = NULL;
    else {
        parse->inms =
            (const char **) R_alloc(sizeof(const char *), parse->imap_n);
        for (int j = 0; j < parse->imap_n; ++j)
            parse->inms[j] = CHAR(STRING_ELT(GET_NAMES(imap), j));
    }

    parse->samp_n = Rf_length(sample);
    parse->gmap_n = Rf_length(gmap);
    parse->gnms =
        (const char **) R_alloc(sizeof(const char *), parse->gmap_n);
    for (int j = 0; j < parse->gmap_n; ++j)
        parse->gnms[j] = CHAR(STRING_ELT(GET_NAMES(gmap), j));
    parse->gmapidx = (int *) R_alloc(sizeof(int), parse->gmap_n);
    parse->chrom = rle_new(parse->vcf_n);
    parse->ref = dna_hash_new(parse->vcf_n);

    return parse;
}

static void _parse_free(struct parse_t *parse)
{
    rle_free(parse->chrom);
    dna_hash_free(parse->ref);
    Free(parse);
}

static void _parse_grow(struct parse_t *parse, int size)
{
    static const double SCALE = 1.6;
    if (0 == size)
        size = parse->vcf_n < 2 ? 2 : parse->vcf_n * SCALE;
    _vcf_grow(parse->vcf, size);
    parse->vcf_n = size;
}

SEXP scan_vcf_connection(SEXP txt, SEXP sample, SEXP fmap, SEXP imap,
                         SEXP gmap)
{
    struct parse_t *parse;

    parse = _parse_new(Rf_length(txt), sample, fmap, imap, gmap);

    /* parse each line */
    for (int irec = 0; irec < parse->vcf_n; irec++) {
        char *line = Strdup(CHAR(STRING_ELT(txt, irec)));
        _parse(line, irec, parse);
        Free(line);
    }

    SEXP result = PROTECT(Rf_allocVector(VECSXP, 1));
    SET_VECTOR_ELT(result, 0, _vcf_as_SEXP(parse, fmap, sample));
    _vcf_types_tidy(parse, result);
    _parse_free(parse);
    UNPROTECT(1);

    return result;
}

SEXP scan_vcf_character(SEXP file, SEXP yield,
                        /* SEXP grow, */
                        SEXP sample, SEXP fmap, SEXP imap, SEXP gmap)
{
    struct parse_t *parse;

    if (!IS_INTEGER(yield) || 1L != Rf_length(yield))
        Rf_error("'yield' must be integer(1)");
    if (!IS_CHARACTER(file) || 1L != Rf_length(file))
        Rf_error("'file' must be character(1) or as on ?scanVcf");

    parse = _parse_new(INTEGER(yield)[0], sample, fmap, imap, gmap);

    const int BUFLEN = 4096;
    char *buf0 = Calloc(BUFLEN, char);
    char *buf = buf0, *end = buf0 + BUFLEN;

    gzFile gz = gzopen(CHAR(STRING_ELT(file, 0)), "rb");
    int irec = 0;
    if (Z_NULL == gz) {
        Free(parse);
        Rf_error("failed to open file");
    }

    while (Z_NULL != gzgets(gz, buf, end - buf)) {
        int  n = strlen(buf);
        if (n == end - buf - 1 && (*(end - 2) != '\n' || *(end - 2) != '\r')) {
            const int len0 = end - buf0, len1 = len0 * 1.6;
            buf0 = Realloc(buf0, len1, char);
            buf = buf0 + len0 - 1;
            end = buf0 + len1;
            continue;
        }
        if ('#' == *buf0 || '\0' == *buf0) {
            buf = buf0;
            continue;
        }

        if (irec == parse->vcf_n) {
            /* if (!grow_b) */
            /*     break; */
            _parse_grow(parse, 0);
        }

        /* trim trailing newlines */
        int last = strlen(buf) - 1;
        while (last >= 0) {
            if (buf[last] == '\n' || buf[last] == '\r')
                buf[last--] = '\0';
            else break;
        }

        _parse(buf0, irec, parse);

        irec += 1;
        buf = buf0;
    }

    gzclose(gz);
    Free(buf0);

    _vcf_grow(parse->vcf, irec);

    SEXP result = PROTECT(Rf_allocVector(VECSXP, 1));
    SET_VECTOR_ELT(result, 0, _vcf_as_SEXP(parse, fmap, sample));
    _vcf_types_tidy(parse, VECTOR_ELT(result, 0));
    _parse_free(parse);
    UNPROTECT(1);

    return result;
}

SEXP tabix_as_vcf(tabix_t *tabix, ti_iter_t iter,
                  const int *keep, const int size,
                  const Rboolean grow, SEXP state)
{
    SEXP sample = VECTOR_ELT(state, 0), fmap = VECTOR_ELT(state, 1);
    struct parse_t *parse =
        _parse_new(size, sample, fmap, VECTOR_ELT(state, 2),
                   VECTOR_ELT(state, 3));

    int BUFLEN = 4096;
    char *buf = Calloc(BUFLEN, char);

    int linelen;
    const char *line;

    int irec = 0, trec=1;
    while (NULL != (line = ti_read(tabix, iter, &linelen))) {

        if (irec == parse->vcf_n) {
            if (!grow) break;
            _parse_grow(parse, 0);
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

        _parse(buf, irec, parse);
        irec += 1;
    }
    Free(buf);

    _vcf_grow(parse->vcf, irec);

    SEXP result = PROTECT(_vcf_as_SEXP(parse, fmap, sample));
    _vcf_types_tidy(parse, result);
    _parse_free(parse);
    UNPROTECT(1);

    return result;
}
