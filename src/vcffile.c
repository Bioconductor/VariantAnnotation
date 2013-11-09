#include <stdarg.h>
#include <zlib.h>
#include "vcffile.h"
#include "vcftype.h"
#include "rle.h"
#include "dna_hash.h"
#include "utilities.h"
#include "IRanges_interface.h"
#include "XVector_interface.h"
#include "samtools/khash.h"
#include "strhash.h"

enum { ROWDATA_IDX = 0, REF_IDX, ALT_IDX, QUAL_IDX, FILTER_IDX,
       INFO_IDX, GENO_IDX };
enum { POS_IDX = 0, ID_IDX };

static const int N_FLDS = 7;
static const int TBX_INIT_SIZE = 32767;

KHASH_SET_INIT_STR(WARNINGS)

static khash_t(WARNINGS) *vcfwarn_new()
{
    return kh_init(WARNINGS);
}

static void vcfwarn(khash_t(WARNINGS) *warnings,
                    const char *fmt, ...)
{
    static const int bufsize = 2048;
    char *buf = Calloc(strlen(fmt), char);
    int ret;

    memcpy(buf, fmt, strlen(fmt) + 1);
    if (kh_get(WARNINGS, warnings, buf) != kh_end(warnings)) {
        Free(buf);
        return;
    }
    kh_put(WARNINGS, warnings, buf, &ret);

    buf = Calloc(bufsize, char);
    va_list argp;
    va_start(argp, fmt);
    (void) vsnprintf(buf, bufsize, fmt, argp);
    va_end(argp);
    Rf_warning(buf);
    Free(buf);

    return;
}

static void vcfwarn_free(khash_t(WARNINGS) *warnings)
{
    khiter_t key;
    for (key = kh_begin(warnings); key != kh_end(warnings); ++key)
        if (kh_exist(warnings, key))
            Free(kh_key(warnings, key));
    kh_destroy(WARNINGS, warnings);
}

struct parse_t {
    struct vcftype_t *vcf;
    struct rle_t *chrom;
    struct dna_hash_t *ref;
    khash_t(strhash) *str;      /* general purpose hash of strings */
    int vcf_n, imap_n, gmap_n, smap_n, *smapidx;
    const char **inms, **gnms, **snms;
    khash_t(WARNINGS) *warnings;
};

static struct vcftype_t *_types_alloc(const int x_n, const int y_n,
                                      Rboolean isInfo, SEXP map,
                                      khash_t(strhash) *str)
{
    SEXP elt;
    const char *n;
    SEXPTYPE type;

    struct vcftype_t *types;
    const int map_n = Rf_length(map);
    const char *dot = _strhash_put(str, ".");

    if (map_n == 0)            /* no INFO or GENO in header */
        return _vcftype_new(VECSXP, NILSXP, '\0', NULL, 0, 0, 0, 0);

    types = _vcftype_new(VECSXP, NILSXP, '\0', NULL, map_n, 1, 1, 0);

    for (int j = 0; j < map_n; ++j) {
        elt = VECTOR_ELT(map, j);
        n = CHAR(STRING_ELT(VECTOR_ELT(elt, 0), 0));
        type = TYPEOF(VECTOR_ELT(elt, 1));

        if (type == NILSXP) {   /* skip */
            types->u.list[j] =
                _vcftype_new(NILSXP, NILSXP, *n, NULL, 0, 0, 0, 0);
        } else if (*n == '.' || *n == 'A' || *n == 'G') { /* ragged array */
            types->u.list[j] =
                _vcftype_new(VECSXP, type, *n, dot, x_n, y_n, 1, 2);
        } else {                /* array */
            int z_n = atoi(n);
            int dim = (z_n == 1) ? (isInfo ? 1 : 2) : 3;
            types->u.list[j] =
                _vcftype_new(type, NILSXP, *n, dot, x_n, y_n, z_n, dim);
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

static struct vcftype_t *_vcf_alloc(const int vcf_n, SEXP smap,
                                    SEXP fmap, SEXP imap, SEXP gmap,
                                    khash_t(strhash) *str)
{
    SEXP elt;
    SEXPTYPE type;
    const char *n;
    struct vcftype_t *vcf, *rowData;

    vcf = _vcftype_new(VECSXP, NILSXP, '\0', NULL, N_FLDS, 1, 1, 0);

    /* FIXED fields */
    rowData = _vcftype_new(VECSXP, VECSXP, '\0', NULL, 2, 1, 1, 0);
    rowData->u.list[POS_IDX] =
        _vcftype_new(INTSXP, NILSXP, '\0', NULL, vcf_n, 1, 1, 0);
    rowData->u.list[ID_IDX] =
        _vcftype_new(STRSXP, NILSXP, '\0', NULL, vcf_n, 1, 1, 0);
    vcf->u.list[ROWDATA_IDX] = rowData;

    const char *blank = _strhash_put(str, ""), *dot = _strhash_put(str, ".");
    for (int i = 2; i < Rf_length(fmap); ++i) {
        const char *nm = CHAR(STRING_ELT(GET_NAMES(fmap), i));
        elt = VECTOR_ELT(fmap, i);
        n = CHAR(STRING_ELT(VECTOR_ELT(elt, 0), 0));
        type = TYPEOF(VECTOR_ELT(elt, 1));
        if (0 == strcmp(nm, "ALT")) {
            vcf->u.list[ALT_IDX] =
                _vcftype_new(VECSXP, type, *n, blank, vcf_n, 1, 1, 0);
        } else if (0 == strcmp(nm, "QUAL")) {
            vcf->u.list[QUAL_IDX] =
                _vcftype_new(type, NILSXP, *n, dot, vcf_n, 1, 1, 0);
        } else if (0 == strcmp(nm, "FILTER")) {
            vcf->u.list[FILTER_IDX] =
                _vcftype_new(type, NILSXP, *n, dot, vcf_n, 1, 1, 0);
        } else
            Rf_error("(internal) unknown 'fixed' field '%s'", nm);
    }

    /* info, geno */
    int nonzero = 0;
    for (int i = 0; i < Rf_length(smap); ++i)
        if (INTEGER(smap)[i] != 0)
            nonzero += 1;
    vcf->u.list[INFO_IDX] = _types_alloc(vcf_n, 1, TRUE, imap, str);
    vcf->u.list[GENO_IDX] = _types_alloc(vcf_n, nonzero, FALSE, gmap, str);

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
    khash_t(strhash) *str = parse->str;
    struct vcftype_t *vcf = parse->vcf, *rowData, *elt, *info;
    const int imap_n = parse->imap_n, gmap_n = parse->gmap_n;
    const int smap_n = parse->smap_n;
    const char **inms = parse->inms, **gnms = parse->gnms;
    const char **snms = parse->snms;
    const int *smapidx = parse->smapidx;
    int fmtidx, imapidx;

    int j;
    struct it_t it0, it1, it2;
    char *sample, *ifld, *ikey;

    /* fixed fields */
    char *chrom, *pos, *id, *ref, *alt, *field;
    int alt_n;

    chrom = it_init(&it0, line, '\t'); /* CHROM */
    rle_append(parse->chrom, chrom);
    rowData = vcf->u.list[ROWDATA_IDX];

    pos = it_next(&it0);      /* POS */
    rowData->u.list[POS_IDX]->u.integer[irec] = atoi(pos);

    id = it_next(&it0);       /* ID */

    ref = it_next(&it0);      /* REF */
    dna_hash_append(parse->ref, ref);

    alt = it_next(&it0);      /* ALT */
    alt_n = _vcftype_ragged_n(alt);
    _vcftype_setarray(vcf->u.list[ALT_IDX], irec, 0, alt, alt_n, str);
    _vcftype_set(vcf->u.list[QUAL_IDX], irec,
                 _strhash_put(str, it_next(&it0))); /* QUAL */
    _vcftype_set(vcf->u.list[FILTER_IDX], irec,
                 _strhash_put(str, it_next(&it0))); /* FILTER */

    if ('.' == *id && '\0' == *(id + 1)) {
        /* construct ID if missing:
           chrom\0pos\0ID\0ref\0alt\0 ==> chrom:pos_ref/alt\0 */
        *(pos - 1) = ':';
        *(id - 1) = '_';
        *(alt - 1) = '/';
        while (*ref != '\0')
            *id++ = *ref++;
        *id = '\0';
        id = chrom;
    }
    rowData->u.list[ID_IDX]->u.character[irec] = _strhash_put(str, id);

    /* INFO */
    field = it_next(&it0);
    info = vcf->u.list[INFO_IDX];
    if (1 == imap_n && NULL == inms) { /* no header; parse as char */
        elt = info->u.list[0];
        elt->u.character[irec] = _strhash_put(str, field);
    } else if (0 != imap_n) {
        for (ifld = it_init(&it1, field, ';'); '\0' != *ifld;
             ifld = it_next(&it1)) {
            ikey = it_init(&it2, ifld, '=');
            for (imapidx = 0; imapidx < imap_n; ++imapidx)
                if (0L == strcmp(ikey, inms[imapidx])) {
                    elt = info->u.list[imapidx];
                    _vcftype_setarray(elt, irec, 0, it_next(&it2), alt_n, str);
                    break;
                }
        }
        /* type 'A' with no data need to be padded w/ NA's */
        for (imapidx = 0; imapidx < imap_n; ++imapidx) {
            elt = info->u.list[imapidx];
            if (elt->number == 'A' || elt->number == 'G')
                _vcftype_padarray(elt, irec, 0, str, alt_n);
        }
    }

    /* FORMAT */
    if (0 == gmap_n)
        return;                 /* early exit */
    field = it_init(&it2, it_next(&it0), ':');
    int n_fld = it_nfld(&it2);
    int *gmapidx = Calloc(n_fld, int);
    for (fmtidx = 0; '\0' != *field; field = it_next(&it2), fmtidx++) {
        for (j = 0; j < gmap_n; ++j)
            if (0L == strcmp(field, gnms[j]))
                break;
        gmapidx[fmtidx] = j;    /* gmap_n to ignore */
    }

    /* SAMPLE */
    struct vcftype_t *geno = vcf->u.list[GENO_IDX];
    const int max_fmtidx = fmtidx;
    for (int j = 0; j < smap_n; ++j) {
        sample = it_next(&it0);
        if (0 == smapidx[j])
            continue;
        for (field = it_init(&it2, sample, ':'), fmtidx = 0;
             '\0' != *field; field = it_next(&it2), fmtidx++) {
            if (fmtidx >= max_fmtidx) {
                vcfwarn(parse->warnings,
                     "record %d sample %s: fewer FORMAT fields than GENO fields",
                     irec + 1, snms[j]);
                continue;
            }
            if (gmap_n == gmapidx[fmtidx])
                continue;    /* unknown FORMAT */
            elt = geno->u.list[ gmapidx[fmtidx] ];
            _vcftype_setarray(elt, irec, smapidx[j] - 1, field, alt_n, str);
        }
        /* type 'A' and 'G' need to be padded w/ NA's */
        for (fmtidx = 0; fmtidx < gmap_n; ++fmtidx) {
            elt = geno->u.list[fmtidx];
            if (elt->number == 'A' || elt->number == 'G')
                _vcftype_padarray(elt, irec, smapidx[j] - 1, str, alt_n);
        }
    }

    Free(gmapidx);
}

static SEXP _vcf_as_SEXP(struct parse_t *parse, SEXP fmap, SEXP smap)
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
    SET_STRING_ELT(nms, ROWDATA_IDX, mkChar("rowData"));
    SET_STRING_ELT(nms, REF_IDX, mkChar("REF"));
    SET_STRING_ELT(nms, ALT_IDX, mkChar("ALT"));
    SET_STRING_ELT(nms, QUAL_IDX, mkChar("QUAL"));
    SET_STRING_ELT(nms, FILTER_IDX, mkChar("FILTER"));
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

    SEXP samplenms;
    int nonzero = 0;
    for (int i = 0; i < Rf_length(smap); ++i)
        if (INTEGER(smap)[i] != 0)
            nonzero += 1;
    PROTECT(samplenms = Rf_allocVector(STRSXP, nonzero));
    for (int i = 0; i < parse->smap_n; ++i)
        if (INTEGER(smap)[i] != 0)
            SET_STRING_ELT(samplenms, INTEGER(smap)[i] - 1, 
                           mkChar(parse->snms[i]));
    PROTECT(nms = Rf_allocVector(VECSXP, 2));
    SET_VECTOR_ELT(nms, 0, R_NilValue);
    SET_VECTOR_ELT(nms, 1, samplenms);

    sxp = VECTOR_ELT(result, GENO_IDX);
    for (int i = 0; i < Rf_length(sxp); ++i) {
        elt = VECTOR_ELT(sxp, i);
        if (R_NilValue != elt)
            Rf_dimnamesgets(elt, nms);
    }
    UNPROTECT(2);

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

static struct parse_t *_parse_new(int vcf_n, SEXP smap, SEXP fmap,
                                  SEXP imap, SEXP gmap)
{
    struct parse_t *parse = Calloc(1, struct parse_t);

    parse->vcf_n = vcf_n;
    parse->str = _strhash_new();
    parse->vcf =
        _vcf_alloc(parse->vcf_n, smap, fmap, imap, gmap, parse->str);

    /* FIXED */
    parse->chrom = rle_new(parse->vcf_n);
    parse->ref = dna_hash_new(parse->vcf_n);

    /* INFO */
    parse->imap_n = Rf_length(imap);
    if (1 == parse->imap_n && R_NilValue == GET_NAMES(imap))
        parse->inms = NULL;
    else {
        parse->inms =
            (const char **) R_alloc(sizeof(const char *), parse->imap_n);
        for (int j = 0; j < parse->imap_n; ++j)
            parse->inms[j] = CHAR(STRING_ELT(GET_NAMES(imap), j));
    }

    /* FORMAT */
    parse->gmap_n = Rf_length(gmap);
    parse->gnms =
        (const char **) R_alloc(sizeof(const char *), parse->gmap_n);
    for (int j = 0; j < parse->gmap_n; ++j)
        parse->gnms[j] = CHAR(STRING_ELT(GET_NAMES(gmap), j));

    /* SAMPLES */
    parse->smap_n = Rf_length(smap);
    parse->snms =
        (const char **) R_alloc(sizeof(const char *), parse->smap_n);
    for (int j = 0; j < parse->smap_n; ++j)
        parse->snms[j] = CHAR(STRING_ELT(GET_NAMES(smap), j));
    parse->smapidx = INTEGER(smap);

    parse->warnings = vcfwarn_new();

    return parse;
}

static void _parse_free(struct parse_t *parse)
{
    rle_free(parse->chrom);
    dna_hash_free(parse->ref);
    vcfwarn_free(parse->warnings);
    _strhash_free(parse->str);
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

SEXP scan_vcf_connection(SEXP txt, SEXP smap, SEXP fmap, SEXP imap,
                         SEXP gmap)
{
    struct parse_t *parse;

    parse = _parse_new(Rf_length(txt), smap, fmap, imap, gmap);

    /* parse each line */
    for (int irec = 0; irec < parse->vcf_n; irec++) {
        char *line = Strdup(CHAR(STRING_ELT(txt, irec)));
        _parse(line, irec, parse);
        Free(line);
    }

    SEXP result = PROTECT(Rf_allocVector(VECSXP, 1));
    SET_VECTOR_ELT(result, 0, _vcf_as_SEXP(parse, fmap, smap));
    _vcf_types_tidy(parse, result);
    _parse_free(parse);
    UNPROTECT(1);

    return result;
}

SEXP scan_vcf_character(SEXP file, SEXP yield,
                        SEXP smap, SEXP fmap, SEXP imap, SEXP gmap)
{
    struct parse_t *parse;

    if (!IS_INTEGER(yield) || 1L != Rf_length(yield))
        Rf_error("'yield' must be integer(1)");
    if (!IS_CHARACTER(file) || 1L != Rf_length(file))
        Rf_error("'file' must be character(1) or as on ?scanVcf");

    parse = _parse_new(INTEGER(yield)[0], smap, fmap, imap, gmap);

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

        if (irec == parse->vcf_n)
            _parse_grow(parse, 0);

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
    SET_VECTOR_ELT(result, 0, _vcf_as_SEXP(parse, fmap, smap));
    _vcf_types_tidy(parse, VECTOR_ELT(result, 0));
    _parse_free(parse);
    UNPROTECT(1);

    return result;
}

SEXP tabix_as_vcf(tabix_t *tabix, ti_iter_t iter, const int yield,
                  SEXP state)
{
    const ti_conf_t *conf = ti_get_conf(tabix->idx);
    SEXP sample = VECTOR_ELT(state, 0), fmap = VECTOR_ELT(state, 1);
    const int nrec = NA_INTEGER == yield ? TBX_INIT_SIZE : yield;
    struct parse_t *parse =
        _parse_new(nrec, sample, fmap, VECTOR_ELT(state, 2),
                   VECTOR_ELT(state, 3));

    int BUFLEN = 4096;
    char *buf = Calloc(BUFLEN, char);

    int linelen;
    const char *line;

    int irec = 0;
    while (NULL != (line = ti_read(tabix, iter, &linelen))) {

        if (conf->meta_char == *line)
            continue;

        if (irec == parse->vcf_n)
            _parse_grow(parse, 0);

        if (linelen + 1 > BUFLEN) {
            Free(buf);
            BUFLEN = 2 * linelen;
            buf = Calloc(BUFLEN, char);
        }
        memcpy(buf, line, linelen);
        buf[linelen] = '\0';

        _parse(buf, irec, parse);
        irec += 1;
        if (NA_INTEGER != yield && irec == parse->vcf_n)
            break;
    }

    if (tabix->fp->errcode) {
        Free(buf);
        _parse_free(parse);
        Rf_error("read line failed, corrupt or invalid file?");
    }
        
    Free(buf);

    _vcf_grow(parse->vcf, irec);

    SEXP result = PROTECT(_vcf_as_SEXP(parse, fmap, sample));
    _vcf_types_tidy(parse, result);
    _parse_free(parse);
    UNPROTECT(1);

    return result;
}
