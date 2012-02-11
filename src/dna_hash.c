#include "dna_hash.h"

/* DNAStringSet -- hash */

#include "Biostrings_interface.h"
#include "IRanges_interface.h"
#include "samtools/khash.h"

KHASH_MAP_INIT_STR(ref, int);

struct dna_hash_t {
    khash_t(ref) *hash;
    int len, size, hash_idx;
    SEXP offset;
};

static const double DNA_GROW = 1.6;

struct dna_hash_t *dna_hash_new(const int size)
{
    struct dna_hash_t *dna = Calloc(1, struct dna_hash_t);
    dna->hash = kh_init(ref);
    dna->offset = Rf_allocVector(INTSXP, size);
    dna->size = size;
    dna->len = dna->hash_idx = 0;
    R_PreserveObject(dna->offset);
    return dna;
}

void dna_hash_free(struct dna_hash_t *dna)
{
    khiter_t key;
    for (key = kh_begin(dna->hash); key != kh_end(dna->hash); ++key) {
        if (kh_exist(dna->hash, key))
            Free(kh_key(dna->hash, key));
    }
    kh_destroy(ref, dna->hash);
    R_ReleaseObject(dna->offset);
    Free(dna);
}

void dna_hash_grow(struct dna_hash_t *dna, int size)
{
    SEXP updt = Rf_lengthgets(dna->offset, size);
    R_PreserveObject(updt);
    R_ReleaseObject(dna->offset);
    dna->offset = updt;
    dna->size = size;
}

void dna_hash_append(struct dna_hash_t *dna, const char *value)
{
    khiter_t key;
    key = kh_get(ref, dna->hash, value);
    if (key == kh_end(dna->hash)) {
        int ret;
        char *buf = Calloc(strlen(value) + 1, char);
        strcpy(buf, value);
        key = kh_put(ref, dna->hash, buf, &ret);
        kh_value(dna->hash, key) = dna->hash_idx++;
    }

    if (dna->len == dna->size)
        dna_hash_grow(dna, DNA_GROW * dna->size);
    INTEGER(dna->offset)[dna->len] = kh_value(dna->hash, key);
    dna->len++;
}

SEXP dna_hash_as_DNAStringSet(struct dna_hash_t *dna)
{
    int *iwidth, *istart, twidth, i;
    SEXP tag, start, width, ranges, xstringset;
    Rbyte *tagp;
    khiter_t key;

    istart = Calloc(dna->hash_idx, int);
    iwidth = Calloc(dna->hash_idx, int);

    for (key = kh_begin(dna->hash); key != kh_end(dna->hash); ++key)
    {
        if (!kh_exist(dna->hash, key))
            continue;
        kh_cstr_t cstr = kh_key(dna->hash, key);
        int idx = kh_value(dna->hash, key);
        iwidth[idx] = strlen(cstr);
    }

    istart[0] = 1;
    for (i = 1; i < dna->hash_idx; ++i)
        istart[i] = istart[i-1] + iwidth[i - 1];
    twidth = istart[dna->hash_idx -1 ] + iwidth[dna->hash_idx - 1];

    /* RAW */
    PROTECT(tag = NEW_RAW(twidth)); tagp = RAW(tag);
    for (key = kh_begin(dna->hash); key != kh_end(dna->hash); ++key)
    {
        if (!kh_exist(dna->hash, key))
            continue;
        kh_cstr_t cstr = kh_key(dna->hash, key);
        int idx = kh_value(dna->hash, key);
        for (int j = 0; j < iwidth[idx]; ++j) {
            *tagp++ = DNAencode(cstr[j]);
        }
    }

    /* ranges */
    PROTECT(start = NEW_INTEGER(dna->len));
    PROTECT(width = NEW_INTEGER(dna->len));
    for (i = 0; i < dna->len; ++i) {
        int idx = INTEGER(dna->offset)[i];
        INTEGER(start)[i] = istart[idx];
        INTEGER(width)[i] = iwidth[idx];
    }
    PROTECT(ranges = new_IRanges("IRanges", start, width, R_NilValue));

    /* DNAStringSet */
    PROTECT(xstringset = new_XRawList_from_tag(
                "DNAStringSet", "DNAString", tag, ranges));

    free(iwidth);
    free(istart);
    UNPROTECT(5);

    return xstringset;
}
