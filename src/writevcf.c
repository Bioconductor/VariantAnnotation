#include "writevcf.h"
#include "samtools/kstring.h"

/* write all genotype fields for a single sample */
void write_geno_sample(int i, int j, int k_last, SEXP k_valid, SEXP geno_zdim,
                       int n_rows, int n_samples, int n_fields, SEXP geno, 
                       const char *f_sep, const char *mv_sep, kstring_t *bufp)
{
    SEXP field, c_elt;
    SEXPTYPE type;
    double d_elt;
    int k, z, z_dim, z_max, index, i_elt;
 
    for (k = 0; k < n_fields; ++k) {
        field = VECTOR_ELT(geno, k);
        type = TYPEOF(field);
        z_dim = INTEGER(geno_zdim)[k];
        z_max = (NA_INTEGER != z_dim) ? z_dim : 1;
        if (LOGICAL(k_valid)[k]) {
            for (z = 0; z < z_max; ++z) {
                index = i + j*n_rows + z*n_rows*n_samples; 
                switch (type) {
                    case NILSXP:
                        break;
                    case LGLSXP:
                        Rf_warning("'logical' is not a valid FORMAT data type");
                        break;
                    case INTSXP:
                        i_elt = INTEGER(field)[index];
                        if (NA_INTEGER != i_elt) 
                            ksprintf(bufp, "%d", i_elt);
                        else
                            ksprintf(bufp, "%s", ".");
                        break;
                    case REALSXP:
                        d_elt = REAL(field)[index];
                        if (NA_REAL != d_elt)
                            ksprintf(bufp, "%f", d_elt);
                        else
                            ksprintf(bufp, "%s", ".");
                        break;
                    case STRSXP:
                        c_elt = STRING_ELT(field, index);
                        if (NA_STRING != c_elt)
                            ksprintf(bufp, "%s", CHAR(c_elt));
                        else
                            ksprintf(bufp, "%s", ".");
                        break;
                }
                /* multi-value separator */
                if (z < z_max - 1 && LOGICAL(k_valid)[k])
                    ksprintf(bufp, "%s", mv_sep);
            }
            /* field separator */
            if (k < n_fields - 1 && k < k_last)
                ksprintf(bufp, "%s", f_sep);
        }
    }
}

/* return TRUE if genotype element is NA */
Rboolean valid_geno_elt(SEXP field, SEXPTYPE type, int index)
{
    Rboolean valid = TRUE;

    switch (type) {
        case NILSXP:
            break;
        case LGLSXP:
            if (NA_INTEGER == LOGICAL(field)[index])
                valid = FALSE;
            break;
        case INTSXP:
            if (NA_INTEGER == INTEGER(field)[index]) 
                valid = FALSE;
            break;
        case REALSXP:
            if (NA_REAL == REAL(field)[index])
                valid = FALSE;
            break;
        case STRSXP:
            if (NA_STRING == STRING_ELT(field, index))
                valid = FALSE;
            break;
    }
    return valid;
} 

/* --- .Call ENTRY POINT --- 
 * 'format'     : character vector of FORMAT names
 * 'geno'       : list of genotype data
 * 'separators' : character vector of length 2 consisting of a field
 *                separator (first) and multi-value per field 
 *                separator (second)
 * 'vcf_dim'    : integer vector of length 2 (n rows, n cols)
 * 'geno_zdim'  : integer vector of z dimension of genotype data 
*/
SEXP make_vcf_geno(SEXP format, SEXP geno, SEXP separators,
                   SEXP vcf_dim, SEXP geno_zdim)
{
    const char *f_sep = CHAR(STRING_ELT(separators, 0));
    const char *mv_sep = CHAR(STRING_ELT(separators, 1));
    int n_rows = INTEGER(vcf_dim)[0];
    int n_samples = INTEGER(vcf_dim)[1];
    int n_fields = length(format);
    int i, j, k, k_last, z, z_dim, z_max;
    int index;
    Rboolean search;

    kstring_t buf;
    buf.l = buf.m = 0; buf.s = NULL;

    SEXP k_valid, field, result; 
    SEXPTYPE type; 

    if (n_fields != length(geno))
      Rf_error("length(format) must equal length(geno)");
    if (length(geno_zdim) != length(geno))
      Rf_error("length(geno_zdim) must equal length(geno)");

    PROTECT(k_valid = NEW_LOGICAL(n_fields));
    PROTECT(result = allocVector(STRSXP, n_rows));

    for (i = 0; i < n_rows; ++i) {
        /* reset buffer */
        buf.l = 0;
        if (NULL != buf.s) *buf.s = '\0';

        k_last = 0;
        /* write format names */
        for (k = 0; k < n_fields; ++k) {
            search = TRUE;
            field = VECTOR_ELT(geno, k);
            type = TYPEOF(field);
            z_dim = INTEGER(geno_zdim)[k];
            z_max = (NA_INTEGER != z_dim) ? z_dim : 1;
            for (j = 0; j < n_samples && search; ++j) {
                for (z = 0; z < z_max; ++z) {
                    index = i + j*n_rows + z*n_rows*n_samples; 
                    if (valid_geno_elt(field, type, index)) {
                        ksprintf(&buf, "%s%s", CHAR(STRING_ELT(format, k)),
                                 (k < n_fields - 1) ? f_sep : "\t");
                        LOGICAL(k_valid)[k] = TRUE;
                        k_last = k;
                        search = FALSE;
                        break;
                    } else if (k == n_fields - 1 && z == z_max - 1) {
                        ksprintf(&buf, "%s", "\t");
                        LOGICAL(k_valid)[k] = FALSE;
                        search = FALSE;
                        break;
                    } else {
                        LOGICAL(k_valid)[k] = FALSE;
                    }
                }
            }
        }
        /* write genotype fields */
        for (j = 0; j < n_samples; ++j) {
            write_geno_sample(i, j, k_last, k_valid, geno_zdim, n_rows, 
                              n_samples, n_fields, geno, f_sep, mv_sep, 
                              &buf);
            if (j < n_samples - 1) {
                if (0 != buf.l)
                    /* sample separator */
                    ksprintf(&buf, "%s", "\t");
            } else {
                SET_STRING_ELT(result, i, mkChar(buf.s));
            }
        }
    }

    free(buf.s);
    UNPROTECT(2);
    return result;
}

/* paste-collapse character matrix by row, ignoring NAs */
SEXP matrix_pasteCollapseRows(SEXP x, SEXP sep) {
  /* VO: looks like 'nc' is not used?
   * int nc = ncols(x), nr = nrows(x);
   */
  int nr = nrows(x);
  char c_sep = CHAR(STRING_ELT(sep, 0))[0];

  SEXP ans = allocVector(STRSXP, nr);
  PROTECT(ans);

  for (int r = 0; r < nr; r++) {
    int len = 0;
    for (int i = r; i < length(x); i += nr) {
      SEXP str = STRING_ELT(x, i);
      if (str != NA_STRING)
        len += length(str) + 1;
    }
    char *collapsed = R_alloc(sizeof(char), len);
    char *dest = collapsed;
    for (int i = r; i < length(x); i += nr) {
      SEXP str = STRING_ELT(x, i);
      if (str != NA_STRING) {
        strcpy(dest, CHAR(str));
        dest[length(str)] = c_sep;
        dest += length(str) + 1;
      }
    }
    SET_STRING_ELT(ans, r, mkCharLen(collapsed, len - 1 * (len > 0)));
  }

  UNPROTECT(1);
  return ans;
}
