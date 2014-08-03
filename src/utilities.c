#include "utilities.h"

/* get_namespace */

SEXP get_namespace(const char *pkg)
{
    SEXP fun = PROTECT(findFun(install("getNamespace"), R_GlobalEnv));
    SEXP nmspc = PROTECT(mkString(pkg));
    nmspc = eval(lang2(fun, nmspc), R_GlobalEnv);
    UNPROTECT(2);
    return nmspc;
}

/* iterator to return null-terminated delimited fields */

char *it_init(struct it_t *it, char *str, char delim)
{
    it->str = str;
    it->delim = delim;
    it->n_fld = (*str == '\0') ? 0 : 1;
    while (*str != '\0')
        if (*str++ == delim) it->n_fld += 1;

    return it_next(it);
}

inline char *it_next(struct it_t *it)
{
    const char delim = it->delim, *curr = it->str;
    char *start = it->str;

    while ('\0' != *curr && delim != *curr)
        ++curr;

    it->str += curr - it->str;
    if ('\0' != *it->str)
        *it->str++ = '\0';
    return start;
}

inline int it_nfld(const struct it_t *it)
{
    return it->n_fld;
}

/* paste-collapse character matrix by row, ignoring NAs */

SEXP matrix_pasteCollapseRows(SEXP x, SEXP sep) {
  int nc = ncols(x), nr = nrows(x);
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
