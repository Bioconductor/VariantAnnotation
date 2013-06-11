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
