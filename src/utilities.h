#ifndef _UTILITIES_H
#define _UTILITIES_H

#include <Rdefines.h>

#define Strdup(x) strcpy(Calloc(strlen(x) + 1, char), x)

SEXP get_namespace(const char *pkg);

struct it_t {
    char *str;
    char delim;
    int n_fld;
};

char *it_init(struct it_t *it, char *str, char delim);
char *it_next(struct it_t *it);
int it_nfld(const struct it_t *it);

SEXP matrix_pasteCollapseRows(SEXP x, SEXP sep);
#endif
