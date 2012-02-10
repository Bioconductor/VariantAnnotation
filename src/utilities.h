#ifndef _UTILITIES_H
#define _UTILITIES_H

#include <Rdefines.h>

SEXP get_namespace(const char *pkg);

struct it_t {
    char *str;
    char delim;
};

char *it_init(struct it_t *it, char *str, char delim);
static inline char *it_next(struct it_t *it)
{
    char *curr = it->str;
    while ('\0' != *it->str && it->delim != *it->str)
        it->str++;
    if ('\0' != *it->str)
        *it->str++ = '\0';
    return curr;
}


#endif
