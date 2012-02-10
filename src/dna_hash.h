#ifndef _DNA_HASH_H
#define _DNA_HASH_H

#include <Rdefines.h>

struct dna_hash_t *dna_hash_new(const int size);
void dna_hash_free(struct dna_hash_t *dna);
void dna_hash_grow(struct dna_hash_t *dna, int size);
void dna_hash_append(struct dna_hash_t *dna, const char *value);
SEXP dna_hash_as_DNAStringSet(struct dna_hash_t *dna);

#endif
