#ifndef _VCFFILE_H_
#define _VCFFILE_H_

#include <Rdefines.h>
#include "tabix/tabix.h"

SEXP scan_vcf_character(SEXP file, SEXP yield, SEXP sample, SEXP fmap,
                        SEXP imap, SEXP gmap, SEXP rownames);
SEXP scan_vcf_connection(SEXP txt, SEXP sample, SEXP fmap, SEXP imap,
                         SEXP gmap, SEXP rownames);
SEXP tabix_as_vcf(tabix_t *tabix, ti_iter_t iter, const int yield,
                  SEXP state, SEXP rownames);

#endif                          /* _VCFFILE_H_ */
