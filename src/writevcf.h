#ifndef _WRITEVCF_H
#define _WRITEVCF_H

#include <Rdefines.h>
#include "vcftype.h"

SEXP matrix_pasteCollapseRows(SEXP x, SEXP sep);

SEXP make_vcf_geno(SEXP fixed, SEXP format, SEXP geno, SEXP separators, 
                   SEXP vcf_dim, SEXP geno_zdim); 

#endif
