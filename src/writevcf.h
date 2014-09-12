#ifndef _WRITEVCF_H
#define _WRITEVCF_H

#include "vcftype.h"
#include <Rdefines.h>

SEXP matrix_pasteCollapseRows(SEXP x, SEXP sep);

void make_vcf_geno(SEXP conn, SEXP fixed, SEXP format, SEXP geno, 
	               SEXP separators, SEXP vcf_dim, SEXP geno_zdim); 

#endif
