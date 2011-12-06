#include <R.h>
#include <Rinternals.h>
#include "writevcf.h"

SEXP code_allele_observations(SEXP options, SEXP observed) {
  
  /* options and observed are both vectors of character arrays where 
     each character array is a collection of alleles at that position */

  SEXP obs_string, opt_string;
  const char *obs_char, *opt_char;
  SEXP coded = PROTECT(allocVector(STRSXP, length(options)));
  int max_num_alleles = 10;  /* World will end if index more than one character long */
  int buffer_index = 0;
  char coded_buffer[max_num_alleles * 2]; /* One int char and one comma for each observation at a position */
  
  for (int i=0; i < length(options); i++) {
    opt_string   = VECTOR_ELT(options, i);
    obs_string   = VECTOR_ELT(observed, i);
    if (length(obs_string) > 10) { error("No more than 10 alleles allowed\n"); }
    buffer_index = 0;
    for (int observed_index=0; observed_index < length(obs_string); observed_index++) {
      coded_buffer[ buffer_index ] = '.';  /* Set to NA to begin with */
      coded_buffer[buffer_index + 1] = '/';
      obs_char = CHAR(STRING_ELT(obs_string,observed_index));
    
      for (int option_index=0; option_index < length(opt_string); option_index++) {
	opt_char = CHAR(STRING_ELT(opt_string,option_index));
	if (obs_char == opt_char) {
	  sprintf(coded_buffer + buffer_index, "%i/", option_index);
	  break;
	}
      }
      buffer_index += 2;
    }
    SET_STRING_ELT(coded, i, mkCharLenCE(coded_buffer, buffer_index-1, CE_UTF8));
  }
  UNPROTECT(1);
  return(coded);
}
