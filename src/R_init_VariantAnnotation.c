#include <R_ext/Rdynload.h>
#include "vcffile.h"

static const R_CallMethodDef callMethods[] = {
    /* vcffile.c */
    {".scan_vcf_character", (DL_FUNC) & scan_vcf_character, 6},
    {".scan_vcf_connection", (DL_FUNC) & scan_vcf_connection, 5},
    {".tabix_as_vcf", (DL_FUNC) & tabix_as_vcf, 5},
    {NULL, NULL, 0}
};

void R_init_VariantAnnotation(DllInfo * info)
{
    R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}
