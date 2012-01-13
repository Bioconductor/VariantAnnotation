#include <R_ext/Rdynload.h>
#include "writevcf.h"
#include "vcffile.h"


static const R_CallMethodDef callMethods[] = {
    /* vcffile.c */
    {".scan_vcf", (DL_FUNC) & scan_vcf, 4},
    {".scan_vcf_connection", (DL_FUNC) & scan_vcf_connection, 4},
    /* writevcf.c */
    {".code_allele_observations", (DL_FUNC) & code_allele_observations, 2},
    {NULL, NULL, 0}
};

void R_init_VariantAnnotation(DllInfo * info)
{
    R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}

