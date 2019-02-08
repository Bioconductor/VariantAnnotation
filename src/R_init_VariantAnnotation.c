#include <R_ext/Rdynload.h>
#include "vcffile.h"
#include "writevcf.h"
#include "utilities.h"

static const R_CallMethodDef callMethods[] = {
    /* vcffile.c */
    {".scan_vcf_character", (DL_FUNC) & scan_vcf_character, 7},
    {".scan_vcf_connection", (DL_FUNC) & scan_vcf_connection, 6},
    {".tabix_as_vcf", (DL_FUNC) & tabix_as_vcf, 6},
    {"matrix_pasteCollapseRows", (DL_FUNC) & matrix_pasteCollapseRows, 2},
    {".make_vcf_geno", (DL_FUNC) & make_vcf_geno, 7},
    {NULL, NULL, 0}
};

void R_init_VariantAnnotation(DllInfo * info)
{
    R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}
