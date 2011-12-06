#include <R_ext/Rdynload.h>
#include "writevcf.h"


static const R_CallMethodDef callMethods[] = {
    /* writevcf.c */
    {".code_allele_observations", (DL_FUNC) & code_allele_observations, 2},
    {NULL, NULL, 0}
};

void R_init_VariantAnnotation(DllInfo * info)
{
    R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}

