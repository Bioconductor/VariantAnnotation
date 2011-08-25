#include <R_ext/Rdynload.h>
#include "vcf_split.h"

static const R_CallMethodDef callMethods[] = {
    {".vcf_split", (DL_FUNC) &vcf_split, 3},
    {NULL, NULL, 0}
};

void
R_init_VariantAnnotation(DllInfo *info)
{
    R_registerRoutines(info, NULL, callMethods, NULL, NULL );
}
