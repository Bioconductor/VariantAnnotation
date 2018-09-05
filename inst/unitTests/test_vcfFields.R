test_vcfFields <- function(){

    ## empty
    checkException(vcfFields(), silent = TRUE)
    checkException(vcfFields(NA), silent = TRUE)
    checkException(vcfFields("/path/to/file"), silent = TRUE)

    ## structure
    fl <- system.file("extdata", "ex2.vcf", package="VariantAnnotation")
    flds <- vcfFields(fl)

    checkTrue(validObject(flds))
    checkTrue(is(flds, "CharacterList"))
    checkIdentical(names(flds), c("fixed", "info", "geno", "samples"))

    ## signatures
    hdr <- scanVcfHeader(fl)
    flds.hdr <- vcfFields(hdr)

    vf <- VcfFile(fl)
    flds.vf <- vcfFields(vf)

    vcf <- readVcf(fl, genome = "hg19")
    flds.vcf <- vcfFields(vcf)
    
    checkIdentical(flds, flds.hdr)
    checkIdentical(flds, flds.vf)
    checkIdentical(flds, flds.vcf)
    
}
