test_vcfFields <- function(){
    ## invalid
    checkException(vcfFields(NA_character_), silent = TRUE)
    checkException(vcfFields(tempfile()), silent = TRUE)

    ## empty
    target <- CharacterList(
        fixed = character(), info = character(), geno = character(),
        samples = character()
    )
    checkIdentical(target, vcfFields())

    ## structure
    fl <- system.file("extdata", "ex2.vcf", package="VariantAnnotation")
    flds <- vcfFields(fl)

    checkTrue(validObject(flds))
    checkTrue(is(flds, "CharacterList"))
    target <- c(fixed = 4L, info = 6L, geno = 4L, samples = 3L)
    checkIdentical(target, lengths(flds))

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

    ## VCFFileList
    vfl <- VcfFileList(c(fl, fl))
    checkIdentical(vcfFields(vfl[[1]]), vcfFields(vfl))
}
