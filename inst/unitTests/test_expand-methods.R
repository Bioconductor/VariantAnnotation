expand <- VariantAnnotation::expand 

test_expand_info_geno <- function()
{
    fl <- system.file("unitTests", "cases", "expand.vcf", 
        package="VariantAnnotation")
    vcf <- suppressWarnings(readVcf(fl, "hg19"))
    exp <- expand(vcf)
    names(mcols(rowRanges(exp)))
    cnms <- c("paramRangeID", "REF", "ALT", "QUAL", "FILTER")
    checkIdentical(cnms, names(mcols(rowRanges(exp))))
    checkTrue(nrow(exp) == 22L)
    checkIdentical(info(exp)$DP[5:6], c(10L, 10L))
    checkIdentical(info(exp)$AF[5:6], c(0.333, 0.667))
    checkIdentical(geno(exp)$AAP[6:11], c(2L, rep(1L, 4L), 3L))
}

test_expand_structural <- function()
{
    fl <- system.file("extdata", "structural.vcf", package="VariantAnnotation")
    vcf <- readVcf(fl, "hg19")
    exp <- expand(vcf)
    checkTrue(nrow(vcf) == nrow(exp))
    checkTrue(ncol(vcf) == ncol(exp))
    checkTrue(is.character(alt(exp)))
}

test_expand_structural <- function()
{
    fl <- system.file("unitTests", "cases", "mixedStructural.vcf", 
                      package="VariantAnnotation")
    vcf <- readVcf(fl, "hg19")
    exp <- expand(vcf)
    checkTrue(nrow(exp) == 8L)
}

test_expand_multiple_info <- function()
{
    fl <- system.file("unitTests", "cases", "multiple_INFO_fields.vcf", 
                      package="VariantAnnotation")
    vcf <- readVcf(fl, "GRCh37")
    exp <- expand(vcf)
    checkTrue(nrow(exp) == 3L)

    ## single row, subset of INFO fields selected 
    exp <- expand(vcf[1,])
    checkTrue(nrow(exp) == 2L)

    ## single row, single (only) INFO field selected
    vcf <- readVcf(fl, "GRCh37", param=ScanVcfParam(info="AF"))
    exp <- expand(vcf[1,])
    checkTrue(nrow(exp) == 2L)
}

test_expand_gvcf <- function()
{
    fl <- system.file("unitTests", "cases", "banded_gvcf.vcf",
                      package="VariantAnnotation")
    vcf <- suppressWarnings(readVcf(fl, ""))
    exp <- expand(vcf)
    checkIdentical(dim(exp), c(7L, 1L))
    checkIdentical(dim(geno(exp)$AD), c(7L, 1L, 2L))
    checkIdentical(geno(exp)$AD[,,1], as.integer(c(NA, 17, 17, NA, 20, 20, NA)))
}



test_expand_adr_adf <- function()
{
	fl <- system.file("unitTests", "cases", "ex1-seq1-90.vcf", package = "VariantAnnotation")
	vcf <- readVcf(fl, genome = "")
	exp <- expand(vcf)
	seqName <- "seq1:90_N/G";sampleName <- colnames(vcf)[1]
	
	# for all/forward/revert strands
	tmp <- sapply(c("AD", "ADF", "ADR"), function(var){
				
				# original counts
				countsVcf <- geno(vcf)[[var]][seqName, sampleName][[1]]
				
				# in VCF: counts reported first for reference, then alternative alleles
				# check counts for reference allele
				checkTrue(all(geno(exp)[[var]][, sampleName, 1] == countsVcf[1]))
				
				# check if counts for alternative alleles
				idxMatch <- match(alt(vcf)[[1]], alt(exp))
				countsExpand <- geno(exp)[[var]][idxMatch, sampleName, 2]
				checkIdentical(countsExpand, countsVcf[-1]) 
				
			})
	
}
