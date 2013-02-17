library(VariantAnnotation)
library(RUnit)
#-------------------------------------------------------------------------------
CGI_DIR <- "/shared/silo_researcher/Morgan_M/BioC/vcf/completeGenomics"

      SMALL_VCF_FILE <- file.path(CGI_DIR, "h1187-10k.vcf.gz")
     MEDIUM_VCF_FILE <- file.path(CGI_DIR, "h1187-750k.vcf.gz")
 COMPLETE_VCF_FILE_1 <- file.path(CGI_DIR, 
                           "somaticVcfBeta-HCC1187-H-200-37-ASM-T1-N1.vcf.gz")
 COMPLETE_VCF_FILE_2 <- file.path(CGI_DIR, 
                           "somaticVcfBeta-HCC2218-H-200-37-ASM-T1-N1.vcf.gz")

stopifnot(file.exists(SMALL_VCF_FILE))
stopifnot(file.exists(MEDIUM_VCF_FILE))
stopifnot(file.exists(COMPLETE_VCF_FILE_1))
stopifnot(file.exists(COMPLETE_VCF_FILE_2))
#-------------------------------------------------------------------------------
runTests <- function()
{
    #test_allelicDepthFilter()
    #test_germlineFilter()
    test_germlinePrefilter()
    test_dbsnpPrefilter()
    #test_snpFilter()
    #test_combinedFilters()
    
} # runTests
#-------------------------------------------------------------------------------
isGermlinePrefilter <- function(x) {
    grepl("Germline", x, fixed=TRUE)
    }

notInDbsnpPrefilter <- function(x) {
    !(grepl("dbsnp", x, fixed=TRUE))
    }

isGermlineFilter <- function(x) {
    result <-info(x)$SS=="Germline"
    result[is.na(result)] <- FALSE
    result
    }

isSnp <- function(x) {
    refSnp <- nchar(ref(x)) == 1L
    a <- alt(x)
    altSnp <- elementLengths(a) == 1L
    ai <- unlist(a[altSnp])    # all length 1, so unlisting is 1:1 map
    altSnp[altSnp] <- nchar(ai) == 1L & (ai %in% c("A", "C", "G", "T"))
    refSnp & altSnp
    }

allelicDepth <- function(x) {
    ##  ratio of AD of the 'alternate allele' for the tumor sample
    ##  OR 'reference allele' for normal samples to total reads for
    ##  the sample should be greater than some threshold (say 0.1,
    ##  that is: at least 10% of the sample should have the allele
    ##  of interest)
    ad <- geno(x)$AD
    tumorPct <- ad[,1,2,drop=FALSE] / rowSums(ad[,1,,drop=FALSE])
    normPct <- ad[,2,1, drop=FALSE] / rowSums(ad[,2,,drop=FALSE])
    test <- (tumorPct > 0.1) | (normPct > 0.1)
    !is.na(test) & test
    }



#-------------------------------------------------------------------------------
recordTiming <- function(tbl, filename, filterName, elapsedTime)
{

   rbind(tbl, data.frame(file=basename(filename), cmd=filterName,
                          elapsed=elapsedTime))

} # recordTiming
#-------------------------------------------------------------------------------
runGermlineFilter <- function (gzippedFilename, tbl, yieldSize=10000, 
                               verbose=FALSE)
{
    tabix.file <- TabixFile(gzippedFilename, yieldSize=yieldSize)
    filters <- FilterRules(list(germline=isGermlineFilter))
    destination.file <- tempfile()
    filtered.filename <- filterVcf(gzippedFilename, "hg19", destination.file, 
                                   filters=filters, verbose=verbose)
    filtered.filename

} # runGermlineFilter
#-------------------------------------------------------------------------------
runGermlinePrefilter <- function (gzippedFilename, tbl, yieldSize=10000, 
                                  verbose=FALSE)
{
    tabix.file <- TabixFile(gzippedFilename, yieldSize=yieldSize)
    filters <- FilterRules(list(germline=isGermlinePrefilter))
    destination.file <- tempfile()
    filtered.filename <- filterVcf(gzippedFilename, "hg19", destination.file, 
                                   prefilters=filters, verbose=verbose)
    filtered.filename

} # runGermlinePrefilter
#-------------------------------------------------------------------------------
runDbsnpPrefilter <- function (gzippedFilename, tbl, yieldSize=10000, 
                               verbose=FALSE)
{
    tabix.file <- TabixFile(gzippedFilename, yieldSize=yieldSize)
    filters <- FilterRules(list(germline=notInDbsnpPrefilter))
    destination.file <- tempfile()
    filtered.filename <- filterVcf(gzippedFilename, "hg19", destination.file, 
                                   prefilters=filters, verbose=verbose)
    filtered.filename

} # runDbsnpPrefilter
#-------------------------------------------------------------------------------
runSnpFilter <- function (gzippedFilename, yieldSize=10000, verbose=FALSE)
{
    tabix.file <- TabixFile(gzippedFilename, yieldSize=yieldSize)
    filteringFunctions <- FilterRules(list(isSnp=isSnp))
    destination.file <- tempfile()
    filtered.filename <- filterVcf(tabix.file, "hg19", destination.file,
                                   filters=filteringFunctions,
                                   verbose=verbose)
    filtered.filename

} # runSnpFilter
#-------------------------------------------------------------------------------
runAllelicDepthFilter <- function (gzippedFilename, yieldSize=10000, 
                                    verbose=FALSE)
{
    tabix.file <- TabixFile(gzippedFilename, yieldSize=yieldSize)
    filteringFunctions <- FilterRules(list(AD=allelicDepth))
    destination.file <- tempfile()
    filtered.filename <- filterVcf(tabix.file, "hg19", destination.file,
                                   filters=filteringFunctions,
                                   verbose=verbose)
    filtered.filename

} # runSnpFilter
#-------------------------------------------------------------------------------
runCombinedFilters <- function (gzippedFilename, yieldSize=10000, verbose=FALSE)
{
    tabix.file <- TabixFile(gzippedFilename, yieldSize=yieldSize)
    prefilters <- FilterRules(list(germline=isGermlinePrefilter,
                              dbsnp=notInDbsnpPrefilter))
    filters <- FilterRules(list(isSnp=isSnp, AD=allelicDepth))
    destination.file <- tempfile()
    filtered.filename <- filterVcf(tabix.file, "hg19", destination.file,
                                   prefilters=prefilters,
                                   filters=filters,
                                   verbose=verbose)
    filtered.filename

} # runSnpFilter
#-------------------------------------------------------------------------------
# use the SMALL_VCF_FILE to ensure that filtering (not prefiltering) by 
#   somatic status == Germline
# produces the expected count
test_germlineFilter <- function()
{
    print("---  test_germlineFilter")
    vcf <- readVcf(SMALL_VCF_FILE, "hg19")
    somaticStatus <- as.list(table(info(vcf)$SS))
    checkEquals(somaticStatus, list(Germline=103, LOH=1, Somatic=4))

        # filter the in-memory data structure
    filters <- FilterRules(list(filter.1=isGermlineFilter))
    checkEquals(dim(subsetByFilter(vcf, filters)), c(103, 2))

        # use another filtering capability: create a new filtered file
    filtered.filename <- runGermlineFilter(SMALL_VCF_FILE)
    return(checkEquals(nrow(readVcf(filtered.filename, "hg19")), 103))

} # test_germlineFilter
#-------------------------------------------------------------------------------
test_germlinePrefilter <- function()
{
    print("--- test_germlinePrefilter")
    filtered.filename <- runGermlinePrefilter(SMALL_VCF_FILE)
    #tabix.file <- TabixFile(SMALL_VCF_FILE, yieldSize=1000)

    #prefilteringFunctions <- FilterRules(list(isGermline=isGermlinePrefilter))
    #destination.file <- tempfile()
    #filtered.filename <- filterVcf(tabix.file, "hg19", destination.file,
    #                               prefilters=prefilteringFunctions,
    #                               verbose=FALSE)

    checkEquals(nrow(readVcf(filtered.filename, "hg19")), 103)

} # test_germlinePrefilter
#-------------------------------------------------------------------------------
test_dbsnpPrefilter <- function()
{
    print("--- test_dbsnpPrefilter")
    filtered.filename <- runDbsnpPrefilter(SMALL_VCF_FILE)
    checkEquals(nrow(readVcf(filtered.filename, "hg19")), 9833) # 166)

} # test_germlinePrefilter
#-------------------------------------------------------------------------------
# of the 9999 rows in the SMALL_VCF_FILE, 186 describe snps
test_snpFilter <- function()
{
    print("--- test_snpFilter")
    filtered.filename <- runSnpFilter(SMALL_VCF_FILE)
    checkEquals(nrow(readVcf(filtered.filename, "hg19")), 186)

} # test_snpFilter
#-------------------------------------------------------------------------------
# of the 9999 rows in the SMALL_VCF_FILE, 186 describe snps
test_allelicDepthFilter <- function()
{
    print("--- test_allelicDepthFilter")
    filtered.filename <- runAllelicDepthFilter(SMALL_VCF_FILE)
    checkEquals(nrow(readVcf(filtered.filename, "hg19")), 160)

} # test_allelicDepthFilter
#-------------------------------------------------------------------------------
# read the 10k gzipped complete genomics vcf file, apply isGermline prefilter and
# isSnp filter, finding rows which meet both criteria, which should be 98
test_combinedFilters <- function()
{
    print("--- test_combinedFilters")
    filtered.filename <- runCombinedFilters(SMALL_VCF_FILE)
        # now check the results
    return(checkEquals(nrow(readVcf(filtered.filename, "hg19")), 96))

} # test_combinedFilters
#-------------------------------------------------------------------------------
collectTimings <- function()
{
    print("--- collectTimings")

    prefilteringFunctions <- FilterRules(list(isGermline=isGermlinePrefilter))
    filteringFunctions <- FilterRules(list(isSnp=isSnp))

    files <- c(SMALL_VCF_FILE, MEDIUM_VCF_FILE, COMPLETE_VCF_FILE_1)
    options(stringsAsFactors=FALSE)
    tbl <- data.frame()

    for(file in files[3]) {
        time <- system.time(runGermlinePrefilter(file, tbl))
        tbl <- recordTiming(tbl, file, "germline prefilter", time[["elapsed"]])
        print(tbl)

        time <- system.time(runDbsnpPrefilter(file, tbl))
        tbl <- recordTiming(tbl, file, "dbsnp prefilter", time[["elapsed"]])
        print(tbl)

        #time <- system.time(runGermlineFilter(file, tbl))
        #tbl <- recordTiming(tbl, file, "germline filter", time[["elapsed"]])
        #print(tbl)

        time <- system.time(runSnpFilter(file))
        tbl <- recordTiming(tbl, file, "snp filter", time[["elapsed"]])
        print(tbl)

        time <- system.time(runAllelicDepthFilter(file))
        tbl <- recordTiming(tbl, file, "allelic depth filter", time[["elapsed"]])
        print(tbl)

        time <- system.time(runCombinedFilters(file))
        tbl <- recordTiming(tbl, file, "combined filters", time[["elapsed"]])
        print(tbl)

        filename <- sprintf("tbl-%s.tsv", format (Sys.time(), "%a.%b.%d.%Y-%H:%M:%S"))
        write.table (tbl, file=filename, sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)
        } # for

} # collectTimings
#-------------------------------------------------------------------------------
#filtered.filename <- runCombinedFilters(SMALL_VCF_FILE)
#filtered.filename <- runCombinedFilters(COMPLETE_VCF_FILE_1)
#print(filtered.filename)
#cmd = sprintf("cp %s %s", filtered.filename, "filtered.vcf")
#system(cmd)
# ps -aux | grep  16098
