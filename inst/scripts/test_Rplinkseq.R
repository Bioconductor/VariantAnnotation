### --------------------------------------------------------------------------
### TESTING
### --------------------------------------------------------------------------
### August 2013
### rhino03 Ubuntu server with 387 Gb RAM

library(microbenchmark)
library(VariantAnnotation)
library(Rplinkseq)

## 494328 variants, 1092 samples, 22 INFO, 3 GENO
## ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20101123/
path <- "/home/vobencha/Download/"
fl <- paste0(path,
"ALL.chr22.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz")

### --------------------------------------------------------------------------
### Testing Rplinkseq.
### --------------------------------------------------------------------------

hdr <- scanVcfHeader(fl)
info <- rownames(info(hdr))
geno <- rownames(geno(hdr))
mask <- "reg=22:20000000..30000000" ## 133449 records

## roughly equivalent to output of scanVcf()
fun0 <- function()
{
    xx <- load.vcf(fl, mask=mask, limit=200000)
    ## Collate INFO and GENO
    lapply(info, function(i) x.consensus.meta(xx$VAR, i))
    lapply(geno, function(i) x.consensus.genotype(xx$VAR, i))
}

microbenchmark(fun0(), times=3)
#Unit: seconds
#   expr      min      lq   median       uq      max neval
# fun0() 1659.694 1666.31 1672.927 1808.175 1943.423     3

## single run memory usage (creating list of lists only)
system.time(xx <- load.vcf(fl, mask=mask, limit=200000))
#   user  system elapsed 
#798.654  46.878 848.497 
gc()
#            used   (Mb) gc trigger    (Mb)  max used   (Mb)
#Ncells 153940428 8221.4  222628720 11889.7 153941308 8221.4
#Vcells 962320366 7342.0 1050991935  8018.5 962563098 7343.8
print(object.size(xx), units="Gb")
#13.1 Gb

## Process took 16-20% of 387 Gb memory during runtime.
.20 * 387
#[1] 77.4

### --------------------------------------------------------------------------
### Testing scanVcf.
### --------------------------------------------------------------------------

tf <- TabixFile(fl)
which <- GRanges("22", IRanges(2e7, 3e7)) ## 133449 records
param <- ScanVcfParam(which=which)

fun0 <- function()
    scn <- scanVcf(tf, param=param)

microbenchmark(fun0(), times=3)
#Unit: seconds
#   expr      min       lq median       uq      max neval
# fun0() 815.7059 828.9979 842.29 1032.187 1222.084     3

## single run memory usage
system.time(scn <- scanVcf(tf, param=param))
#    user   system  elapsed 
#1193.947   26.225 1220.582 
gc()
#             used   (Mb) gc trigger    (Mb)   max used    (Mb)
#Ncells  149149655 7965.5  463692252 24763.9  441938172 23602.1
#Vcells 1026383299 7830.7 3609168658 27535.8 3074722381 23458.3
print(object.size(scn), units="Gb")
#13.1 Gb

## Process took 16-20% of 387 Gb memory during runtime.
.20 * 387
#[1] 77.4

### --------------------------------------------------------------------------
### Testing scalability of INFO fields.
### --------------------------------------------------------------------------

tf <- TabixFile(fl)
hdr <- scanVcfHeader(fl)
info <- rownames(info(hdr))

## 1 -> 22 fields 
res <- sapply(1:length(info), function(i) {
              param <- ScanVcfParam(which=which, info=info[1:i], geno=NA)
              system.time(scn <- scanVcf(tf, param=param)["elapsed"])})
res[rownames(res) == "elapsed"]
# [1] 20.428 19.053 18.659 19.054 19.764 19.724 19.613 19.565 19.739 19.704
#[11] 19.787 19.708 20.281 20.027 20.559 20.735 20.710 20.983 20.453 20.411
#[21] 20.882 21.001


### --------------------------------------------------------------------------
### Testing scalability of GENO fields.
### --------------------------------------------------------------------------

tf <- TabixFile(fl)
hdr <- scanVcfHeader(fl)
geno <- rownames(geno(hdr))

## 1 -> 3 fields 
res <- sapply(1:length(geno), function(i) {
              param <- ScanVcfParam(which=which, info=NA, geno=geno[1:i])
              system.time(scn <- scanVcf(tf, param=param)["elapsed"])})
res[rownames(res) == "elapsed"]
res[rownames(res) == "elapsed"]
#[1]  135.811  174.274 1133.909

sessionInfo()
#R version 3.0.0 Patched (2013-04-04 r62492)
#Platform: x86_64-unknown-linux-gnu (64-bit)
#
#locale:
# [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
# [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
# [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
# [7] LC_PAPER=C                 LC_NAME=C                 
# [9] LC_ADDRESS=C               LC_TELEPHONE=C            
#[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
#
#attached base packages:
#[1] parallel  stats     graphics  grDevices utils     datasets  methods  
#[8] base     
#
#other attached packages:
#[1] Rplinkseq_0.08           VariantAnnotation_1.7.38 Rsamtools_1.13.27       
#[4] Biostrings_2.29.14       GenomicRanges_1.13.35    XVector_0.1.0           
#[7] IRanges_1.19.19          BiocGenerics_0.7.3       microbenchmark_1.3-0    
#
#loaded via a namespace (and not attached):
# [1] AnnotationDbi_1.23.18   Biobase_2.21.6          biomaRt_2.17.2         
# [4] bitops_1.0-5            BSgenome_1.29.1         DBI_0.2-7              
# [7] GenomicFeatures_1.13.26 RCurl_1.95-4.1          RSQLite_0.11.4         
#[10] rtracklayer_1.21.9      stats4_3.0.0            tools_3.0.0            
#[13] XML_3.98-1.1            zlibbioc_1.7.0     
