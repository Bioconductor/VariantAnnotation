### --------------------------------------------------------------------------
### Testing Rplinkseq and VariantAnnotation
### --------------------------------------------------------------------------
###
### September 2013
### 'rhino03' Ubuntu server, 387 Gb RAM
### 16 processors with the following configuration:
###
### vendor_id	: GenuineIntel
### cpu family	: 6
### model		: 45
### model name	: Intel(R) Xeon(R) CPU E5-2690 0 @ 2.90GHz
### stepping	: 7
### microcode	: 0x70d
### cpu MHz		: 2900.142
### cache size	: 20480 KB
### physical id	: 1
### siblings	: 8
### core id		: 0
### cpu cores	: 8
### apicid		: 32
### initial apicid	: 32
### fpu		: yes
### fpu_exception	: yes
### cpuid level	: 13
### wp		: yes
### flags		: fpu vme de pse tsc msr pae mce cx8 apic sep mtrr pge mca cmov
### pat pse36 clflush dts acpi mmx fxsr sse sse2 ss ht tm pbe syscall nx pdpe1gb
### rdtscp lm constant_tsc arch_perfmon pebs bts rep_good nopl xtopology nonstop_tsc
### aperfmperf pni pclmulqdq dtes64 monitor ds_cpl vmx smx est tm2 ssse3 cx16 xtpr
### pdcm pcid dca sse4_1 sse4_2 x2apic popcnt tsc_deadline_timer aes xsave avx
### lahf_lm ida arat xsaveopt pln pts dtherm tpr_shadow vnmi flexpriority ept vpid
### bogomips	: 5799.87
### clflush size	: 64
### cache_alignment	: 64
### address sizes	: 46 bits physical, 48 bits virtual
### power management:


### --------------------------------------------------------------------------
### Set up.
### --------------------------------------------------------------------------
library(microbenchmark)
library(VariantAnnotation)
library(Rplinkseq)

## 494328 variants, 1092 samples, 22 INFO, 3 GENO
## ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20101123/
path <- "/loc/no-backup/vobencha/"
fl <- paste0(path,
"ALL.chr22.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz")


### --------------------------------------------------------------------------
### Iterate through complete file.
### --------------------------------------------------------------------------

### Rplinkseq:
### Rplinkseq has no 'yield' or iterating capability. Attempting to 
### read the full file takes > 250 GIG. The file can be read in chunks 
### by genomic position. Looping through positions would require previous
### knowledge of the number of records in each range.

xx <- load.vcf(fl, limit=500000) ## requires > 250 GIG


### scanVcf():
tf <- TabixFile(fl, yieldSize=10000)
open(tf)
system.time( while(length(scanVcf(tf)[[1]]$rowData) > 0)
                 cat("yield\n") )
close(tf)
###     user   system  elapsed 
### 2235.707   27.394 2263.733 

### max memory used (0.053)*(396101596k) = 20 GIG


### readGT():
tf <- TabixFile(fl, yieldSize=10000)
open(tf)
system.time( while(length(readGT(tf)) > 0)
                 cat("yield\n") )
close(tf)
###    user  system elapsed 
### 512.220   1.172 513.525 


### --------------------------------------------------------------------------
### Targeted queries by field and range.
### --------------------------------------------------------------------------

### No comparisions for Rplinkseq because load.vcf() cannot read a 
### specific field from a VCF file.

### Define a range of 133499 records.
tf <- TabixFile(fl)
which <- GRanges("22", IRanges(2e7, 3e7))
param <- ScanVcfParam(which=which)

### (i) 'DS' geno field
fun0 <- function() readGeno(tf, "DS", param=param)
> microbenchmark(fun0(), times=5)
### Unit: seconds
###    expr      min       lq   median       uq      max neval
###  fun0() 44.84508 45.10242 45.15887 45.18403 45.29007     5

### (ii) 'GT' geno field
fun0 <- function() readGT(tf, param=param)
> microbenchmark(fun0(), times=5)
### Unit: seconds
###    expr      min      lq   median       uq      max neval
###  fun0() 107.1373 107.604 108.2043 108.3288 110.3359     5

### (iii) 'THETA' info field
fun0 <- function() readInfo(tf, "THETA", param=param)
> microbenchmark(fun0(), times=5)
### Unit: seconds
###    expr      min       lq   median       uq      max neval
###  fun0() 18.31676 18.32141 18.32868 18.33132 18.34386     5


### --------------------------------------------------------------------------
### Linear scaling with number of variants.
### --------------------------------------------------------------------------

### Testing is done with scanVcf() because it is called by all other 
### read* functions.

tf <- TabixFile(fl)
yieldSize <- c(100000, 200000, 300000, 400000, 500000) 
res <- lapply(yieldSize, function(i) {
              tf <- TabixFile(fl, yieldSize=i)
              param <- ScanVcfParam(info=NA, geno=NA)
              system.time(scanVcf(tf, param=param))})
> res
## [[1]]
##    user  system elapsed 
##   9.485   0.144   9.634 
## 
## [[2]]
##    user  system elapsed 
##  19.453   0.176  19.636 
## 
## [[3]]
##    user  system elapsed 
##  29.150   0.296  29.454 
## 
## [[4]]
##    user  system elapsed 
##  39.502   0.424  39.936 
## 
## [[5]]
##    user  system elapsed 
##  48.703   0.500  49.214 


### --------------------------------------------------------------------------
### Targeted queries by range only.
### --------------------------------------------------------------------------

### Rplinkseq cannot read specific fields from the VCF file and
### reading the entire file with load.vcf() requires > 250 GIG. 
### To get a direct comparision with scanVcf() we query by genomic 
### range only.

### scanVcf() takes approximately 14 minutes to read in the
### subset while the equivalent set of functions in Rplinkseq
### takes approximately 28 minutes.

### Rplinkseq:
hdr <- scanVcfHeader(fl)
info <- rownames(info(hdr))
geno <- rownames(geno(hdr))
mask <- "reg=22:20000000..30000000" ## 133449 records

### load.vcf() returns a list of lists; more processing is
### needed to get at the rough equivalent of scanVcf() output.
fun0 <- function()
{
    xx <- load.vcf(fl, mask=mask, limit=200000)
    ## Collate INFO and GENO
    lapply(info, function(i) x.consensus.meta(xx$VAR, i))
    lapply(geno, function(i) x.consensus.genotype(xx$VAR, i))
}

microbenchmark(fun0(), times=3)
### Unit: seconds
###    expr      min      lq   median       uq      max neval
###  fun0() 1659.694 1666.31 1672.927 1808.175 1943.423     3

### single run memory usage (creating list of lists only)
system.time(xx <- load.vcf(fl, mask=mask, limit=200000))
###    user  system elapsed 
### 798.654  46.878 848.497 
gc()
###             used   (Mb) gc trigger    (Mb)  max used   (Mb)
### Ncells 153940428 8221.4  222628720 11889.7 153941308 8221.4
### Vcells 962320366 7342.0 1050991935  8018.5 962563098 7343.8
print(object.size(xx), units="Gb")
### 13.1 Gb

### max memory during runtime (0.2)*(396101596k) = 75 GIG 


### scanVcf():
tf <- TabixFile(fl)
which <- GRanges("22", IRanges(2e7, 3e7)) ## 133449 records
param <- ScanVcfParam(which=which)

fun0 <- function()
    scn <- scanVcf(tf, param=param)

microbenchmark(fun0(), times=3)
### Unit: seconds
###    expr      min       lq median       uq      max neval
###  fun0() 815.7059 828.9979 842.29 1032.187 1222.084     3

### single run memory usage
system.time(scn <- scanVcf(tf, param=param))
###     user   system  elapsed 
### 1193.947   26.225 1220.582 
gc()
###              used   (Mb) gc trigger    (Mb)   max used    (Mb)
### Ncells  149149655 7965.5  463692252 24763.9  441938172 23602.1
### Vcells 1026383299 7830.7 3609168658 27535.8 3074722381 23458.3
print(object.size(scn), units="Gb")
### 13.1 Gb

### max memory during runtime (0.2)*(396101596k) = 75 GIG 

