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

### (i) Rplinkseq:
### Rplinkseq has no 'yield' or iterating capability. Attempting to 
### read the full file takes > 250 GIG. The file can be read in chunks 
### by genomic position. Looping through positions would require previous
### knowledge of the number of records in each range.

xx <- load.vcf(fl, limit=500000) ## requires > 250 GIG


### (ii) scanVcf():
tf <- TabixFile(fl, yieldSize=10000)
open(tf)
system.time( while(length(scanVcf(tf)[[1]]$rowData) > 0) {} )
close(tf)
###     user   system  elapsed 
### 2235.707   27.394 2263.733 

### max memory used (0.053)*(396101596k) = 20 GIG


### (iii) readGT():
tf <- TabixFile(fl, yieldSize=10000)
open(tf)
system.time( while(length(readGT(tf)) > 0) {} )
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
### Scaling with number of variants and samples.
### --------------------------------------------------------------------------

### Testing done with scanVcf() because it underlies all read* functions.

### Linear scaling with number of variants:
tf <- TabixFile(fl)
yieldSize <- c(100000, 200000, 300000, 400000, 500000) 
res <- lapply(yieldSize, function(i) {
              tf <- TabixFile(fl, yieldSize=i)
              param <- ScanVcfParam(info=NA, geno=NA)
              system.time(scanVcf(tf, param=param))})
res
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

### Approximately linear scaling through ~ n=600 
### and nlog(n) for n > 600
tf <- TabixFile(fl)
ids <- samples(scanVcfHeader(fl))
sampleSize <- c(200, 400, 600, 800, 1000) 
which <- GRanges("22", IRanges(2e7, 2.5e7)) ## 63088 records
res <- lapply(sampleSize, function(i) {
              param <- ScanVcfParam(which=which, samples=ids[1:i])
              system.time(scanVcf(tf, param=param))})

res
### [[1]]
###    user  system elapsed 
###  48.539   0.068  48.617 
### 
### [[2]]
###    user  system elapsed 
### 101.903   0.512 102.439 
### 
### [[3]]
###    user  system elapsed 
### 180.211   0.084 180.365 
### 
### [[4]]
###    user  system elapsed 
### 320.100   1.456 321.686 
### 
### [[5]]
###    user  system elapsed 
### 378.144   4.816 383.275 

### --------------------------------------------------------------------------
### Targeted queries by range only.
### --------------------------------------------------------------------------

### Rplinkseq cannot read specific fields from the VCF file and
### reading the entire file with load.vcf() requires > 250 GIG. 
### To get a direct comparision with scanVcf() we query by genomic 
### range only.

### A subset of 63088 records was chosen as a sufficient number
### to compare runtime and memory use within a reasonable
### time.

### scanVcf() took approximately 5.4 minutes to read all data in
### the subset while the equivalent set of functions in Rplinkseq
### took approximately 12.6 minutes.

### (i) Rplinkseq:
hdr <- scanVcfHeader(fl)
info <- rownames(info(hdr))
geno <- rownames(geno(hdr))
mask <- "reg=22:20000000..25000000" ## 63088 record-dev

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
###    expr     min       lq   median       uq      max neval
###  fun0() 729.048 743.8944 758.7409 823.2939 887.8469     3

### single run memory usage (list of lists only)
xx <- load.vcf(fl, mask=mask, limit=200000)
gc()
###             used   (Mb) gc trigger   (Mb)  max used   (Mb)
### Ncells  74163160 3960.8  106673040 5697.0  74164017 3960.8
### Vcells 455995498 3479.0  503020304 3837.8 455998211 3479.0

print(object.size(xx), units="Mb")
### 6322.6 Mb


### (ii) scanVcf():
tf <- TabixFile(fl)
which <- GRanges("22", IRanges(2e7, 2.5e7)) ## 63088 records
param <- ScanVcfParam(which=which)
fun0 <- function()
    scn <- scanVcf(tf, param=param)
microbenchmark(fun0(), times=3)
### Unit: seconds
###    expr      min       lq  median       uq      max neval
###  fun0() 292.5651 309.4915 326.418 422.9329 519.4477     3

### single run memory usage
scn <- scanVcf(tf, param=param)
gc()
###             used   (Mb) gc trigger    (Mb)   max used    (Mb)
### Ncells  71896387 3839.7  211989258 11321.5  211989258 11321.5
### Vcells 486275159 3710.0 1601058158 12215.2 1667942560 12725.4

print(object.size(scn), units="Mb")
### 6343.3 Mb
