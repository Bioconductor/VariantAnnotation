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
library(Rplinkseq)  ## version 0.08

## 494328 variants, 1092 samples, 22 INFO, 3 GENO
## ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20110521/ALL.chr22.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz
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
> gc(reset=TRUE)
          used  (Mb) gc trigger  (Mb) max used  (Mb)
Ncells 2600229 138.9    3289005 175.7  2600229 138.9
Vcells 1951494  14.9    3124965  23.9  1951494  14.9


tf <- TabixFile(fl, yieldSize=10000)
open(tf)
system.time( while(length(scanVcf(tf)[[1]]$rowData) > 0) {} )
close(tf)
###     user   system  elapsed 
### 2235.707   27.394 2263.733 

> gc(reset=TRUE)
          used  (Mb) gc trigger   (Mb) max used  (Mb)
Ncells 2620925 140.0   24923264 1331.1  2620925 140.0
Vcells 1971533  15.1  182613726 1393.3  1971533  15.1

## memory used : (138.9+14.9) - (140+15.1) = 1.3 Mb 


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

### Linear scaling by variants:
tf <- TabixFile(fl)
yieldSize <- c(100000, 200000, 300000, 400000, 500000) 
param <- ScanVcfParam(info=NA, geno=NA)
fun0 <- function(tf, param) scanVcf(tf, param=param)

res <- lapply(yieldSize, function(i) {
              tf <- TabixFile(fl, yieldSize=i)
              microbenchmark(fun0(tf, param), times=5)})
res
### [[1]]
### Unit: seconds
###             expr      min       lq   median       uq      max neval
###  fun0(tf, param) 9.883719 10.13574 10.15165 10.17318 10.70124     5
### 
### [[2]]
### Unit: seconds
###             expr    min      lq   median       uq      max neval
###  fun0(tf, param) 19.385 19.6715 19.71968 19.72206 20.23362     5
### 
### [[3]]
### Unit: seconds
###             expr      min       lq   median       uq      max neval
###  fun0(tf, param) 29.24268 29.33945 29.36927 29.38207 30.01012     5
### 
### [[4]]
### Unit: seconds
###             expr      min       lq   median      uq      max neval
###  fun0(tf, param) 38.71675 38.73608 38.77965 38.8696 39.27382     5
### 
### [[5]]
### Unit: seconds
###             expr      min       lq   median       uq      max neval
###  fun0(tf, param) 47.94605 48.76884 49.03272 49.03859 49.08618     5


### Linear scaling by sample:
tf <- TabixFile(fl)
ids <- samples(scanVcfHeader(fl))
sampleSize <- c(200, 400, 600, 800, 1000)
which <- GRanges("22", IRanges(2e7, 2.5e7)) ## 63088 records
fun0 <- function(tf, param) scanVcf(tf, param=param) 
 
res <- lapply(sampleSize, function(i) {
              param <- ScanVcfParam(which=which, samples=ids[1:i])
              microbenchmark(fun0(tf, param), times=5)})
res
### [[1]]
### Unit: seconds
###             expr      min       lq   median       uq      max neval
###  fun0(tf, param) 61.80377 61.91141 69.62993 70.10218 103.1949     5
### 
### [[2]]
### Unit: seconds
###             expr      min       lq   median      uq      max neval
###  fun0(tf, param) 128.3595 132.1819 135.7386 135.966 150.7588     5
### 
### [[3]]
### Unit: seconds
###             expr      min       lq   median       uq     max neval
###  fun0(tf, param) 204.0067 211.0546 219.2069 232.3388 249.188     5
### 
### [[4]]
### Unit: seconds
###             expr      min       lq   median       uq      max neval
###  fun0(tf, param) 266.7383 269.8461 273.7986 299.7805 330.6505     5
### 
### [[5]]
### Unit: seconds
###             expr      min       lq   median       uq      max neval
###  fun0(tf, param) 329.6348 337.1841 348.7968 359.5972 364.2781     5


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
gc(reset=TRUE)
###           used  (Mb) gc trigger  (Mb) max used  (Mb)
### Ncells 2600229 138.9    3289005 175.7  2600229 138.9
### Vcells 1951494  14.9    3124965  23.9  1951494  14.9

xx <- load.vcf(fl, mask=mask, limit=200000)

gc(reset=TRUE)
###             used   (Mb) gc trigger   (Mb)  max used   (Mb)
### Ncells  74151020 3960.1  106673040 5697.0  74151020 3960.1
### Vcells 455979993 3478.9  502969959 3837.4 455979993 3478.9

### memory used: (138.9 + 14.9) - (3960.1 + 3478.9) = 7285.2 Mb

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
gc(reset=TRUE)
###           used  (Mb) gc trigger  (Mb) max used  (Mb)
### Ncells 2616411 139.8    3289005 175.7  2616411 139.8
### Vcells 1965798  15.0    3124965  23.9  1965798  15.0

scn <- scanVcf(tf, param=param)

gc(reset=TRUE)
###             used   (Mb) gc trigger    (Mb)  max used   (Mb)
### Ncells  71895814 3839.7  211989258 11321.5  71895814 3839.7
### Vcells 486274427 3710.0 1601055918 12215.1 486274427 3710.0

### memory used: (139.8 + 15.0) - (3839.7 + 3710) = 7394.9 Mb

print(object.size(scn), units="Mb")
### 6343.3 Mb
