### ======================================================================
### Testing Rplinkseq and VariantAnnotation
### ======================================================================
###
### February 2014
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
### Objectives
### --------------------------------------------------------------------------

### Compare runtimes between Rplinkseq and VariantAnnotation packages.
###
###
### Test functions:
###
### load.vcf:    Data are loaded from a vcf file into a list of lists 
###              and accessed with x.consensus* functions. (Rplinkseq) 
###
### var.fetch:   Data are loaded from a PLINK/Seq 'project' into a list of
###              lists and accessed with x.consensus* functions. (Rplinkseq) 
###
### meta.fetch:  Data are loaded from a PLINK/Seq 'project' and parsed into
###              a data.frame. (Rplinkseq)
###
### var.iterate: Applies a function to data from a PLINK/Seq 'project'
###              while iterating. Result of function is returned.
###              (Rplinkseq)
###
### scanVcf:     Data are loaded from a vcf file. Info and geno fields are 
###              parsed into a list of lists; other 'core' fields are
###              returned as a GRAnges object. (VariantAnnotation)
###
###
### Test cases:
###
### Test I:      Query by range.
###              Import data by randomly chosen genomic range. 
###              (Functions: load.vcf, var.fetch, scanVcf)
###
### Test II:     Query by range and fields.
###              Import data by genomic range and 4 randomly
###              chosen info (2) and geno (2) fields.
###              (Functions: meta.fetch, scanVcf)
### 
### Test III.    Iterate through file.
###              A simple function is applied to all records in the
###              file in an iterative fashion.
###              (Functions: var.iterate, scanVcf)
###

### --------------------------------------------------------------------------
### Data
### --------------------------------------------------------------------------

### Test file:
### Size on disk: 1.8G compressed
### Contents: 494328 variants, 1092 samples, 22 INFO, and 3 GENO fields
### Download location: 
### ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20110521/ALL.chr22.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz
###
### Test subset:
### An arbitrary subset of data was chosen for testing because the full file 
### requires >250G in memory. The subset range is 2e7 to 2.5e7 and contains 
### 63088 records. This subset is the 'range' in the mask and param objects. 
###
### PLINK/Seq 'project':
### From the raw vcf on disk we creted a PLINK/Seq 'project' so we could use 
### the var.iterate function. A 'project' creates compressed, indexed, SQLite
### database from input files. This one took ~30 minutes to create.
### pseq proj new-project
### pseq proj load-vcf --vcf 'ALL.chr22.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz'

### --------------------------------------------------------------------------
### Set up 
### --------------------------------------------------------------------------
library(microbenchmark)
library(VariantAnnotation)
library(Rplinkseq)  ## version 0.08
fl <- "/loc/no-backup/vobencha/ALL.chr22.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz"
### Attach the 'project' to the R session:
pseq.project("proj") 

### --------------------------------------------------------------------------
### Test I. Query by range
### --------------------------------------------------------------------------

### load.vcf:
mask <- "reg=22:20000000..25000000"
loadvcf_1 <- function()
    load.vcf(fl, mask=mask, limit=200000)

### var.fetch:
varfetch_1 <- function()
     var.fetch(mask=mask, limit=200000)

### scanVcf:
tfile <- TabixFile(fl)
which <- GRanges("22", IRanges(2e7, 2.5e7))
param <- ScanVcfParam(which=which)
scanvcf_1 <- function() 
    scanVcf(tfile, param=param)

micro <- microbenchmark(loadvcf_1(), varfetch_1(), scanvcf_1(), times=5)
micro
###Unit: seconds
###         expr      min       lq   median       uq      max neval
###  loadvcf_1() 313.1729 335.5738 359.7979 368.4054 369.9397     5
### varfetch_1() 264.4101 283.2388 291.7533 305.7028 318.7879     5
###  scanvcf_1() 300.5585 308.5507 359.0814 400.2049 678.2296     5

### --------------------------------------------------------------------------
### Test II. Query by range and fields
### --------------------------------------------------------------------------

### Randomly select 2 INFO and 2 GENO fields:
info_var <- c("RSQ", "THETA")
geno_var <- c("GT", "DS")
### Other fields parsed by scanVcf:
other_var <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER")

### load.vcf: To the best of our knowledge, this function can import
###           by range but not by select fields. To isolate a field,
###           all fields are read in and parsed with the x.consensus 
###           functions. There is no equivalent comparision for 
###           load.vcf in this test case.

### meta.fetch:
mask <- "reg=22:20000000..25000000"
metafetch_2 <- function()
    meta.fetch(c(info_var, geno_var, other_var), mask=mask) 

### scanVcf:  NOTE: Fields in 'other_var' are imported and parsed by
###           scanVcf() under names of 'rowRanges' and 'fixed'.
tfile <- TabixFile(fl)
which <- GRanges("22", IRanges(2e7, 2.5e7))
param <- ScanVcfParam(which=which, info=info_var, geno=geno_var)
scanvcf_2 <- function()
    scanVcf(tfile, param=param)

micro <- microbenchmark(metafetch_2(), scanvcf_2(), times=5)
micro
###Unit: seconds
###          expr       min        lq    median        uq       max neval
### metafetch_2() 118.14090 120.65009 120.87346 120.88615 121.13340     5
###   scanvcf_2()  35.45512  35.46179  35.51305  35.59804  37.41809     5


### --------------------------------------------------------------------------
### Test III. Iterate through file 
### --------------------------------------------------------------------------

### load.vcf: To the best of our knowledge Rplinkseq cannot iterate on 
###           raw vcf files. Iteration must be done on a 'project' with
###           var.iterate. There is no equivalent comparision for load.vcf
###           in this test case.

### var.iterate:
simple_ct <- function(v)
    ct <<- ct + length(v$ID) 

variterate_3 <- function()
{
    ct <<- 0
    var.iterate(simple_ct)
}

### scanVcf:
tfile <- TabixFile(fl, yieldSize=10000)
param <- ScanVcfParam(fixed=c("ALT"), info=NA, geno=NA)
scanvcf_3 <- function()
{
    ct2 <<- 0
    open(tfile)
    while((len <- length(scanVcf(tfile, param=param)[[1]]$rowRanges)) > 0) 
        ct2 <<- ct2 + len 
    close(tfile)
}

micro <- microbenchmark(variterate_3(), scanvcf_3(), times=5)
micro
##Unit: seconds
##           expr        min         lq     median         uq        max neval
## variterate_3() 1552.88520 1576.55856 1583.13120 1585.65281 1591.96375     5
##    scanvcf_3()   50.04566   50.06945   50.27194   50.56145   50.70112     5


### --------------------------------------------------------------------------
### Scaling with number of variants and samples (scanVcf only)
### --------------------------------------------------------------------------

### Linear scaling by variants:
tfile <- TabixFile(fl)
yieldSize <- c(100000, 200000, 300000, 400000, 500000) 
param <- ScanVcfParam(info=NA, geno=NA)
fun0 <- function(tfile, param) scanVcf(tfile, param=param)

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
