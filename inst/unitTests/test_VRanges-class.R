.TARGET_seqnames <- Rle(factor(c("chr1", "chr2")))
.TARGET_ranges <- IRanges(c(1, 10), c(5, 20))
.TARGET_strand <- Rle(strand("+"), 2)
.TARGET_gr <- GRanges(.TARGET_seqnames, .TARGET_ranges, .TARGET_strand)
.TARGET_ref <- c("T", "A")
.TARGET_alt <- c("C", "T")
.TARGET_refDepth <- as.integer(c(5, 10))
.TARGET_altDepth <- as.integer(c(7, 6))
.TARGET_totalDepth <- as.integer(c(12, 17))
.TARGET_sampleNames <- factor(letters[1:2])
.TARGET_hardFilters <- FilterRules(list(a = function(x) softFilterMatrix(x)[,1]))
.TARGET_softFilterMatrix <- FilterMatrix(matrix = cbind(a = c(TRUE, FALSE)),
                                         filterRules = .TARGET_hardFilters)
.TARGET_tumorSpecific <- c(FALSE, TRUE)

integerRle <- function(...) as(Rle(...), "integerRle")
characterRle <- function(...) as(Rle(...), "characterRle")
factorRle <- function(...) as(Rle(...), "factorRle")

make_TARGET_VRanges_simple <- function() {
  new("VRanges", .TARGET_gr, ref = .TARGET_ref,
      alt = .TARGET_alt, totalDepth = integerRle(NA_integer_, 2),
      refDepth = integerRle(NA_integer_, 2),
      altDepth = integerRle(NA_integer_, 2),
      sampleNames = factorRle(NA_character_, 2),
      softFilterMatrix =
        FilterMatrix(matrix(nrow = length(.TARGET_gr), ncol = 0L),
                     filterRules = FilterRules()),
      hardFilters = FilterRules())
}

make_TARGET_VRanges <- function(i = seq_len(length(.TARGET_seqnames))) {
  new("VRanges",
      GRanges(.TARGET_seqnames, .TARGET_ranges, .TARGET_strand,
              tumorSpecific = .TARGET_tumorSpecific)[i],
      ref = .TARGET_ref[i],
      alt = .TARGET_alt[i], totalDepth = .TARGET_totalDepth[i],
      refDepth = .TARGET_refDepth[i], altDepth = .TARGET_altDepth[i],
      sampleNames = .TARGET_sampleNames[i],
      softFilterMatrix = .TARGET_softFilterMatrix[i,,drop=FALSE],
      hardFilters = .TARGET_hardFilters)
}

make_TARGET_VRanges_empty <- function() {
  new("VRanges", GRanges(), ref = character(),
      alt = characterRle(character()), totalDepth = integerRle(integer()),
      refDepth = integerRle(integer()),
      altDepth = integerRle(integer()),
      sampleNames = factorRle(factor()),
      softFilterMatrix = FilterMatrix(matrix(nrow = 0L, ncol = 0L),
        filterRules = FilterRules()),
      hardFilters = FilterRules())
}

test_VRanges_constructor <- function() {

  vr <- make_TARGET_VRanges_simple()

  ## CHECK: 'ref' NA
  ranges <- IRanges(c(1, 5), c(10, 20))
  checkException(VRanges(.TARGET_seqnames, .TARGET_ranges, NA_character_))
  
  ## CHECK: depth < 0
  checkException(VRanges(.TARGET_seqnames, .TARGET_ranges, .TARGET_ref,
                         refDepth = -1))
  checkException(VRanges(.TARGET_seqnames, .TARGET_ranges, .TARGET_ref,
                         altDepth = -1))
  
  ## CHECK: invalid totalDepth sum
  options(warn=2)
  checkException(VRanges(.TARGET_seqnames, .TARGET_ranges,
                         .TARGET_ref, .TARGET_alt, 0, 1))
  options(warn=0)
  
  ## CHECK: DNAStringSet handling (ref and alt); also, coercion depth=>integer
  test.vr <- VRanges(.TARGET_seqnames, .TARGET_ranges,
                     DNAStringSet(.TARGET_ref), DNAStringSet(.TARGET_alt),
                     totalDepth = .TARGET_totalDepth)
  vr@totalDepth <- .TARGET_totalDepth
  checkIdentical(test.vr, vr)
  
  ## CHECK: recycling (every column)
  test.vr <- VRanges(.TARGET_seqnames[1], .TARGET_ranges,
                     DNAStringSet(.TARGET_ref[1]), .TARGET_alt[1],
                     tumorSpecific = .TARGET_tumorSpecific,
                     totalDepth = .TARGET_totalDepth)
  seqnames(vr) <- .TARGET_seqnames[1]
  vr$tumorSpecific <- .TARGET_tumorSpecific
  vr@ref <- rep(.TARGET_ref[1], 2)
  alt(vr) <- .TARGET_alt[1]
  checkIdentical(test.vr, vr)

  ## CHECK: coercion of 0/1 matrix to logical
  test.vr <- VRanges(.TARGET_seqnames, .TARGET_ranges, .TARGET_ref, .TARGET_alt,
                     softFilterMatrix = cbind(a = c(1, 0)))
  vr <- make_TARGET_VRanges_simple()
  vr@softFilterMatrix <- as(.TARGET_softFilterMatrix, "matrix")
  checkIdentical(test.vr, vr)
  
  ## CHECK: all arguments, with unnamed metadata
  test.vr <- VRanges(.TARGET_seqnames, .TARGET_ranges, .TARGET_ref,
                     .TARGET_alt, .TARGET_totalDepth,
                     .TARGET_refDepth, .TARGET_altDepth,
                     tumorSpecific = .TARGET_tumorSpecific,
                     sampleNames = .TARGET_sampleNames,
                     softFilterMatrix = .TARGET_softFilterMatrix,
                     hardFilters = .TARGET_hardFilters)
  vr <- make_TARGET_VRanges()
  checkIdentical(test.vr, vr)

  ## CHECK: zero-length input
  test.vr <- VRanges()
  empty.vr <- make_TARGET_VRanges_empty()
  checkIdentical(test.vr, empty.vr)
}

test_VRanges_extract <- function() {
  vr <- make_TARGET_VRanges()
  checkIdentical(vr[1], make_TARGET_VRanges(1))
  checkIdentical(vr[0], make_TARGET_VRanges(0))
  checkIdentical(vr[rep(1, 3)], make_TARGET_VRanges(rep(1, 3)))
  checkIdentical(vr[c(2, 1)], make_TARGET_VRanges(c(2, 1)))
  checkIdentical(vr[c(TRUE,FALSE)], make_TARGET_VRanges(c(TRUE,FALSE)))
  checkIdentical(vr[TRUE], make_TARGET_VRanges(TRUE))
  checkIdentical(vr[], make_TARGET_VRanges())
  checkIdentical(vr[logical()], make_TARGET_VRanges(logical()))
  vr$tumorSpecific <- NULL
  checkIdentical(vr, make_TARGET_VRanges()[,integer()])
}

test_VRanges_accessors <- function() {
  vr <- make_TARGET_VRanges()
  
  checkIdentical(seqnames(vr), .TARGET_seqnames)
  checkIdentical(ranges(vr), .TARGET_ranges)
  checkIdentical(strand(vr), .TARGET_strand)
  checkIdentical(ref(vr), .TARGET_ref)
  checkIdentical(alt(vr), .TARGET_alt)
  checkIdentical(refDepth(vr), .TARGET_refDepth)
  checkIdentical(altDepth(vr), .TARGET_altDepth)
  checkIdentical(totalDepth(vr), .TARGET_totalDepth)
  checkIdentical(Biobase::sampleNames(vr), .TARGET_sampleNames)
  checkIdentical(softFilterMatrix(vr), .TARGET_softFilterMatrix)
  checkIdentical(hardFilters(vr), .TARGET_hardFilters)
  checkIdentical(vr$tumorSpecific, .TARGET_tumorSpecific)

  alt(vr) <- NA
  checkIdentical(alt(vr), characterRle(NA, 2))
  twoFilters <- rep(.TARGET_hardFilters, 2)
  hardFilters(vr) <- twoFilters
  checkIdentical(hardFilters(vr), twoFilters)
  twoFilters <- .TARGET_softFilterMatrix[,c(1,1)]
  softFilterMatrix(vr) <- twoFilters
  checkIdentical(softFilterMatrix(vr), twoFilters)
}

test_VRanges_subassign <- function() {
  vr <- make_TARGET_VRanges()
  vr[2:1] <- vr
  checkIdentical(vr, make_TARGET_VRanges(2:1))
  
  vr <- make_TARGET_VRanges()
  sampleNames(vr) <- Rle(sampleNames(vr))
  vr2 <- vr
  vr2[2:1] <- vr
  checkIdentical(vr2, vr[2:1])
}

test_VRanges_combine <- function() {
  vr <- make_TARGET_VRanges()
  vr.rep <- vr[rep(seq_len(length(vr)), 2)]
  hardFilters(vr.rep) <- FilterRules()
  vr.combined <- c(vr, vr)
  checkIdentical(vr.combined, vr.rep)
}

test_VRanges_coerce <- function() {
   df <- data.frame(as.data.frame(.TARGET_gr), ref = as.vector(.TARGET_ref),
                    alt = as.vector(.TARGET_alt),
                    totalDepth = as.vector(.TARGET_totalDepth),
                    refDepth = as.vector(.TARGET_refDepth),
                    altDepth = as.vector(.TARGET_altDepth),
                    sampleNames = as.factor(.TARGET_sampleNames),
                    softFilterMatrix = .TARGET_softFilterMatrix,
                    tumorSpecific = .TARGET_tumorSpecific,
                    stringsAsFactors = FALSE)
   vr <- make_TARGET_VRanges()
   checkIdentical(as.data.frame(vr), df)
}

make_TARGET_VRanges_vcf <- function() {
  gr <- GRanges(c("1", "1", "3", "2", "1", "2"),
                IRanges(c(1, 20, 15, 5, 20, 5), width = 1L), "+",
                TS = c(TRUE, FALSE, FALSE, TRUE, TRUE, NA))
  new("VRanges", gr, ref = c("A", "C", "G", "T", "C", "T"),
      alt = c("T", "G", NA, "A", "G", "C"),
      totalDepth = c(NA, 10L, 15L, 10L, 5L, 3L),
      refDepth = c(0L, 5L, 15L, 10L, NA, 2L),
      altDepth = c(5L, NA, NA, 0L, 3L, 1L),
      sampleNames = factorRle(factor(c("A", "B", "B", "A", "A", "B"))),
      softFilterMatrix =
      FilterMatrix(matrix = cbind(a = c(FALSE, TRUE, FALSE, TRUE, FALSE, TRUE),
                     b = c(FALSE, FALSE, TRUE, TRUE, NA, NA)),
                   filterRules = FilterRules(a = x > 1, b = c != "foo")),
      hardFilters = FilterRules())
}

checkIdenticalVCF <- function(orig, vcf) {
  vcf$QUAL <- NULL
  if (!is.null(orig$TS))
    orig$TS <- as.integer(orig$TS)
  softFilterMatrix(orig) <- as(softFilterMatrix(orig), "matrix")
  ## Information loss due to binary nature of VCF filters
  if (!identical(softFilterMatrix(orig), softFilterMatrix(vcf))) {
    origFilt <- softFilterMatrix(orig)
    origFilt[is.na(origFilt)] <- TRUE
    origFilt[rowSums(origFilt) == 2 &
             rowSums(softFilterMatrix(orig),na.rm=TRUE) != 2,] <- NA
    softFilterMatrix(orig) <- origFilt
  }
  checkIdentical(orig, vcf)
}

test_VRanges_vcf <- function() {
  dest <- tempfile()
  vr <- make_TARGET_VRanges_vcf()
  writeVcf(vr, dest)
  vcf <- readVcf(dest, genome = "hg19")
  perm <- c(1, 7, 8, 4, 2, 10)
  vcf.vr <- as(vcf, "VRanges")[perm]
  checkIdenticalVCF(vr, vcf.vr)
  
  test.vcf <- asVCF(vr, info = "TS", filter = "a")
  writeVcf(test.vcf, dest)
  vcf <- readVcf(dest, genome = "hg19")
  vcf.vr <- as(vcf, "VRanges")[perm]
  ## adjustments necessary, due to folding at INFO and FILTER
  vcf.vr$TS[5:6] <- c(1L, NA)
  softFilterMatrix(vcf.vr)[5,] <- cbind(FALSE, NA)
  checkIdenticalVCF(vr, vcf.vr)
  
  vrA <- vr[sampleNames(vr) == "A"]
  runValue(sampleNames(vrA)) <- factor(runValue(sampleNames(vrA)))
  vrA <- keepSeqlevels(vrA, unique(as.character(seqnames(vrA))))
  writeVcf(vrA, dest)
  vcfA <- readVcf(dest, genome = "hg19")
  vcfA.vr <- as(vcfA, "VRanges")
  checkIdenticalVCF(vrA, vcfA.vr)

  geno(vcfA) <- SimpleList()
  vrA.vcf <- as(vcfA, "VRanges")
  vrA.stripped <- VRanges(seqnames(vrA), ranges(vrA), ref(vrA), alt(vrA),
                          sampleNames = sampleNames(vrA),
                          softFilterMatrix = matrix(nrow = length(vrA),
                            ncol = 0L))
  checkIdenticalVCF(vrA.stripped, vrA.vcf)
}


### GOT TIRED OF TESTING. I mean, seriously, that VCF stuff was insane.
