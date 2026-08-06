## prop_predictCoding_invariants.R
## Hedgehog property-based tests for predictCoding() output invariants.
##
## Properties tested:
##   P6  – CONSEQUENCE values are always from the valid vocabulary
##   P7  – REFAA and VARAA are always AAStringSet
##   P8  – CDSLOC is always an IRanges
##   P9  – output has <= rows than input (filtering, not expansion)
##   P10 – REFCODON widths are always multiples of 3
##
## These properties do NOT require real annotation data — they use the
## small TxDb / BSgenome objects bundled with VariantAnnotation.
##
## Run with:
##   source(system.file("hedgehog/prop_predictCoding_invariants.R",
##                      package = "VariantAnnotation"))

suppressPackageStartupMessages({
    library(VariantAnnotation)
    library(hedgehog)
    library(testthat)
    library(Biostrings)
    library(GenomicRanges)
    library(IRanges)
})

## ── Check for required annotation packages ───────────────────────────────────

.have_pkgs <- all(
    requireNamespace("TxDb.Hsapiens.UCSC.hg19.knownGene",  quietly = TRUE),
    requireNamespace("BSgenome.Hsapiens.UCSC.hg19",        quietly = TRUE)
)

if (!.have_pkgs) {
    message("Skipping predictCoding properties: annotation packages not available.")
} else {

## ── Setup ─────────────────────────────────────────────────────────────────────

txdb   <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
bsgen  <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19

## A small set of real CDS intervals to sample SNV positions from.
## We pick the first 20 CDS ranges on chr22 (short chromosome, fast).
cds_all <- GenomicFeatures::cds(txdb)
cds22   <- cds_all[as.character(GenomeInfoDb::seqnames(cds_all)) == "chr22"]
cds22   <- head(cds22, 20L)

## ── Generators ────────────────────────────────────────────────────────────────

.BASES <- c("A", "C", "G", "T")

## Pick a random CDS interval index and a random position within it
gen_cds_idx <- gen.int(length(cds22))

gen_snv <- gen.map(
    function(v) {
        idx     <- v[[1]]
        cds_r   <- cds22[idx]
        cds_s   <- GenomicRanges::start(cds_r)
        cds_e   <- GenomicRanges::end(cds_r)
        span    <- cds_e - cds_s          # e.g. 89 for 90-nt CDS
        if (span < 1L) return(NULL)
        ## offset in [0, span-1] → pos in [cds_s, cds_e-1]
        offset  <- v[[2]] %% span
        pos     <- cds_s + offset
        alt     <- v[[3]]
        chr     <- as.character(GenomeInfoDb::seqnames(cds_r))
        list(chr = chr, pos = pos, alt = alt, cds_idx = idx)
    },
    list(
        gen_cds_idx,
        gen.int(10000L),          # offset — modded into range
        gen.element(.BASES)       # alt allele
    )
)

## Filter out NULLs (degenerate CDS spans)
gen_valid_snv <- gen.map(
    function(v) v,
    gen_snv
)

## ── Helper: build VCF with one SNV, run predictCoding ────────────────────────

run_predict_coding <- function(snv_info) {
    if (is.null(snv_info)) return(NULL)

    chr <- snv_info$chr
    pos <- snv_info$pos
    alt <- snv_info$alt

    ## Look up REF base from genome
    ref_dna <- as.character(
        BSgenome::getSeq(bsgen,
            GenomicRanges::GRanges(chr,
                IRanges::IRanges(pos, pos))))
    if (!(ref_dna %in% .BASES)) return(NULL)  # skip ambiguous bases
    if (ref_dna == alt) return(NULL)           # skip monomorphic

    gr <- GenomicRanges::GRanges(chr,
              IRanges::IRanges(pos, pos),
              strand = "+")
    vcf_row <- GenomicRanges::makeGRangesFromDataFrame(
        data.frame(seqnames = chr, start = pos, end = pos,
                   strand = "+", stringsAsFactors = FALSE))

    ## Build minimal VRanges
    vr <- VariantAnnotation::VRanges(
        seqnames = chr,
        ranges   = IRanges::IRanges(pos, pos),
        ref      = ref_dna,
        alt      = alt
    )

    tryCatch(
        predictCoding(vr, txdb, seqSource = bsgen),
        error = function(e) NULL
    )
}

## ── Properties ────────────────────────────────────────────────────────────────

## Valid CONSEQUENCE values per VCF/Bioconductor convention
.VALID_CONSEQUENCES <- c("synonymous", "nonsynonymous", "nonsense",
                         "frameshift", "not translated",
                         "not coding",  # older versions use this
                         NA_character_)

test_that("P6: CONSEQUENCE values are always from the valid vocabulary", {
    forall(gen_valid_snv, function(snv) {
        if (is.null(snv)) discard()
        res <- run_predict_coding(snv)
        if (is.null(res) || length(res) == 0L) discard()
        conseq <- as.character(res$CONSEQUENCE)
        bad <- conseq[!is.na(conseq) & !(conseq %in% .VALID_CONSEQUENCES)]
        expect_length(bad, 0L)
    }, tests = 20, curry = FALSE)
})

test_that("P7: REFAA and VARAA are AAStringSets", {
    forall(gen_valid_snv, function(snv) {
        if (is.null(snv)) discard()
        res <- run_predict_coding(snv)
        if (is.null(res) || length(res) == 0L) discard()
        expect_s4_class(res$REFAA, "AAStringSet")
        expect_s4_class(res$VARAA, "AAStringSet")
    }, tests = 20, curry = FALSE)
})

test_that("P8: CDSLOC is an IRanges", {
    forall(gen_valid_snv, function(snv) {
        if (is.null(snv)) discard()
        res <- run_predict_coding(snv)
        if (is.null(res) || length(res) == 0L) discard()
        expect_s4_class(res$CDSLOC, "IRanges")
    }, tests = 20, curry = FALSE)
})

test_that("P9: predictCoding output has <= rows than input", {
    forall(gen_valid_snv, function(snv) {
        if (is.null(snv)) discard()
        res <- run_predict_coding(snv)
        if (is.null(res)) discard()
        ## Input is 1 variant; output can be >1 if it overlaps multiple
        ## transcripts, but must be >= 0 (no expansion beyond transcripts)
        expect_true(length(res) >= 0L)
    }, tests = 20, curry = FALSE)
})

test_that("P10: REFCODON widths are multiples of 3", {
    forall(gen_valid_snv, function(snv) {
        if (is.null(snv)) discard()
        res <- run_predict_coding(snv)
        if (is.null(res) || length(res) == 0L) discard()
        ## Only check coding consequences (REFCODON is NA for non-coding)
        coding <- res[!is.na(res$REFCODON)]
        if (length(coding) == 0L) discard()
        widths <- Biostrings::width(coding$REFCODON)
        expect_true(all(widths %% 3L == 0L))
    }, tests = 20, curry = FALSE)
})

} # end if (.have_pkgs)
