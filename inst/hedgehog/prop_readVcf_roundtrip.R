## prop_readVcf_roundtrip.R
## Hedgehog property-based tests for writeVcf / readVcf round-trip correctness.
##
## Properties tested:
##   P1 – row count is preserved
##   P2 – rownames (CHROM:POS_REF/ALT) match expected pattern
##   P3 – seqnames are always the expected chromosome
##   P4 – ALT alleles survive the round-trip
##   P5 – expand() produces >= as many rows as input (monotonicity)
##
## Run with:
##   source(system.file("hedgehog/prop_readVcf_roundtrip.R",
##                      package = "VariantAnnotation"))

suppressPackageStartupMessages({
    library(VariantAnnotation)
    library(hedgehog)
    library(testthat)
})

## ── Generators ────────────────────────────────────────────────────────────────

.BASES <- c("A", "C", "G", "T")

## Single VCF record generator: named list with ref, alt, pos, qual.
## gen.map(fn, gen): fn is applied to the generated value.
## Named list of generators produces a named list value.
gen_vcf_record <- gen.map(
    function(v) list(ref = v[[1]], alt = v[[2]],
                     pos = v[[3]], qual = v[[4]]),
    list(
        gen.element(.BASES),                              # ref
        gen.element(.BASES),                              # alt
        gen.map(function(i) i + 1L, gen.int(999999L)),   # pos in [1, 1e6]
        gen.element(c(NA_real_, 10, 20, 30, 99))          # qual
    )
)

## List of 1-5 records
gen_records <- gen.list(gen_vcf_record, from = 1L, to = 5L)

## ── Helper ────────────────────────────────────────────────────────────────────

## Build a minimal VCF string and read it back.
## Caller must have already filtered/validated records via discard().
vcf_roundtrip <- function(records) {

    ## Deduplicate by key
    keys <- vapply(records, function(r)
        paste0("chr1:", r$pos, "_", r$ref, "/", r$alt), character(1L))
    records <- records[!duplicated(keys)]
    keys    <- keys[!duplicated(keys)]

    lines <- c(
        "##fileformat=VCFv4.1",
        "##FILTER=<ID=PASS,Description=\"All filters passed\">",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
    )
    for (r in records) {
        qs <- if (is.na(r$qual)) "." else as.character(r$qual)
        lines <- c(lines,
            paste("chr1", r$pos, ".", r$ref, r$alt, qs, "PASS", ".",
                  sep = "\t"))
    }

    tmp <- tempfile(fileext = ".vcf")
    on.exit(unlink(tmp), add = TRUE)
    writeLines(lines, tmp)

    list(vcf = readVcf(tmp, genome = "hg19"),
         n   = length(records),
         keys = keys)
}

## ── Properties ────────────────────────────────────────────────────────────────

## Precondition helper: discard if all records are monomorphic
valid_records <- function(records) {
    recs <- Filter(function(r) r$ref != r$alt, records)
    ## Deduplicate by key
    keys <- vapply(recs, function(r)
        paste0("chr1:", r$pos, "_", r$ref, "/", r$alt), character(1L))
    recs[!duplicated(keys)]
}

test_that("P1: row count preserved through VCF round-trip", {
    forall(gen_records, function(records) {
        recs <- valid_records(records)
        if (length(recs) == 0L) discard()
        rt <- vcf_roundtrip(recs)
        expect_equal(nrow(rt$vcf), rt$n)
    }, tests = 50, curry = FALSE)
})

test_that("P2: rownames match CHROM:POS_REF/ALT pattern", {
    forall(gen_records, function(records) {
        recs <- valid_records(records)
        if (length(recs) == 0L) discard()
        rt <- vcf_roundtrip(recs)
        rn <- rownames(rt$vcf)
        expect_true(all(grepl("^chr1:[0-9]+_[ACGT]/[ACGT]$", rn)))
    }, tests = 50, curry = FALSE)
})

test_that("P3: seqnames are always chr1 for synthetic VCF", {
    forall(gen_records, function(records) {
        recs <- valid_records(records)
        if (length(recs) == 0L) discard()
        rt <- vcf_roundtrip(recs)
        expect_true(all(as.character(seqnames(rt$vcf)) == "chr1"))
    }, tests = 50, curry = FALSE)
})

test_that("P4: ALT alleles are single nucleotides from {A,C,G,T}", {
    forall(gen_records, function(records) {
        recs <- valid_records(records)
        if (length(recs) == 0L) discard()
        rt <- vcf_roundtrip(recs)
        alt_chars <- as.character(unlist(alt(rt$vcf)))
        expect_true(all(alt_chars %in% .BASES))
    }, tests = 50, curry = FALSE)
})

test_that("P5: expand() produces >= nrow(vcf) rows (monotonicity)", {
    forall(gen_records, function(records) {
        recs <- valid_records(records)
        if (length(recs) == 0L) discard()
        rt <- vcf_roundtrip(recs)
        expect_true(nrow(expand(rt$vcf)) >= nrow(rt$vcf))
    }, tests = 50, curry = FALSE)
})
