## Tests for issue #81: locateVariants() and predictCoding() silently drop
## large INDELs that span multiple exons.

## Create a synthetic transcript with 3 exons separated by introns:
##   Exon 1: 1000-1200 (201 bp)
##   Intron: 1201-1499
##   Exon 2: 1500-1700 (201 bp)
##   Intron: 1701-1999
##   Exon 3: 2000-2200 (201 bp)
##
## A 2542bp deletion starting at position 1000 will span all 3 exons.

cdsbytx <- GRangesList(
    tx1 = GRanges(seqnames = "chr1",
                  ranges = IRanges(start = c(1000, 1500, 2000),
                                   end = c(1200, 1700, 2200)),
                  strand = "+",
                  cds_id = c(1L, 2L, 3L),
                  exon_rank = c(1L, 2L, 3L))
)

test_locateVariants_large_indel_not_dropped <- function()
{
    ## A deletion spanning all 3 exons (width > any single exon)
    ## This was previously silently dropped by mapToTranscripts.
    query <- GRanges("chr1", IRanges(start = 1000, width = 2543))

    ## locateVariants with CodingVariants should now include this variant
    ## (previously it returned 0 rows)
    result <- withCallingHandlers(
        VariantAnnotation:::.makeResult(query, cdsbytx, "coding",
                                        ignore.strand = TRUE, asHits = FALSE),
        warning = function(w) {
            ## Verify that a warning about multi-exon spanning is emitted
            checkTrue(grepl("span multiple exons", conditionMessage(w)))
            invokeRestart("muffleWarning")
        }
    )

    ## The variant should appear in results
    checkTrue(length(result) > 0L)
    ## QUERYID should reference our variant (index 1)
    checkTrue(1L %in% mcols(result)$QUERYID)
    ## LOCATION should be "coding"
    checkTrue(all(mcols(result)$LOCATION == "coding"))
    ## LOCSTART/LOCEND should be NA (can't map to transcript coords)
    rescued <- result[mcols(result)$QUERYID == 1L]
    checkTrue(all(is.na(mcols(rescued)$LOCSTART)))
    checkTrue(all(is.na(mcols(rescued)$LOCEND)))
}

test_locateVariants_small_indel_still_works <- function()
{
    ## A small variant within a single exon should still work normally
    query <- GRanges("chr1", IRanges(start = 1050, width = 10))
    result <- VariantAnnotation:::.makeResult(query, cdsbytx, "coding",
                                              ignore.strand = TRUE,
                                              asHits = FALSE)
    ## Should find the variant
    checkTrue(length(result) > 0L)
    ## LOCSTART/LOCEND should NOT be NA (normal mapping works)
    checkTrue(!any(is.na(mcols(result)$LOCSTART)))
}

test_predictCoding_large_indel_warning <- function()
{
    ## A large deletion spanning multiple exons should emit a warning
    ## from .localCoordinates, not silently return empty.
    query <- GRanges("chr1", IRanges(start = 1000, width = 2543))

    warned <- FALSE
    result <- withCallingHandlers(
        VariantAnnotation:::.localCoordinates(query, cdsbytx,
                                              ignore.strand = TRUE),
        warning = function(w) {
            if (grepl("span multiple exons", conditionMessage(w))) {
                warned <<- TRUE
                invokeRestart("muffleWarning")
            }
        }
    )

    ## A warning should have been emitted
    checkTrue(warned)
}

test_locateVariants_mixed_variants <- function()
{
    ## Mix of a small variant (fits in one exon) and a large spanning deletion.
    ## Both should appear in results.
    query <- GRanges("chr1",
                     IRanges(start = c(1050, 1000), width = c(10, 2543)))

    result <- suppressWarnings(
        VariantAnnotation:::.makeResult(query, cdsbytx, "coding",
                                        ignore.strand = TRUE, asHits = FALSE)
    )

    ## Both variants should be in results
    checkTrue(1L %in% mcols(result)$QUERYID)
    checkTrue(2L %in% mcols(result)$QUERYID)
    checkIdentical(length(unique(mcols(result)$QUERYID)), 2L)
}
