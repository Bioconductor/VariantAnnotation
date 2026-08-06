# False Positive Example

## What this is

A **false positive PR** is one where an automated agent "fixes" something that
was never actually broken.  This directory contains a labelled example so that
human reviewers have a calibration reference: *this is what a false positive
looks like*.

## The fabricated "bug"

The agent noticed that `isSNV()` returns `FALSE` for variants where `REF` and
`ALT` are the same nucleotide, and filed a PR titled:

> *"fix: isSNV() incorrectly classifies REF==ALT variants as non-SNV"*

The proposed change added a special-case:

```r
# PROPOSED (wrong) change in R/methods-isSNV.R
isSNV <- function(x, ...) {
    ref <- ref(x); alt <- unlist(alt(x))
    # NEW: treat REF==ALT as SNV
    same <- as.character(ref) == as.character(alt)
    nchar(ref) == 1L & nchar(alt) == 1L | same
}
```

## Why it is wrong

`REF == ALT` is a **monomorphic site** — by definition not a variant at all.
Every VCF specification and bioinformatics tool agrees: a SNV requires
`REF != ALT`.  The existing behaviour (`isSNV` returns FALSE) is correct.

The agent was fooled because:
1. It found a user report saying "isSNV returns FALSE for my variants" without
   reading that the user's file had malformed REF==ALT entries.
2. It fixed the symptom (the return value) without asking whether the input
   was valid.

## Hedgehog property that catches it

```r
test_that("isSNV: REF==ALT is never a SNV", {
    forall(gen.element(c("A","C","G","T")), function(b) {
        vcf <- make_vcf(ref = b, alt = b)   # monomorphic
        expect_false(isSNV(vcf))
    })
})
```

This property runs 100 random monomorphic variants and fails immediately on
the proposed change.  The existing unit tests did not cover this case because
they only tested clearly-variant inputs.

## Lesson

Property P: *"isSNV(x) is TRUE only when REF != ALT and both have width 1"*
is a logical invariant derivable from the VCF spec alone — no domain knowledge
of the specific codebase needed.  Automated PRs that violate spec-level
invariants are almost always false positives.
