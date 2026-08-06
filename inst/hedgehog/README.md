# VariantAnnotation — Hedgehog Property-Based Tests

This directory contains [hedgehog](https://cran.r-project.org/package=hedgehog)
property-based tests for `VariantAnnotation`.  They complement the
example-based unit tests in `inst/unitTests/` by testing **invariants** —
properties that must hold for *any* valid input, not just the handful of
hand-crafted cases covered by unit tests.

## Why property-based testing?

Example-based tests check that `f(a) == b` for specific `a`.  Property-based
tests check that `∀ a ∈ A: P(f(a))` — a logical invariant over the whole input
space.  This matters for evaluating automated PRs: a patch can be tuned to pass
6 unit tests while still failing a hedgehog property in 3 of 100 random draws.

The `hedgehog` package provides:
- **Integrated shrinking** — when a property fails, hedgehog automatically
  shrinks the failing input to the smallest possible counterexample, making
  bugs much easier to diagnose.
- **testthat integration** — properties run inside `test_that()` blocks and
  appear in standard R CMD check output.

## Files

| File | Properties |
|------|-----------|
| `prop_readVcf_roundtrip.R` | Row count, rownames, seqnames, ALT alleles, expand row count all survive a write→read round-trip |
| `prop_predictCoding_invariants.R` | CONSEQUENCE is always a valid value; REFAA/VARAA are always AAStringSet; CDSLOC is IRanges; output ≤ input; REFCODON width is a multiple of 3 |
| `false_positive/README.md` | Annotated example of an AI-generated PR that "fixes" correct behaviour |
| `false_negative/README.md` | Annotated example of an AI-generated PR that correctly finds a bug but applies the wrong fix |

## Running

```r
library(VariantAnnotation)
library(hedgehog)
library(testthat)

source(system.file("hedgehog/prop_readVcf_roundtrip.R",
                   package = "VariantAnnotation"))
source(system.file("hedgehog/prop_predictCoding_invariants.R",
                   package = "VariantAnnotation"))
```

Or run all properties at once:

```r
hedgehog_dir <- system.file("hedgehog", package = "VariantAnnotation")
for (f in list.files(hedgehog_dir, pattern = "^prop_.*\\.R$",
                     full.names = TRUE))
    source(f)
```

## Adding new properties

A good property for `VariantAnnotation` answers one of:
1. **Round-trip**: does writing and re-reading preserve information?
2. **Taxonomy**: is every output value drawn from a known finite set?
3. **Type contract**: does the output always have the right S4 class?
4. **Monotonicity**: does `expand()` always produce ≥ as many rows as input?
5. **Spec compliance**: does the output satisfy a VCF spec requirement?

Properties that come from the VCF specification are particularly valuable
because they are independent of the implementation — an automated agent cannot
accidentally tune a fix to pass them without actually being correct.

## Relationship to automated PR review (issue #113)

The `false_positive/` and `false_negative/` directories provide labelled
calibration examples for evaluating AI-generated PRs.  Before approving any
automated PR, run the relevant property tests: a PR that passes the unit tests
but fails a property test is almost certainly a false positive or false
negative.
