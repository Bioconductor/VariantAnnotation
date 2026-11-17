# False Negative Example

## What this is

A **false negative PR** is one where an automated agent correctly identifies a
real bug but its fix is wrong or incomplete — it patches the symptom without
curing the disease, or introduces a new bug while fixing the original.

## The real bug

`expand()` on a `CollapsedVCF` with multi-allelic sites was silently dropping
INFO fields whose `Number` header tag was `"A"` (one value per ALT allele)
when those fields contained `NA` for some alleles.  Issue #79.

## The agent's (incomplete) fix

The agent found the right function (`expand,CollapsedVCF`) and added:

```r
# PROPOSED (incomplete) fix
info_A_fields <- names(info(header(x)))[
    info(header(x))$Number == "A"]
for (fld in info_A_fields) {
    if (anyNA(info(x)[[fld]]))
        info(x)[[fld]][is.na(info(x)[[fld]])] <- NA  # no-op!
}
```

The fix is a no-op — it assigns NA to already-NA positions and never actually
expands the per-allele INFO values alongside the ALT column.  All existing
unit tests still pass because they used complete (non-NA) INFO fields.

## Why the existing tests didn't catch it

The unit tests in `test_expand-methods.R` built VCFs with fully-populated INFO
fields.  The hedgehog property below generates VCFs with random NA patterns
and immediately finds the failure:

```r
test_that("P: expand() preserves INFO/Number=A length relative to ALT", {
    forall(gen_multiallelic_vcf(), function(vcf) {
        ex <- expand(vcf)
        ## After expansion each row has exactly 1 ALT;
        ## every Number=A INFO field must have length 1 (not NA-collapsed)
        a_fields <- names(info(header(ex)))[info(header(ex))$Number == "A"]
        for (fld in a_fields) {
            vals <- info(ex)[[fld]]
            expect_equal(length(vals), nrow(ex))
        }
    })
})
```

## The correct fix

The correct fix (merged in PR #XXX) uses `unlist()` on the `List` column so
that each expanded row gets the corresponding allele's value, and explicitly
propagates `NA` for missing entries without collapsing them.

## Lesson

The agent correctly diagnosed *which* function was broken but wrote a patch
based on a surface-level reading of the stack trace.  A property test exposing
the *output contract* — "expanded VCF has one INFO value per row for Number=A
fields" — would have immediately shown the fix was wrong, because the no-op
patch still fails 100/100 random draws with NA-containing inputs.

This illustrates why property-based testing is more valuable for evaluating
automated PRs than re-running the original unit tests: the agent's fix was
specifically tuned to pass those exact tests.
