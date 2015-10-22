### =========================================================================
### All classes 
### =========================================================================


### -------------------------------------------------------------------------
### VCF (VIRTUAL) 
###

setClass("VCF",
    contains=c("VIRTUAL", "RangedSummarizedExperiment"),
    representation( 
        fixed="DataFrame",
        info="DataFrame")
)

.valid.VCF.fixed <- function(object)
{
    xlen <- dim(object)[1]
    ffld <- slot(object, "fixed")
    nms <- colnames(ffld)

    if (nrow(ffld) != 0) {
        if (nrow(ffld) != xlen)
            return(paste("'fixed(object)' and 'rowRanges(object) must have the same ",
                   "number of rows", sep=""))
        if (!all(nms %in% c("REF", "ALT", "QUAL", "FILTER")))
            return(paste("'fixed(object)' colnames must be ",
                   "'REF', 'ALT', 'QUAL' and 'FILTER'", sep=""))
        if ("REF" %in% nms)
            if (!is(ffld$REF, "DNAStringSet"))
                return("'ref(object)' must be a DNAStringSet")
        if ("QUAL" %in% nms)
            if (!is(ffld$QUAL, "numeric"))
                return("'qual(object)' must be numeric")
        if ("FILTER" %in% nms)
            if (!is(ffld$FILTER, "character"))
                return("'filt(object)' must be a character")
    }
    NULL
}

.valid.VCF.info <- function(object)
{
    xlen <- nrow(object)
    i <- slot(object, "info")
    if (nrow(i) != xlen)
        return("'info' must have the same number of rows as 'rowRanges'")
    NULL
}

.valid.VCF.alt <- function(object)
{
    if (length(alt <- alt(object))) {
        if (is(object, "ExpandedVCF")) {
            if (!is(alt, "DNAStringSet") && !is.character(alt))
                return(paste("'alt(object)' must be a DNAStringSet or a ",
                       "character", sep=""))
        } else if (is(object, "CollapsedVCF")) {
            if (!is(alt, "DNAStringSetList") && !is(alt, "CharacterList"))
                return(paste("'alt(object)' must be a DNAStringSetList or a ",
                       "CharacterList", sep=""))
        }
    }
    NULL
}

.valid.VCF <- function(object)
{
    c(.valid.VCF.fixed(object),
      .valid.VCF.info(object),
      .valid.VCF.alt(object))
}

### -------------------------------------------------------------------------
### CollapsedVCF 
###

setClass("CollapsedVCF", 
    contains="VCF",
    prototype=prototype(
      fixed=DataFrame(REF=DNAStringSet(), ALT=DNAStringSetList(),
                      QUAL=numeric(), FILTER=character())),
    validity=.valid.VCF)

### Automatically generated "coerce<-" method is broken so we fix it.
### See S4Vectors/R/S4-utils.R in the S4Vectors package for more information.
S4Vectors:::setReplaceAs("CollapsedVCF", "RangedSummarizedExperiment",
    S4Vectors:::canonical_replace_as_2
)

### -------------------------------------------------------------------------
### ExpandedVCF 
###

setClass("ExpandedVCF", 
    contains="VCF",
    prototype=prototype(
      fixed=DataFrame(REF=DNAStringSet(), ALT=DNAStringSet(),
                      QUAL=numeric(), FILTER=character())),
    validity=.valid.VCF)

### Automatically generated "coerce<-" method is broken so we fix it.
### See S4Vectors/R/S4-utils.R in the S4Vectors package for more information.
S4Vectors:::setReplaceAs("ExpandedVCF", "RangedSummarizedExperiment",
    S4Vectors:::canonical_replace_as_2
)

### -------------------------------------------------------------------------
### VCFHeader
###

setClass("VCFHeader",
    representation(
        reference="character",
        samples="character",
        header="SimpleDataFrameList"
    )
)

.valid.VCFHeader.colnames <- function(value, slotname)
{
    if (length(value)) {
        col <- c("Number", "Type", "Description")
        if (ncol(value) != 3 || !all(names(value) %in% col))
            return(paste0("'", slotname, "(VCFHeader)' must be a ",
                   "3 column DataFrame with names ", 
                   paste(col, collapse=", ")))
        else
            return(NULL)
    }
    NULL
}

.valid.VCFHeader.info <- function(object)
    .valid.VCFHeader.colnames(info(object), "info")

.valid.VCFHeader.geno <- function(object)
    .valid.VCFHeader.colnames(geno(object), "geno")

.valid.VCFHeader.meta <- function(object)
{
    msg <- NULL
    len <- sapply(names(meta(object)), nchar)
    if (any(len == 0))
        msg <- "all elements of 'meta(VCFHeader)' must be named"

    msg
}

## Test VCFHeader only 
.valid.VCFHeader <- function(object)
{
    return(c(.valid.VCFHeader.info(object),
             .valid.VCFHeader.geno(object),
             .valid.VCFHeader.meta(object)))
}
setValidity("VCFHeader", .valid.VCFHeader, where=asNamespace("VariantAnnotation"))

## These validity checks test the fields in the VCF against 
## the VCFHeader. Called from header<-,VCF-method. 
.valid.VCFHeadervsVCF.fields <- function(object, slotname)
{
    diff <- setdiff(names(slotname(object)), rownames(slotname(header(object))))
    if (length(diff)) {
        warning(paste0(attributes(slotname)$generic[1],
                " fields with no header: ", paste(diff, collapse=",")), 
                call.=FALSE)
    }
    NULL
}

.valid.VCFHeadervsVCF <- function(object)
{
    validObject(header(object))
    c(.valid.VCFHeadervsVCF.fields(object, info),
      .valid.VCFHeadervsVCF.fields(object, geno))
} 

### -------------------------------------------------------------------------
### VcfFile, VcfFileList
###

.VcfFile = setRefClass("VcfFile", contains="TabixFile") 


setClass("VcfFileList", contains="TabixFileList",
         prototype=prototype(elementType="VcfFile"))

### -------------------------------------------------------------------------
### ScanVcfParam
###

setClass("ScanVcfParam", contains="ScanBVcfParam")

### -------------------------------------------------------------------------
### SIFT, PROVEAN and PolyPhen
###

.SIFTDb <- setRefClass("SIFTDb", contains = "AnnotationDb")

.PROVEANDb <- setRefClass("PROVEANDb", contains = "AnnotationDb")

.PolyPhenDb <- setRefClass("PolyPhenDb", contains = "AnnotationDb")

### -------------------------------------------------------------------------
### VariantType class hierarchy
###

setClass("VariantType",
    representation(
        "VIRTUAL"
    )
)

setMethod("show", "VariantType",
    function(object)
    {
        cat("class:", class(object), "\n")
    }
)

setClass("CodingVariants", contains="VariantType")

setClass("IntronVariants", contains="VariantType")

setClass("ThreeUTRVariants", contains="VariantType")

setClass("FiveUTRVariants", contains="VariantType")

setClass("SpliceSiteVariants", contains="VariantType")

setClass("IntergenicVariants",
    contains="VariantType",
    representation(upstream="integer",
                   downstream="integer")
)

setClass("PromoterVariants",
    contains="VariantType",
    representation(upstream="integer",
                   downstream="integer")
)

setClass("AllVariants",
    contains="VariantType",
    representation(promoter="PromoterVariants",
                   intergenic="IntergenicVariants")
)

### -------------------------------------------------------------------------
### VRanges
###

### This is useful when variants are the unit of analysis, especially
### when diagnosis the effects of filters. Thus, this corresponds to
### an expansion of VCF, where each variant only has one alt
### allele. If the alt is NA, it describes a reference/WT call.
###

.valid.VRanges.ref <- function(object) {
  c(if (any(is.na(object@ref)))
      "'ref' must not contain 'NA' values",
    if (length(which(object@ref == object@alt)) > 0L)
      "'ref' must not match 'alt'")
}

.valid.VRanges.alt <- function(object) {
  if (any(is.na(object@alt) & !is.na(object@altDepth)))
    paste("if 'alt' is 'NA', then 'altDepth' should be 'NA'")
}

.valid.VRanges.strand <- function(object) {
  if (any(object@strand == "-"))
    paste("'strand' must always be '+' or '*'")
}

naToZero <- function(x) {
    x[is.na(x)] <- 0L
    x
}

.valid.VRanges.depth <- function(object) {
  checkDepth <- function(name) {
    if (any(slot(object, name) < 0, na.rm = TRUE))
      paste0("'", name, "' must be non-negative")
  }
  depth.slots <- c("refDepth", "altDepth", "totalDepth")
  do.call(c, lapply(depth.slots, checkDepth))
}

### FIXME: we are not yet checking for redundant observations, i.e.,
### rows with the same seqname, start, end, alt, and sample.
.valid.VRanges <- function(object)
{
  c(.valid.VRanges.ref(object),
    .valid.VRanges.alt(object),
    .valid.VRanges.strand(object),
    .valid.VRanges.depth(object))
}

setTypedRle <- function(type) {
  cname <- paste0(type, "Rle")
  setClass(cname, representation(values = type), contains = "Rle")
  setValidity(cname,
              function(object) {
                if (!is(runValue(object), type))
                  paste("run values in Rle must be", type)
              })
  coercer <- get(paste0("as.", type)) # as(from, type) is NOT the same
  setAs("Rle", cname, function(from) {
    if (!is(runValue(from), type))
      slot(from, "values", check = FALSE) <- coercer(runValue(from))
    new(cname, from)
  })
  setAs("vectorORfactor", cname, function(from) {
    new(cname, Rle(coercer(from)))
  })
  setMethod("[", cname,
            function(x, i, j, ..., drop = getOption("dropRle", FALSE)) {
              as(callNextMethod(), cname)
            })
  setReplaceMethod("[", cname,
                   function(x, i, j, ..., value) {
                     if (missing(i))
                       ans <- callGeneric(x = as(x, "Rle"), value = value)
                     else
                       ans <- callGeneric(x = as(x, "Rle"), i = i, value = value)
                     as(ans, paste0(class(runValue(ans)), "Rle"))
                   })
  setMethod("window", cname,
            function(x, ...) {
              as(callNextMethod(), cname)
            })
  setMethod("c", cname,
            function(x, ...) {
              ### FIXME: cannot use callNextMethod, because broken for
              ### primitive generics where the actual primitive object
              ### is part of the call, e.g., do.call(generic, args)
              as(do.call(c, lapply(list(x, ...), as, "Rle")), cname)
            })
  cname
}

### TODO: These could land in IRanges at some point. They are useful
### for constraining class representations, and are convenient
### coercion targets. In theory, we could also use them more broadly,
### i.e., for Rle method dispatch, but that would require work in
### Rle() and runValue<- first, and would break serialized objects.

setTypedRle("raw")
setTypedRle("logical")
setTypedRle("integer")
setTypedRle("numeric")
setTypedRle("complex")
setTypedRle("character")
setTypedRle("factor")

setAtomicRleUnion <- function(type) {
  cname <- paste0(type, "OrRle")
  rlename <- paste0(type, "Rle")
  setClassUnion(cname, c(type, rlename))
  coercer <- get(paste0("as.", type)) # as(from, type) is NOT the same
  setAs("ANY", cname, function(from) {
    coercer(from)
  })
  setAs("Rle", cname, function(from) {
    as(from, rlename)
  })
  cname
}

.integerOrRle <- setAtomicRleUnion("integer")
.characterOrRle <- setAtomicRleUnion("character")
.factorOrRle <- setAtomicRleUnion("factor")

setClass("VRanges",
         representation(ref = "character", # or XStringSet?
                        alt = .characterOrRle,
                        totalDepth = .integerOrRle,
                        refDepth = .integerOrRle,
                        altDepth = .integerOrRle,
                        sampleNames = .factorOrRle,
                        softFilterMatrix = "matrix",
                        hardFilters = "FilterRules"),
         contains = "GRanges",
         validity = .valid.VRanges)

### -------------------------------------------------------------------------
### VRangesList
###

setClass("VRangesList", representation("VIRTUAL"),
         prototype = prototype(elementType = "VRanges"),
         contains = "List")

setClass("SimpleVRangesList",
         contains = c("VRangesList", "SimpleGenomicRangesList"))

setClass("CompressedVRangesList",
         representation(elementMetadata = "DataFrame"),
         prototype = prototype(unlistData = new("VRanges")),
         contains = c("VRangesList", "GRangesList"))
