### =========================================================================
### vcfHeader class methods 
### =========================================================================


## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Constructor 
##

VCFHeader <-
    function(reference=character(), samples=character(), 
             header=DataFrameList(), ...)
{
    new("VCFHeader", reference=reference, samples=samples, header=header, ...) 
}
 
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Getters and Setters
##

setMethod("reference", "VCFHeader", 
    function(x) 
{
    slot(x, "reference")
})

setMethod("samples", "VCFHeader", 
    function(object) 
{
    slot(object, "samples")
})

setMethod("header", "VCFHeader", 
    function(x) 
{
    slot(x, "header")
})

## meta
setMethod("meta", "VCFHeader", 
    function(x) 
{
    slot(x, "header")$META
})

## fixed (QUAL, FILTER, ALT, REF)
setMethod("fixed", "VCFHeader", 
    function(x) 
{
    fixed <- c("REF", "ALT", "QUAL", "FILTER") 
    lst <- slot(x, "header")
    lst[names(lst) %in% fixed]
})

## info 
setMethod("info", "VCFHeader", 
    function(x) 
{
    info <- slot(x, "header")$INFO
    if (is.null(info))
        info <- DataFrame()
    info
})

## geno
setMethod("geno", "VCFHeader",
    function(x)
{
    geno <- slot(x, "header")$FORMAT
    if (is.null(geno))
        geno <- DataFrame()
    geno 
})

setMethod("seqinfo", "VCFHeader",
function(x)
{
  contig <- slot(x, "header")$contig
  if (is.null(contig))
    Seqinfo()
  else Seqinfo(rownames(contig),
               seqlengths =
               if (is.null(contig$length)) NA else as.integer(contig$length),
               genome = if (is.null(contig$assembly)) NA else contig$assembly)
})

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Show
##

setMethod(show, "VCFHeader",
    function(object)
{
    selectSome <- BiocGenerics:::selectSome
    scat <- function(fmt, vals=character(), exdent=2, ...)
    {
        vals <- ifelse(nzchar(vals), vals, "''")
        lbls <- paste(selectSome(vals), collapse=" ")
        txt <- sprintf(fmt, length(vals), lbls)
        cat(strwrap(txt, exdent=exdent, ...), sep="\n")
    }
    cat("class:", class(object), "\n")

    samples <- samples(object) 
    scat("samples(%d): %s\n", samples)

    meta <- rownames(meta(object)) 
    scat("meta(%d): %s\n", meta)

    fixed <- names(fixed(object)) 
    scat("fixed(%d): %s\n", fixed)

    info <- rownames(info(object)) 
    scat("info(%d): %s\n", info)

    geno <- rownames(geno(object))
    scat("geno(%d): %s\n", geno)
})
