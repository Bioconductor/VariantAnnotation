### =========================================================================
### VariantType class methods 
### =========================================================================

### 'show' methods

setMethod("show", "VariantType",
    function(object) 
    {
        cat("class:", class(object), "\n")
    }
)

setMethod("show", "AllVariants",
    function(object) 
    {
        cat("class:", class(object), "\n")
        cat("promoter: \n")
        cat("  upstream:", upstream(promoter(object)), "\n")
        cat("  downstream:", downstream(promoter(object)), "\n")
        cat("intergenic: \n")
        cat("  upstream:", upstream(intergenic(object)), "\n")
        cat("  downstream:", downstream(intergenic(object)), "\n")
        cat("  idType:", idType(intergenic(object)), "\n")
    }
)

setMethod("show", "PromoterVariants",
    function(object) 
    {
        cat("class:", class(object), "\n")
        cat("upstream:", upstream(object), "\n")
        cat("downstream:", downstream(object), "\n")
    }
)

setMethod("show", "IntergenicVariants",
    function(object) 
    {
        cat("class:", class(object), "\n")
        cat("upstream:", upstream(object), "\n")
        cat("downstream:", downstream(object), "\n")
        cat("idType:", idType(object), "\n")
    }
)

### constructors

CodingVariants <- function() new("CodingVariants")

IntronVariants <- function() new("IntronVariants")

UTRVariants <- function() new("UTRVariants")

ThreeUTRVariants <- function() new("ThreeUTRVariants")

FiveUTRVariants <- function() new("FiveUTRVariants")

SpliceSiteVariants <- function() new("SpliceSiteVariants")

### .checkArgs() takes the place of a validity method.
### 'upstream' and 'downstream' are forced to integers
### which changes the object and can't be done a validity
### method. 
.checkArgs <- function(x, argName)
{
    if (!isSingleNumber(x) | x < 0)
        stop(paste0("'", argName, "'", " must be a single integer >= 0"))
    if (!is.integer(x))
        as.integer(x)
    else
        x
}

IntergenicVariants <- function(upstream=1e+06L, downstream=1e+06L,
                               idType=c("gene", "tx"))
{
    upstream <- .checkArgs(upstream, "upstream") 
    downstream <- .checkArgs(downstream, "downstream")
    if (missing(idType)) {
        idType <- match.arg(idType)
    } else {
        idType <- tolower(idType)
        if (!idType %in% c("gene", "tx")) 
            stop("'idType' must be one of 'gene' or 'tx'")
    }
    new("IntergenicVariants", upstream=upstream, 
        downstream=downstream, idType=idType) 
}

PromoterVariants <- function(upstream=2000L, downstream=200L)
{
    upstream <- .checkArgs(upstream, "upstream") 
    downstream <- .checkArgs(downstream, "downstream") 
    new("PromoterVariants", upstream=upstream, 
        downstream=downstream) 
}

AllVariants <- function(promoter=PromoterVariants(),
                        intergenic=IntergenicVariants())
{
    new("AllVariants", promoter=promoter, intergenic=intergenic)
}

### getters and setters

setMethod("upstream", "PromoterVariants",
    function(x) slot(x, "upstream"))

setReplaceMethod("upstream", "PromoterVariants",
    function(x, value)
{
    slot(x, "upstream") <- .checkArgs(value, "upstream") 
    x
})

setMethod("downstream", "PromoterVariants",
    function(x) slot(x, "downstream"))

setReplaceMethod("downstream", "PromoterVariants",
    function(x, value)
{
    slot(x, "downstream") <- .checkArgs(value, "downstream") 
    x
})

setMethod("upstream", "IntergenicVariants",
    function(x) slot(x, "upstream"))

setReplaceMethod("upstream", "IntergenicVariants",
    function(x, value)
{
    slot(x, "upstream") <- .checkArgs(value, "upstream") 
    x
})

setMethod("downstream", "IntergenicVariants",
    function(x) slot(x, "downstream"))

setReplaceMethod("downstream", "IntergenicVariants",
    function(x, value)
{
    slot(x, "downstream") <- .checkArgs(value, "downstream") 
    x
})

setMethod("idType", "IntergenicVariants",
    function(x) slot(x, "idType"))

setReplaceMethod("idType", "IntergenicVariants",
    function(x, value)
{
    if (!value %in% c("gene", "tx"))
        stop("'idType' must be one of 'gene' or 'tx'")
    slot(x, "idType") <- value 
    x
})

setMethod("downstream", "IntergenicVariants",
    function(x) slot(x, "downstream"))

setReplaceMethod("downstream", "IntergenicVariants",
    function(x, value)
{
    slot(x, "downstream") <- .checkArgs(value, "downstream") 
    x
})

setMethod("promoter", "AllVariants",
    function(x) slot(x, "promoter"))

setReplaceMethod("promoter", "AllVariants",
    function(x, value)
{
    if (class(value) != "PromoterVariants")
        stop("'value' must be a 'PromoterVariants' object")
    slot(x, "promoter") <- value 
    x
})

setMethod("intergenic", "AllVariants",
    function(x) slot(x, "intergenic"))

setReplaceMethod("intergenic", "AllVariants",
    function(x, value)
{
    if (class(value) != "IntergenicVariants")
        stop("'value' must be a 'IntergenicVariants' object")
    slot(x, "intergenic") <- value 
    x
})
