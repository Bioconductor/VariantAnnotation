### =========================================================================
### VariantType class methods 
### =========================================================================

## 'show' methods

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
        cat("upstream:", upstream(object), "\n")
        cat("downstream:", downstream(object), "\n")
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

## Classes with constructors only

CodingVariants <- function() new("CodingVariants")

IntronVariants <- function() new("IntronVariants")

IntergenicVariants <- function() new("IntergenicVariants")

UTRVariants <- function() new("UTRVariants")

ThreeUTRVariants <- function() new("ThreeUTRVariants")

FiveUTRVariants <- function() new("FiveUTRVariants")

SpliceSiteVariants <- function() new("SpliceSiteVariants")

## Classes with getters and setters

PromoterVariants <- function(upstream=2000, downstream=200)
{
    if (!isSingleNumber(upstream))
        stop("'upstream' must be a single integer")
    if (!is.integer(upstream))
        upstream <- as.integer(upstream)
    if (!isSingleNumber(downstream))
        stop("'downstream' must be a single integer")
    if (!is.integer(downstream))
        downstream <- as.integer(downstream)
    new("PromoterVariants", upstream=upstream, 
        downstream=downstream) 
}

setMethod("upstream", "PromoterVariants",
    function(x)
{
    slot(x, "upstream")
})

setReplaceMethod("upstream", "PromoterVariants",
    function(x, value)
{
    slot(x, "upstream") <- value
    x
})

setMethod("downstream", "PromoterVariants",
    function(x)
{
    slot(x, "downstream")
})

setReplaceMethod("downstream", "PromoterVariants",
    function(x, value)
{
    slot(x, "downstream") <- value
    x
})

AllVariants <- function(upstream=2000L, downstream=200L)
{
    if (!isSingleNumber(upstream))
        stop("'upstream' must be a single integer")
    if (!is.integer(upstream))
        upstream <- as.integer(upstream)
    if (!isSingleNumber(downstream))
        stop("'downstream' must be a single integer")
    if (!is.integer(downstream))
        downstream <- as.integer(downstream)
    new("AllVariants", upstream=upstream, 
        downstream=downstream) 
}

setMethod("upstream", "AllVariants",
    function(x)
{
    slot(x, "upstream")
})

setReplaceMethod("upstream", "AllVariants",
    function(x, value)
{
    slot(x, "upstream") <- value
    x
})

setMethod("downstream", "AllVariants",
    function(x)
{
    slot(x, "downstream")
})

setReplaceMethod("downstream", "AllVariants",
    function(x, value)
{
    slot(x, "downstream") <- value
    x
})

