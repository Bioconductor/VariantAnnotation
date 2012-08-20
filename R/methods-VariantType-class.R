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

setMethod("show", "FlankingVariants",
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

FlankingVariants <- function(upstream = 500, downstream = 500)
{
    if (any((upstream < 1)  | (downstream < 1)))
        stop("'upstream' and 'downstream' must be integers > 0")
    new("FlankingVariants", upstream=as.integer(upstream), 
        downstream=as.integer(downstream))
}

setMethod("upstream", "FlankingVariants",
    function(x)
{
    slot(x, "upstream")
})

setReplaceMethod("upstream", "FlankingVariants",
    function(x, value)
{
    slot(x, "upstream") <- value
    x
})

setMethod("downstream", "FlankingVariants",
    function(x)
{
    slot(x, "downstream")
})

setReplaceMethod("downstream", "FlankingVariants",
    function(x, value)
{
    slot(x, "downstream") <- value
    x
})

AllVariants <- function(upstream = 500, downstream = 500)
{
    if (any((upstream < 1)  | (downstream < 1)))
        stop("'upstream' and 'downstream' must be integers > 0")
    new("AllVariants", upstream=as.integer(upstream), 
        downstream=as.integer(downstream))
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

