### =========================================================================
### VariantType class methods 
### =========================================================================

## 'show' method

setMethod("show", "VariantType",
    function(object) 
    {
        cat("class:", class(object), "\n")
    }
)

## constructors

AllVariants <- function() new("AllVariants")

CodingVariants <- function() new("CodingVariants")

IntronVariants <- function() new("IntronVariants")

IntergenicVariants <- function() new("IntergenicVariants")

UTRVariants <- function() new("UTRVariants")

ThreeUTRVariants <- function() new("ThreeUTRVariants")

FiveUTRVariants <- function() new("FiveUTRVariants")

SpliceSiteVariants <- function() new("SpliceSiteVariants")
