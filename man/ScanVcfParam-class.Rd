\name{ScanVcfParam-class}
\Rdversion{1.1}
\docType{class}

\alias{ScanVcfParam}
\alias{ScanVcfParam-class}
\alias{ScanVcfParam,missing-method}
\alias{ScanVcfParam,ANY-method}
\alias{vcfFixed}
\alias{vcfInfo}
\alias{vcfGeno}
\alias{vcfTrimEmpty}
\alias{vcfWhich}
\alias{vcfAsGRanges}

\title{Parameters for scanning VCF files}

\description{
  Use \code{ScanVcfParam()} to create a parameter object influencing 
  which records and fields are imported from a VCF file. Record
  parsing is based on genomic coordinates and requires a Tabix index
  file. Individual elements can be selected from the \sQuote{INFO} and 
  \sQuote{GENO} fields.
}

\usage{
  ScanVcfParam(fixed=character(), info=character(), geno=character(), 
               trimEmpty=TRUE, which, asGRanges=FALSE, ...)

  ## Accessors
 
  vcfFixed(object)
  vcfInfo(object)
  vcfGeno(object)
  vcfTrimEmpty(object)
  vcfWhich(object)
  vcfAsGRanges(object)
}

\arguments{
  \item{info}{A character() vector of \sQuote{INFO} fields (see
    \code{\link{scanVcfHeader}}) to be returned. If not specified,
    all fields are returned; if \code{info=NA} no fields 
    are returned.
  }
  \item{geno}{A character() vector of \sQuote{GENO} fields (see
    \code{\link{scanVcfHeader}}) to be returned. If not specified,
    all fields are returned; if \code{geno=NA} no fields 
    are returned.
  }
  \item{trimEmpty}{A logical(1) indicating whether \sQuote{GENO} fields
    with no values should be returned.
  }
  \item{which}{An object, for which a method is defined (see usage,
    above), describing the sequences and ranges to be queried. Variants
    whose \code{POS} lies in the interval(s) \code{[start, end]} are
    returned. 
  }
  \item{asGRanges}{A logical(1) indicating whether a long form (expanded) 
    \code{GRanges} should be returned instead of a \code{VCF} object. For 
    use with \code{readVcf} only; ignored by \code{scanVcf}.

    The return value is a \code{GRanges} expanded to match the dimensions
    of the unlisted elements specified in \sQuote{info} or \sQuote{geno}.
    The \sQuote{info} or \sQuote{geno} elements are included as elementMetadata
    columns. When \code{asGRanges=TRUE} multiple elements can be provided to 
    the \sQuote{info} argument but only a single element can be provided to 
    \sQuote{geno}. This restriction is due to the dimension of the \sQuote{geno} 
    data.
  }
  \item{object}{An instance of class \code{ScanVcfParam}.
  }
  \item{\dots}{Arguments passed to methods.
  }
  \item{fixed}{Internal use only.
  }
}

\section{Objects from the Class}{

  Objects can be created by calls of the form \code{ScanVcfParam()}.

}
\section{Slots}{
  \describe{
    \item{\code{which}:}{Object of class \code{"RangesList"} indicating
      which reference sequence and coordinate variants must overlap.
    }
    \item{\code{info}:}{Object of class \code{"character"} indicating
      portions of \sQuote{INFO} to be returned.
    }
    \item{\code{geno}:}{Object of class \code{"character"} indicating
      portions of \sQuote{GENO} to be returned. 
    }
    \item{\code{trimEmpty}:}{Object of class \code{"logical"} indicating
      whether empty \sQuote{GENO} fields are to be returned.
    }
    \item{\code{asGRanges}:}{Object of class \code{"character"} indicating
      the return of a long form \code{GRanges} with either
      the \sQuote{INFO} or \sQuote{GENO} fields expanded in the elementMetadata.
    }
    \item{\code{fixed}:}{Internal use only.
    }
  }
}

\section{Functions and methods}{

  See 'Usage' for details on invocation.

  Constructor:
  \describe{
    \item{ScanVcfParam:}{Returns a \code{ScanVcfParam} object. 
      The \code{which} argument to the constructor can be one of several types, 
      as documented above.}  
  }

  Accessors:
  \describe{
    \item{vcfFixed, vcfInfo, vcfGeno, vcfTrimEmpty, vcfWhich, vcfAsGRanges:}{
      Return the corresponding field from \code{object}.
    }
  }

  Methods:
  \describe{
    \item{show}{Compactly display the object.
    }
  }
}
\author{
  Martin Morgan and Valerie Obenchain <vobencha@fhcrc.org>
}
\seealso{
  \code{\link{readVcf}}
}

\examples{
  ScanVcfParam()

  ## -----------------------------
  ## 'which' argument
  ## -----------------------------
  ## To subset on genomic coordinates, create a GRanges, RangedData or 
  ## RangesList with the ranges of interest. This object is supplied
  ## to ScanVcfParam() as the 'which' argument.
  which <- RangesList(seq1=IRanges(1000, 2000), 
      seq2=IRanges(c(100, 1000), c(1000, 2000)))
  ScanVcfParam(which=which)
 
  ## -----------------------------
  ## 'info' and 'geno' arguments
  ## -----------------------------
  ## Subsets of INFO or FORMAT field elements are specified in the 'info'
  ## and 'geno' arguments. This call specifies the return of the "GT" 'geno' 
  ## field and the subset of ranges in 'which'. All 'info' fields will be 
  ## returned. 
  ScanVcfParam(geno="GT", which=which)
 
  ## This param specifies the geno 'GT' and none of the 'info' fields.
  ScanVcfParam(geno="GT", info=NA) 

  ## -----------------------------
  ## 'asGRanges' argument 
  ## -----------------------------
  ## When asGRanges=TRUE only one of 'info' or 'geno' can be provided.
  ## Specifying values for both is wrong and will generate an error :
  ## ScanVcfParam(info="DP", geno="GT", asGranges=TRUE)

  ## When asGRanges=TRUE only a single element can be provided to 'geno' 
  ScanVcfParam(geno="GT", asGranges=TRUE)
  ## but multiple can be provided to 'info'.
  ScanVcfParam(info=c("DP", "NS"), asGranges=TRUE)
} 
\keyword{classes}
