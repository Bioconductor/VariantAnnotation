\name{ScanVcfParam-class}
\docType{class}

\alias{ScanVcfParam}
\alias{ScanVcfParam-class}
\alias{ScanVcfParam,missing-method}
\alias{ScanVcfParam,ANY-method}
\alias{vcfFixed}
\alias{vcfFixed<-}
\alias{vcfInfo}
\alias{vcfInfo<-}
\alias{vcfGeno}
\alias{vcfGeno<-}
\alias{vcfTrimEmpty}
\alias{vcfTrimEmpty<-}
\alias{vcfWhich}
\alias{vcfWhich<-}

\title{Parameters for scanning VCF files}

\description{
  Use \code{ScanVcfParam()} to create a parameter object influencing 
  which records and fields are imported from a VCF file. Record
  parsing is based on genomic coordinates and requires a Tabix index
  file. Individual VCF elements can be specified in the \sQuote{fixed},
  \sQuote{info} and \sQuote{geno} arguments.
}

\usage{
  ScanVcfParam(fixed=character(), info=character(), geno=character(), 
               trimEmpty=TRUE, which, ...)

  ## Getters and Setters 
 
  vcfFixed(object)
  vcfFixed(object) <- value
  vcfInfo(object)
  vcfInfo(object) <- value
  vcfGeno(object)
  vcfGeno(object) <- value
  vcfTrimEmpty(object)
  vcfTrimEmpty(object) <- value
  vcfWhich(object)
  vcfWhich(object) <- value
}

\arguments{
  \item{fixed}{A character() vector of fixed fields to be returned. Possible
    values are ALT, QUAL and FILTER. The CHROM, POS, ID and REF fields are 
    needed to create the \code{GRanges} of variant locations. Because these 
    are essential fields there is no option to request or omit them. If not 
    specified, all fields are returned; if \code{fixed=NA} only REF is 
    returned. 
  }
  \item{info}{A character() vector naming the \sQuote{INFO} fields to return.
    \code{scanVcfHeader()} returns a vector of available fields. If not 
    specified, all fields are returned; if \code{info=NA} no fields are returned.
  }
  \item{geno}{A character() vector naming the \sQuote{GENO} fields to return.
    \code{scanVcfHeader()} returns a vector of available fields. If not specified, 
    all fields are returned; if \code{geno=NA} no fields are returned.
  }
  \item{trimEmpty}{A logical(1) indicating whether \sQuote{GENO} fields
    with no values should be returned.
  }
  \item{which}{A \code{\linkS4class{GRanges}}, \code{\linkS4class{RangesList}}, 
    \code{\linkS4class{RangedData}} or missing object, describing the sequences 
    and ranges to be queried. Variants whose \code{POS} lies in the interval(s) 
    \code{[start, end]} are returned. If \code{which} is not specified all
    ranges are returned. 
  }
  \item{object}{An instance of class \code{ScanVcfParam}.
  }
  \item{value}{An instance of the corresponding slot, to be assigned to
    \code{object}. 
  }
  \item{\dots}{Arguments passed to methods.
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
    \item{\code{fixed}:}{Object of class \code{"character"} indicating
      fixed fields to be returned.
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
    \item{vcfFixed, vcfInfo, vcfGeno, vcfTrimEmpty, vcfWhich:}{
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

  ## --------------------------------------
  ## 'which' argument
  ## --------------------------------------
  ## To subset on genomic coordinates, create a GRanges, RangedData or 
  ## RangesList with the ranges of interest. This object is supplied
  ## to ScanVcfParam() as the 'which' argument.
  which <- RangesList(seq1=IRanges(1000, 2000), 
      seq2=IRanges(c(100, 1000), c(1000, 2000)))
  ScanVcfParam(which=which)
 
  ## --------------------------------------
  ## 'fixed', 'info' and 'geno' arguments
  ## --------------------------------------
  ## This param specifies the return of the "GT" 'geno' field and the 
  ## subset of ranges in 'which'. All 'fixed' and 'info' fields will be 
  ## returned. 
  ScanVcfParam(geno="GT", which=which)
 
  ## Here two 'fixed' and one 'geno' field are specified
  ScanVcfParam(fixed=c("ALT", "QUAL"), geno="GT", info=NA) 

  ## Return only the 'fixed' fields
  ScanVcfParam(geno=NA, info=NA) 
}

\keyword{classes}
