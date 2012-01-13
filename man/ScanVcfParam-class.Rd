\name{ScanVcfParam-class}
\Rdversion{1.1}
\docType{class}

\alias{ScanVcfParam}
\alias{ScanVcfParam-class}
\alias{ScanVcfParam,missing-method}
\alias{ScanVcfParam,ANY-method}
\alias{vcfInfo}
\alias{vcfGeno}
\alias{vcfTrimEmpty}
\alias{vcfWhich}
\alias{vcfAsGRanges}

\title{Parameters for scanning VCF files}

\description{
  Use \code{ScanVcfParam()} to create a parameter object influencing the
  \sQuote{INFO} and \sQuote{GENO} fields parsed, and which recrods are
  imported from a VCF file. Use of \code{which} requires that a Tabix 
  index file exists.
}

\usage{
  ScanVcfParam(info=character(), geno=character(), trimEmpty=TRUE,
               which, asGRanges=character(), ...)

  ## Accessors
 
  vcfInfo(object)
  vcfGeno(object)
  vcfTrimEmpty(object)
  vcfWhich(object)
  vcfAsGRanges(object)
}

\arguments{
  \item{info}{A character() vector of \sQuote{INFO} fields (see
    \code{\link{scanVcfHeader}}) to be returned.
  }
  \item{geno}{A character() vector of \sQuote{GENO} fields (see
    \code{\link{scanVcfHeader}}) to be returned. \code{character(0)}
    returns all fields, \code{NA_character_} returns none.
  }
  \item{trimEmpty}{A logical(1) indicating whether \sQuote{GENO} fields
    with no values should be returned.
  }
  \item{which}{An object, for which a method is defined (see usage,
    above), describing the sequences and ranges to be queried. Variants
    whose \code{POS} lies in the interval(s) \code{[start, end]} are
    returned. 
  }
  \item{asGRanges}{A character(1) with the name \sQuote{info} or \sQuote{geno}.
    When specified, a long form \code{GRanges} is returned with the specified 
    field expanded as the elementMetadata columns. For use with \code{readVcf} 
    only; ignored by \code{scanVcf}.
  }
  \item{object}{An instance of class \code{ScanVcfParam}.}

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
    \item{vcfInfo, vcfGeno, vcfTrimEmpty, vcfWhich, vcfAsGRanges:}{Return the
      corresponding field from \code{object}.
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
  p0 <- ScanVcfParam()

  ## subset of reads based on genomic coordinates
  which <- RangesList(seq1=IRanges(1000, 2000),
                      seq2=IRanges(c(100, 1000), c(1000, 2000)))
  p1 <- ScanVcfParam(which=which)
  
  ## return only specified GENO field(s)
  p2 <- ScanVcfParam(geno="GT", which=which)
}

\keyword{classes}
