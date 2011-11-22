\name{VAFilter-class}
\docType{class}

% Class:
\alias{VAFilter-class}
\alias{class:VAFilter}
\alias{VAFilter}

% Constructors:
\alias{vaFilter,VAFilter-method}
\alias{vaFilter,missing-method}
\alias{vaFilter,function-method}
\alias{vaFilter}
\alias{dbSNPFilter}
\alias{regionFilter}
\alias{compose}

% Other methods:
\alias{name,VAFilter-method}
\alias{name}
\alias{show,VAFilter-method}


\title{VAFilter objects}

\description{
  Objects of the \code{VAFilter} class are functions for filtering variants.
  A \code{VAFilter} instance applied to a \code{\linkS4class{GRanges}}
  returns either a subset of the records that passed the filter or a 
  \code{\linkS4class{VAFilterResult}} object. The \code{VAFilterResult} is a
  logical vector indicating which records passed the filter and contains
  summary statistics on the records that passed.
}

\details{
  A \code{VAFilter} is constructed by using one of the built-in filters, 
  \code{dbSNPFilter} or \code{regionFilter}, or a user can create their 
  own with \code{vaFilter}. Once the filter is created it can be applied to a
  \code{\linkS4class{GRanges}} object. The format of the output is dictated by 
  the optional\code{subset} argument. If \code{subset=TRUE} (the default) a 
  \code{\linkS4class{GRanges}} is returned containing only the records that 
  passed the filter. If \code{subset=FALSE} a \code{\linkS4class{VAFilterResult}} 
  object is returned. This object contains a logical vector indicating which 
  records passed the filter and summary statistics on the records that passed.
}

\section{Slots}{
  \describe{
    \item{\code{.Data}:}{Object of class \code{"function"} taking a
      single named argument \code{x} corresponding to the 
      \code{\link{GRanges}} object that the filter will be applied to. 
      The return value of the filter function is expected to be a logical 
      vector that can be used to subset \code{x} to include those elements 
      of \code{x} satisfying the filter.
    }
    \item{\code{name}:}{Object of class \code{"ScalarCharacter"}
      representing the name of the filter. The name is useful for
      suggesting the purpose of the filter, and for debugging failed
      filters.
    }
  }
}

\section{Constructors}{
  Built-in filters :
  \describe{
    \item{}{
      \code{dbSNPFilter(dbSNP=character(0), name="dbSNPFilter")}:
      Overlaps the records with a range width of 1 with snp locations in the 
      specified SNPlocs package.
    }
    \item{}{
      \code{regionFilter(txdb, region="coding", name="regionFilter")}:
      Records will be filtered by the name of the region supplied. Possible
      regions are the same as those in the \code{Location} output from
      \code{locateVariants}: `coding', `intronic', `3UTR', `5UTR', `intergenic',
      `noGeneMatch' and `unknown'. This regions is specified when calling the
      filter on a \code{\linkS4class{GRanges}} object. 
    }
  }
  Combining filters :
  \describe{
    \item{}{
      \code{compose(filt, ..., .name)}:
      The \code{compose} function constructs a new filter from one or more 
      existing filters. The result is a filter that returns records that pass 
      the criteria for all filters. If a name is not provided for the filter, a
      name will be constructed by joining the names of all component filters
      separated by \code{" o "}.
    }
  }
  User defined filters :
  A user can construct their own filter using the \code{vaFilter} function. The
  \code{fun} argument must be a function accepting a single argument \code{x}
  and returning a logical vector indicating which records passed the filter.
  \describe{
    \item{}{
      \code{vaFilter(fun="VAFilter", name=NA_character_, ...)}:
      Returns the function representing the underlying filter. This is primarily 
      for viewing the filter function.
    }
    \item{}{
      \code{vaFilter(fun="missing", name=NA_character_, ...)}:
      Creates a default filter that returns a vector of \code{TRUE} values with 
      length equal to \code{length(x)}.
    }
  }
}

\section{Other Methods}{
  \describe{
    \item{name}{\code{signature(x = "VAFilter")}: Return, as a
      \code{ScalarCharacter}, the name of the function.
    }
    \item{show}{\code{signature(object = "VAFilter")}: display a brief
      summary of the filter
    }
  }
}

\arguments{
  \item{fun}{An object of class \code{function} to be used as a
    filter. \code{fun} must accept a single named argument \code{x}, and
    is expected to return a logical vector such that \code{x[fun(x)]}
    selects only those elements of \code{x} satisfying the conditions of
    \code{fun}
  }
  \item{name}{A \code{character(1)} object to be used as the name of the
    filter.  The \code{name} is useful for debugging and reference.
  }
  \item{filt}{A \code{\linkS4class{VAFilter}} object, to be used with
    additional arguments to create a composite filter.
  }
  \item{.name}{An optional \code{character(1)} object used to over-ride
    the name applied to default filters.
  }
  \item{dbSNP}{The \code{character(0)} name of the dbSNP package
    to be used. The default (\code{character(0)})
    performs no filtering
  }
  \item{txdb}{The \link[GenomicFeatures]{TranscriptDb} object used
    to identifying gene regions.
  }
  \item{region}{A \code{character(1)} object specifying the region
    on which to filter the results. Possible values are `coding',
    `intron', `3'UTR', `5'UTR', `intergenic', `noGeneMatch' and
    `unknown'. See ?\code{locateVariants} for regions descriptions.
  }
  \item{\dots}{Additional arguments to methods.
  }
}

\author{Valerie Obenchain <vobencha@fhcrc.org>}

\seealso{
  \code{\link{VAFilterResult}}
}

\examples{
  data(variants)

  ## Create a filter that will identify variants present in dbSNP : 
  snpFilt <- dbSNPFilter("SNPlocs.Hsapiens.dbSNP.20110815")
  ## Apply the filter and control the output :
  ## Optional subset argument is TRUE (default) 
  snpFilt(variants)
  ## Optional subset argument is FALSE 
  snpFilt(variants, subset=FALSE)
 
  ## Create a filter to select only records in 3'UTR : 
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  txdb19 <- TxDb.Hsapiens.UCSC.hg19.knownGene
  regionFilt <- regionFilter(txdb19, region="coding")
  regionFilt(variants, subset=FALSE)
}

\keyword{classes}
