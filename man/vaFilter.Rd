\name{vaFilter}
\alias{vaFilter}
\alias{vaFilter,missing-method}
\alias{vaFilter,function-method}
\alias{dbSNPFilter}
\alias{regionFilter}
\alias{compose}

\title{Functions for VariantAnnotation filters}

\description{

  These functions create instances of \code{\linkS4class{VAFilter}}  
  objects. Filters can be applied to genomic variants in a 
  \code{\linkS4class{GRanges}} object. The object returned contains a logical 
  vector indicating which records satisfied the filter and stats on how many 
  elements satisfied the filter.

}

\usage{
vaFilter(fun, name = NA_character_, ...)
\S4method{vaFilter}{missing}(fun, name=NA_character_, ...)
\S4method{vaFilter}{function}(fun, name=NA_character_, ...)

compose(filt, ..., .name)

dbSNPFilter(dbSNP=character(0), .name="dbSNPFilter")
regionFilter(txdb, region="coding", .name="regionFilter")
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
  \item{txdb}{The \code{\linkS4class{TranscriptDb}} object used
    to identifying gene regions.
  }
  \item{region}{A \code{character(1)} object specifying the region
    on which to filter the results. Possible values are "coding",
    "intron", "3UTR", "5UTR" and "intergenic".
  }
  \item{\dots}{Additional arguments
  }
}

\details{

  \code{vaFilter} allows users to construct their own filters. The
  \code{fun} argument to \code{vaFilter} must be a function accepting
  a single argument \code{x} and returning a logical vector indicating
  which records passed the filter. 

  The \code{signature(fun="missing")} method creates a default filter
  that returns a vector of \code{TRUE} values with length equal to
  \code{length(x)}.

  \code{compose} constructs a new filter from one or more existing
  filter. The result is a filter that returns a logical vector with
  indices corresponding to components of \code{x} that pass all
  filters. If not provided, the name of the filter consists of the names
  of all component filters, each separated by \code{" o "}.

  The remaining functions documented on this page are built-in filters
  that accept arguments \code{x, subset=TRUE} and return a logical vector of
  \code{length(x)} indicating which components of \code{x} satisfy the
  filter.

  \code{dbSNPFilter} selects elements present in the specified dbSNP
   package 

  \code{regionFilter} selects elements present in the region specified.
}

\value{

  \code{vaFilter} returns an object of \code{\linkS4class{VAFilter}}.

  Filters return a logical vector of \code{length(x)}, with
  \code{TRUE} indicating components that pass the filter.

}

\author{Valerie Obenchain <vobencha@fhcrc.org>}

\seealso{\code{\linkS4class{VAFilter}}.}

\examples{
data(variants)

# variants present in dbSNP 
snpFilt <- dbSNPFilter("SNPlocs.Hsapiens.dbSNP.20090506")
# subset of variants that passed the filter
snpFilt(variants)
# VAFilterResult with logical vector indicating which variants passed the filter
snpFilt(variants, subset=FALSE)

# variants present in coding regions 
library(TxDb.Hsapiens.UCSC.hg18.knownGene)
txdb <- Hsapiens_UCSC_hg18_knownGene_TxDb
regionFilt <- regionFilter(txdb, region="coding")
regionFilt(variants, subset=FALSE)

# subset of variants that passed both filters
comboFilt <- compose(snpFilt, regionFilt)
comboFilt(variants)

}


\keyword{manip}
