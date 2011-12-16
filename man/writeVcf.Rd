\name{writeVcf}
\alias{writeVcf}
\alias{writeVcf,character,SummarizedExperiment-method}


\title{Write out a VCF (variant call format) file}

\description{Write out a VCF file - under construction - }

\usage{
  \S4method{writeVcf}{character,SummarizedExperiment}(file, obj, ...)
}

\arguments{
  \item{file}{The character() name of the VCF file to be
    written out.
  }
  \item{obj}{Object containing data to be written out. At
    present only accepts \linkS4class{SummarizedExperiment}.
  }
  \item{\dots}{Additional arguments, passed to methods.
  }
}

\details{
}

\value{
}

\references{
  \url{http://vcftools.sourceforge.net/specs.html} outlines the VCF
  specification.

  \url{http://samtools.sourceforge.net/mpileup.shtml} contains
  information on the portion of the specification implemented by
  \code{bcftools}.

  \url{http://samtools.sourceforge.net/} provides information on
  \code{samtools}.
}

\author{Valerie Obenchain <vobencha@fhcrc.org>}

\seealso{
  \code{\link{readVcf}}
}

\examples{
}

\keyword{manip}
