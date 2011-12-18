\name{writeVcf}
\alias{writeVcf}
\alias{writeVcf,SummarizedExperiment,character-method}


\title{Write out a VCF (variant call format) file}

\description{Write out a VCF file - under construction - }

\usage{
  \S4method{writeVcf}{SummarizedExperiment,character}(obj, filename, ...)
}

\arguments{
  \item{obj}{Object containing data to be written out. At
    present only accepts \linkS4class{SummarizedExperiment}.
  }
  \item{filename}{The character() name of the VCF file to be
    written out.
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
#  fl <- system.file("extdata", "ex2.vcf", package="VariantAnnotation")
#  se <- readVcf(fl, "hg19")
#  writeVcf(se, "out.vcf")
}

\keyword{manip}
