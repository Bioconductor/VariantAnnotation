\name{writeVcf}
\alias{writeVcf}
\alias{writeVcf,VCF,character-method}


\title{Write VCF files}

\description{Write Variant Call Format (VCF) files - under construction - }

\usage{
  \S4method{writeVcf}{VCF,character}(obj, filename, ...)
}

\arguments{
  \item{obj}{Object containing data to be written out. At
    present only accepts \linkS4class{VCF}.
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
  fl <- system.file("extdata", "ex2.vcf", package="VariantAnnotation")
 
  out1.vcf <- tempfile()
  out2.vcf <- tempfile() 
  in1 <- readVcf(fl, "hg19")
  writeVcf(in1, out1.vcf)
  in2 <- readVcf(out1.vcf, "hg19")
  writeVcf(in2, out2.vcf)
  in3 <- readVcf(out2.vcf, "hg19")

  identical(in2, in3)
}

\keyword{manip}
