\name{writeVcf}
\alias{writeVcf}
\alias{writeVcf,VCF,character-method}
\alias{writeVcf,VCF,connection-method}

\title{Write VCF files}

\description{Write Variant Call Format (VCF) files to disk}

\usage{
\S4method{writeVcf}{VCF,character}(obj, filename, index = FALSE, ...)
\S4method{writeVcf}{VCF,connection}(obj, filename, index = FALSE, ...)
}

\arguments{
  \item{obj}{Object containing data to be written out. At
    present only accepts \linkS4class{VCF}.
  }

  \item{filename}{The character() name of the VCF file, or a connection
    (e.g., \code{\link{file}()}), to be written out. A connection opened
    with \code{open = "a"} will have header information written only if
    the file does not already exist.}

  \item{index}{Whether to bgzip the output file and generate a tabix index.
  }
  \item{\dots}{Additional arguments, passed to methods.
  }
}

\details{A VCF file can be written out from data in a \code{VCF} file. More
  general methods to write out from other objects may be added in the future.

  \code{writeVcf} conforms to the VCF standards on the 1000 Genomes Project
  web site (see references).
}

\value{VCF file
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

\author{Valerie Obenchain <vobencha@fhcrc.org> and Michael Lawrence}

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
  stopifnot(all(in2 == in3))

  ## write incrementally
  out3.vcf <- tempfile()
  con <- file(out3.vcf, open="a")
  writeVcf(in1[1:2,], con)
  writeVcf(in1[-(1:2),], con)
  close(con)
  readVcf(out3.vcf, "hg19")
}

\keyword{manip}
