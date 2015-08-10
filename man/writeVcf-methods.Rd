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
      \itemize{
        \item{nchunk:} Integer or NA. When provided this argument
        overrides the default chunking behavior of \code{writeVcf}, 
        see Details section. An integer value specifies the number
        of records in each chunk; NA disables chunking.
      }
  }
}

\details{A VCF file can be written out from data in a \code{VCF} file. More
  general methods to write out from other objects may be added in the future.

  \code{writeVcf} conforms to the VCF standards on the 1000 Genomes Project
  web site (see references). When 'fileformat' is not present in the header
  data, format is written out according to the standard at the time of the
  Bioconductor release, e.g., for Bioconductor 3.1 the format is VCFv4.2.

  Large VCF files (i.e., > 1e5 records) are written out in chunks; VCF
  files with < 1e5 records are not chunked. The optimal number of records
  per chunk depends on both the number of records and complexity of the
  data. Currently \code{writeVcf} determines records per chunk based on the 
  total number of records only. To override this behavior or experiment with 
  other values use \code{nchunk} as an integer or NA. An integer value 
  represents the number of records per chunk reguardless of the size of the VCF; 
  NA disables all chunking.
  \itemize{
    \item writeVcf(vcf, tempfile())                ## default chunking
    \item writeVcf(vcf, tempfile(), nchunk = 1e6)  ## chunk by 1e6
    \item writeVcf(vcf, tempfile(), nchunk = NA)   ## no chunking
  } 
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
