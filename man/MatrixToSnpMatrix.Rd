\name{MatrixToSnpMatrix}
\alias{MatrixToSnpMatrix}


\title{Convert a matrix of genotype calls to a SnpMatrix object}

\description{
  Convert a matrix of genotype calls of format "x/x" to a
  \linkS4class{SnpMatrix} object.
}

\usage{MatrixToSnpMatrix(callMatrix, ...)}

\arguments{
  \item{callMatrix}{A \code{matrix} of genotype calls to be
    converted to a \linkS4class{SnpMatrix}. The columns of
    the \code{matrix} should be samples and the rows
    the snps. The format of the calls is expected to be 
    \code{character} of the form "x/x" where "x" is an integer 
    value representing the major (reference) or minor (other)
    allele. Typically values will be 0 = major, 1 = minor.
  }
  \item{\dots}{Additional arguments, passed to methods.
  }
}

\details{
  \code{MatrixToSnpMatrix} converts genotype calls of the form
  "x/x" to a \linkS4class{SnpMatrix} object. Values of "x" are
  expected to be integers >= 0. The calls are converted
  to values consistant with \linkS4class{SnpMatrix} where
  0 = missing, 1 = homozygous major, 2 = heterozygous and 
  3 = homozygous minor.
}

\value{
  The object returned is a \linkS4class{SnpMatrix} where the
  columns are snps and the rows are samples. 
  See the help page for \linkS4class{SnpMatrix} for
  complete details of the class structure.
}

\author{Valerie Obenchain <vobencha@fhcrc.org>}

\seealso{
  \linkS4class{SnpMatrix}
}

\examples{

## import vcf file to SummarizedExperiment
vcfFile <- system.file("extdata", "ex1.vcf", package="VariantAnnotation") 
se <- readVcf(vcfFile, index=character(0))

## genotype calls are stored in the assays slot
assays(se)
calls <- assays(se)$GT
snps <- MatrixToSnpMatrix(calls)
}


\keyword{manip}
