\name{variants}
\alias{variants}
\docType{data}

\title{
  A genetic variant test dataset
}

\description{
  The genetic varaints used in the ANNOVAR examples
}

\usage{data(variants)}

\format{
  The dataset is an object of class \code{GRanges} with 12 variants
  from human. These data are from the ANNOVAR site and used to
  demonstrant the annotation of multiple types of variants.

  Information for each variant includes chromosome, range, strand,
  reference and observed allele, and a comment describing the varaint.
}

\source{
  ANNOVAR web site :
  \url{http://www.openbioinformatics.org/annovar/}
}
\references{
  Wang K., Li M, and Hakonarson H. (2010), \dQuote{ANNOVAR: functional
  annotation of genetic variants from high-throughput sequencing data}.
  \emph{Nucleic Acids Research}, Vol. 38, No. 16, e164.
}

\examples{
data(variants)
variants
}
\keyword{datasets}
