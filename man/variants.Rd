\name{variants}
\alias{variants}
\docType{data}

\title{
  Human genetic variant test dataset
}

\description{
  This dataset contains sample Human genetic variants from the ANNOVAR web 
  site and fabricated variants. 
}

\usage{data(variants)}

\format{
  The dataset is an object of class \code{GRanges} with Human 
  variants from build hg18. Some of the varaints come from the ANNOVAR 
  web site and some are fabricated for example purposes. 

  Information for each variant includes chromosome, range, strand,
  reference and observed allele, and a comment.
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
