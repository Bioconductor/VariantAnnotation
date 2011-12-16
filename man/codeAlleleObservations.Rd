\name{codeAlleleObservations}
\alias{codeAlleleObservations}

\title{Convert genotypes from nucleotide to integer representation}

\usage{codeAlleleObservations(ref, alt, observed, ...)}

\description{Convert genotypes from nucleotide to integer representation}

\arguments{
  \item{ref}{A \code{character} vector of reference alleles. Each reference
    can be one or more nucleotides long. \code{ref} must be the same length
    as \code{alt} and \code{observed}. Missings are represented as "-".
  }
  \item{alt}{A \code{list} of alternate (variant) alleles. Each alternate can 
    be one or more nucleotides long. \code{alt} need not include the reference 
    alleles specified in \code{ref}. This \code{list} should include all 
    alternative alleles that are possible for each variant.
  }
  \item{observed}{A \code{character} vector of alleles reflecting the genotypes
    for a given sample. One genotype per sample.
  }
  \item{\dots}{Additional arguments to methods.
  }
}

\details{
  
  Take references bases (vector of A, C, G, or T) and alleles (list of
  vectors with allowed letters (i.e list( c("A","C") )). Return column
  for observed alleles in one sample.

}

\value{
  A \code{character} vector of integer genotypes where 0 represents the 
  reference allele and integers > 0 represent the alternate. All 
  genotypes are considered unphased and separated with a "/". 
}

\author{Peter M. Haverty <phaverty@gene.com>}

\seealso{
}

\examples{
  refbase <- c("A","C","A","G")
  altbase <- list("T", c("A","G"), c("A"), c("C"))
  observed <- list(c("A","A"), c("A","C"), c("A", "A"), c("G","C"))
  codeAlleleObservations(refbase, altbase, observed)

  ## Missings are specified with "-"
  refbase <- c("G","G")
  altbase <- list("C", c("A","C", "TCAGT"))
  observed <- list(c("-", "G"), c("A","TCAGT","C"))
  codeAlleleObservations(refbase, altbase, observed)
}

\keyword{methods}
