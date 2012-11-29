\name{GLtoGP}

\alias{GLtoGP}

\title{Convert genotype likelihoods to genotype probabilities}

\description{
  Convert an array of genotype likelihoods to
  posterior genotype probabilities.
}

\usage{
GLtoGP(gl)
}

\arguments{
  \item{gl}{Array of genotype likelihoods.  The format can be a matrix
  of lists, or a three-dimensional array in which the third dimension
  corresponds to the probabilities for each genotype.}
}

\details{
  Computes the probability of each genotype as \code{10^x / sum(10^x)}.
}

\value{
  An array of posterior genotype probabilities, in the same format as
  the input (matrix of lists or 3D array).
}

\author{
  Stephanie Gogarten  <sdmorris@u.washington.edu>
}

\seealso{
  \link{readVcf},
  \link{genotypeToSnpMatrix}
}

\examples{
  ## Read a vcf file with a "GL" field.
  vcfFile <- system.file("extdata", "gl_chr1.vcf", package="VariantAnnotation") 
  vcf <- readVcf(vcfFile, "hg19")

  ## extract genotype likelihoods as a matrix of lists
  gl <- geno(vcf)$GL
  class(gl)
  mode(gl)

  # convert to posterior probabilities
  gp <- GLtoGP(gl)
}

\keyword{manip}
