\name{probabilityToSnpMatrix}

\alias{probabilityToSnpMatrix}

\title{Convert posterior genotype probability to a SnpMatrix object}

\description{
  Convert a matrix of posterior genotype probabilites P(AA), P(AB),
  P(BB) to a \link[snpStats:SnpMatrix-class]{SnpMatrix}.
}

\usage{
  probabilityToSnpMatrix(probs)
}

\arguments{
  \item{probs}{Matrix with three columns for the posterior probabilities
  of the three genotypes: "P(A/A)", "P(A/B)", "P(B/B)".  Each row must
  sum to 1.}
}

\details{
  \code{probabilityToSnpMatrix} converts a matrix of posterior
  probabilites of genotype calls into a 
  \link[snpStats:SnpMatrix-class]{SnpMatrix}.
}

\value{
  An object of class \code{"SnpMatrix"} with one row (one sample).
  Posterior probabilities are encoded (approximately) as byte
  values, one per SNP. See 
    the help page for \link[snpStats:SnpMatrix-class]{SnpMatrix} for complete 
    details of the class structure.
}

\author{
  Stephanie Gogarten <sdmorris@u.washington.edu>
}

\seealso{
  \link{genotypeToSnpMatrix},
  \link[snpStats:SnpMatrix-class]{SnpMatrix}
}

\examples{
probs <- matrix(c(1,0,0,
                  0,1,0,
                  0,0,1,
                  NA,NA,NA),
                  ncol=3, byrow=TRUE,
                  dimnames=list(1:4,c("A/A","A/B","B/B")))
sm <- probabilityToSnpMatrix(probs)
as(sm, "character")
}

\keyword{manip}
