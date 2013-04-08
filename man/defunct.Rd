\name{VariantAnnotation-defunct}
\alias{VariantAnnotation-defunct}
\alias{dbSNPFilter}
\alias{regionFilter}
\alias{MatrixToSnpMatrix}
\alias{MatrixToSnpMatrix,matrix,DNAStringSet,DNAStringSetList-method}

\title{Defunct Functions in Package \code{VariantAnnotation}}

\description{
  The functions or variables listed here are no longer part
  of \code{VariantAnnotation}.
}

\section{usage} {
  ## Defunct in 1.7.z
  dbSNPFilter()
  regionFilter()
  MatrixToSnpMatrix(callMatrix, ref, alt, ...)
}

\details{
  \code{dbSNPFilter}, \code{regionFilter} and \code{MatrixToSnpMatrix}
  have stubs which report they are defunct.

  \code{dbSNPFilter} and \code{regionFilter} have been replaced by
  \code{filterVcf}.

  \code{MatrixToSnpMatrix} has been replaced by \code{genotypeToSnpMatrix}.
}

\author{
  Valerie Obenchain <vobencha@fhcrc.org>
}

\seealso{
  \code{\link{filterVcf}}, \code{\link{genotypeToSnpMatrix}}.
}
