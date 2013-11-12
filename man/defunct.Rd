\name{VariantAnnotation-defunct}
\alias{VariantAnnotation-defunct}
\alias{readVcfLongForm}
\alias{dbSNPFilter}
\alias{regionFilter}
\alias{MatrixToSnpMatrix}

\title{Defunct Functions in Package \code{VariantAnnotation}}

\description{
  The functions or variables listed here are no longer part
  of \code{VariantAnnotation}.
}

\section{usage}{
  ## Defunct in 1.9.z
  readVcfLongForm()

  ## Removed in 1.9.z
  dbSNPFilter()
  regionFilter()
  MatrixToSnpMatrix(callMatrix, ref, alt, ...)
}

\details{
  ## Defunct in 1.9.z:
  \code{readVcfLongForm} has been replaced by \code{expand}.

  ## Removed in 1.9.z:
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
