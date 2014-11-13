\name{VariantAnnotation-defunct}
\alias{VariantAnnotation-defunct}
\alias{refLocsToLocalLocs}
\alias{refLocsToLocalLocs,GRanges,TxDb,missing-method}
\alias{refLocsToLocalLocs,GRanges,missing,GRangesList-method}
\alias{readVcfLongForm}
\alias{readVcfLongForm,TabixFile,character,GRanges-method}
\alias{readVcfLongForm,TabixFile,character,missing-method}
\alias{readVcfLongForm,TabixFile,character,RangedData-method}
\alias{readVcfLongForm,TabixFile,character,RangesList-method}
\alias{readVcfLongForm,TabixFile,character,ScanVcfParam-method}
\alias{readVcfLongForm,character,character,ScanVcfParam-method}
\alias{readVcfLongForm,character,character,missing-method}
\alias{readVcfLongForm,character,missing,missing-method}
\alias{dbSNPFilter}
\alias{regionFilter}
\alias{MatrixToSnpMatrix}

\title{Defunct Functions in Package \code{VariantAnnotation}}

\description{
  The functions or variables listed here are no longer part
  of \code{VariantAnnotation}.
}

\section{usage}{
  ## Defunct in 1.13.z
  \code{refLocsToLocalLocs} has been replaced by \code{mapCoords}.

  ## Removed
  readVcfLongForm()
  dbSNPFilter()
  regionFilter()
  MatrixToSnpMatrix(callMatrix, ref, alt, ...)
}

\details{
  ## Defunct in 1.13.z
  \code{refLocsToLocalLocs} has been replaced by \code{mapCoords}.

  ## Removed 
  \code{readVcfLongForm} has been replaced by \code{expand}.
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
