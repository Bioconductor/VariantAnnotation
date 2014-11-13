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
  \itemize{
    \item refLocsToLocalLocs().
  } 
  ## Removed
  \itemize{
    \item readVcfLongForm()
    \item dbSNPFilter()
    \item regionFilter()
    \item MatrixToSnpMatrix()
  }
}

\details{
  ## Defunct in 1.13.z
  \itemize{
    \item \code{refLocsToLocalLocs} has been replaced by \code{mapCoords}.
  }
  ## Removed
  \itemize{
    \item readVcfLongForm() has been replaced by expand().
    \item dbSNPFilter() and regionFilter() have been replaced by filterVcf().
    \item MatrixToSnpMatrix() has been replaced by genotypeToSnpMatrix().
  }
}

\author{
  Valerie Obenchain <vobencha@fhcrc.org>
}

\seealso{
  \itemize{
    \item \code{\link{expand}}
    \item \code{\link{filterVcf}}
    \item \code{\link{genotypeToSnpMatrix}}
    \item \code{\link{mapCoords}}
  }
}
