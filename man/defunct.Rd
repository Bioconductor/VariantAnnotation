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
  ## Defunct in > 1.15.0
  \itemize{
    \item getTranscriptSeqs()
  }

  ## Removed
  \itemize{
    \item refLocsToLocalLocs()
    \item readVcfLongForm()
    \item dbSNPFilter()
    \item regionFilter()
    \item MatrixToSnpMatrix()
  }
}

\details{
  ## Defunct in > 1.15.0
  \itemize{
    \item \code{getTranscriptSeqs} has been replaced by
          \code{extractTranscriptSeqs} in the GenomicFeatures package.
  }

  ## Removed
  \itemize{
    \item{
      \code{refLocsToLocalLocs} has been replaced by \code{mapToTranscripts}
      and \code{pmapToTranscripts}.
    }
    \item{
      \code{readVcfLongForm} has been replaced by \code{expand}.
    }
    \item{
      \code{dbSNPFilter} and \code{regionFilter} have been replaced by 
      \code{filterVcf}.
    }
    \item{
      \code{MatrixToSnpMatrix} has been replaced by \code{genotypeToSnpMatrix}.
    }
  }
}

\author{
  Valerie Obenchain
}

\seealso{
  \itemize{
    \item \code{\link{expand}}
    \item \code{\link{filterVcf}}
    \item \code{\link{genotypeToSnpMatrix}}
    \item \code{\link[GenomicFeatures]{mapToTranscripts}}
    \item \code{\link[GenomicFeatures]{extractTranscriptSeqs}}
  }
}
