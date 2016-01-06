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
\alias{VRangesScanVcfParam}
\alias{restrictToSNV}

\title{Defunct Functions in Package \code{VariantAnnotation}}

\description{
  The functions or variables listed here are no longer part
  of \code{VariantAnnotation}.
}

\section{usage}{
  ## Defunct in > 1.17.0
  \itemize{
    \item VRangesScanVcfParam()
    \item restrictToSNV()
  }

  ## Removed
  \itemize{
    \item refLocsToLocalLocs()
    \item readVcfLongForm()
    \item dbSNPFilter()
    \item regionFilter()
    \item MatrixToSnpMatrix()
    \item getTranscriptSeqs()
  }
}

\details{
  ## Defunct in > 1.17.0
  \itemize{
    \item{
      \code{VRangesScanVcfParam} has been replaced by \code{ScanVcfParam}.
    }
    \item{
      \code{restrictToSVN} has been replaced by \code{isSNV}.
    }
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
      \code{regionFilter} has been replaced by \code{filterVcf}.
    }
    \item{
      \code{MatrixToSnpMatrix} has been replaced by \code{genotypeToSnpMatrix}.
    }
    \item \code{getTranscriptSeqs} has been replaced by
          \code{extractTranscriptSeqs} in the GenomicFeatures package.
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
    \item \code{\link{isSNV}}
    \item \code{\link{ScanVcfParam}}
  }
}
