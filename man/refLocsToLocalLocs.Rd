\name{refLocsToLocalLocs}
\alias{refLocsToLocalLocs}
\alias{refLocsToLocalLocs,GRanges,TranscriptDb,missing,missing-method}
\alias{refLocsToLocalLocs,GRanges,missing,GRangesList,GRangesList-method}

\title{refLocsToLocalLocs}

\description{
  Converts reference locations (aka chromosome-based or genomic)
  to coding regions, and protein based locations
}

\usage{
refLocsToLocalLocs(ranges, txdb, cdsbytx, exbytx, ...)
\S4method{refLocsToLocalLocs}{GRanges,TranscriptDb,missing,missing}(ranges, txdb, cdsbytx, exbytx, ...)
\S4method{refLocsToLocalLocs}{GRanges,missing,GRangesList,GRangesList}(ranges, txdb, cdsbytx, exbytx, ...)
}

\arguments{
  \item{ranges}{A \link[GenomicRanges]{GRanges}es]{GRanges} object containing 
    the variants in reference-based coordinates. 
  }
  \item{txdb}{A \link[GenomicFeatures]{TranscriptDb} object that serves
    as the annotation reference.
  }
  \item{cdsbytx}{A \link[GenomicRanges]{GRangesList} object with transcripts
    as the outer list elements and coding regions as the inner.
  }
  \item{exbytx}{A \link[GenomicRanges]{GRangesList} object with transcripts
    as the outer list elements and exons as the inner.
  }
  \item{\dots}{Additional arguments passed to methods
  }
}

\details{
  This function translates the reference-based coordinates in \code{ranges} to 
  \sQuote{local} coordinates in the coding region (CDS) and protein sequences.

  When a \code{txdb} is suplied the \code{cdsbytx} and \code{exbytx} are created 
  from calls to the \code{cdsBy()} and \code{exonsBy()} functions. When 
  \code{cdsbytx} and \code{exbytx} are provided the outer list elements must be 
  transcripts and the inner list elements represent coding regions and exons, 
  respectively. The \code{GRangesList} objects must have names on the outer list 
  elements (i.e., transcript names).

  Only ranges that fall \sQuote{within} coding sequences are reported in the 
  result. Output is a modified \code{GRanges} of the \code{ranges} input where 
  each row represents a range-transcript match making multiple rows per range 
  posible. The \code{elementMetadata} columns include \code{tx_id}, 
  \code{cdsLoc} and \code{proteinLoc}. When a \code{txdb} is provided the 
  \code{txId} is the internal transcript id from the annotation. When 
  \code{cdsbytx} and \code{exbytx} are provided \code{tx_id} are the names 
  on the outer list elements of the \code{GRangesList} objects.
}

\value{
  A \link[GenomicRanges]{GRanges} with the following \code{elementMetadata}
  columns,

  \describe{
    \item{\code{cdsLoc}}{
      Location in coding region coordinates
    }
    \item{\code{proteinLoc}}{
      Location in protein (codon triplet) coordinates
    }
    \item{\code{queryID}}{
      Character vector mapping to the of rows of the original \code{query}
    }
    \item{\code{txID}}{
      Character vector of internal transcript ids from the 
      \link[GenomicFeatures]{TranscriptDb} or the names of the outer list 
      elements of the \code{cdsbytx} and \code{exbytx} objects.
    }
    \item{\code{cdsID}}{
      Character vector of internal coding region ids from the 
      \link[GenomicFeatures]{TranscriptDb} or the names of the outer list 
      elements of the \code{cdsbytx} and \code{exbytx} objects.
    }
  }
}

\author{Michael Lawrence and Valerie Obenchain <vobencha@fhcrc.org>}

\seealso{
  \link[IRanges]{map}
  \code{predictCoding}
  \code{getTranscriptSeqs}
  \link[GenomicFeatures]{transcriptLocs2refLocs}
  \link[GenomicFeatures]{extractTranscriptsFromGenome}
}

\examples{
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
ranges <- GRanges(c("chr12", "chr1", "chr5"), 
                 IRanges(c(1017956, 881906, 140532), 
                         c(1017956, 881907, 140532)))
refLocsToLocalLocs(ranges, TxDb.Hsapiens.UCSC.hg19.knownGene)
}

\keyword{manip}
