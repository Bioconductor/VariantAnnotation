\name{getTranscriptSeqs}
\alias{getTranscriptSeqs}
\alias{getTranscriptSeqs,GRangesList,ANY-method}
\alias{getTranscriptSeqs,GRangesList,BSgenome-method}
\alias{getTranscriptSeqs,GRangesList,FaFile-method}
\alias{getTranscriptSeqs,GRanges,FaFile-method}

\title{Get transcript sequences}

\description{
  This function is defunct. Use GenomicFeatures::extractTranscriptSeqs()
  instead.
 
  Extract transcript sequences from a \link[BSgenome]{BSgenome} object or 
  an \link[Rsamtools]{FaFile}.
}

\usage{
  \S4method{getTranscriptSeqs}{GRangesList,BSgenome}(query, subject, ...)
  \S4method{getTranscriptSeqs}{GRangesList,FaFile}(query, subject, ...)
  \S4method{getTranscriptSeqs}{GRanges,FaFile}(query, subject, ...)
}

\arguments{
  \item{query}{A \link[GenomicRanges]{GRangesList} object containing
  exons or cds grouped by transcript.
  }
  \item{subject}{A \link[BSgenome]{BSgenome} object or a \link[Rsamtools]{FaFile} from which
    the sequences will be taken.
  }
  \item{\dots}{Additional arguments
  }
}

\details{
  \code{getTranscriptSeqs} is a wrapper for the
  \code{extractTranscriptSeqs} and \code{getSeq} functions. The 
  purpose is to allow sequence extraction from either a
  \link[BSgenome]{BSgenome} or \link[Rsamtools]{FaFile}. Transcript
  sequences are extracted based on the boundaries of the feature
  provided in the \code{query} (i.e., either exons or cds regions).
}

\value{
  A \link[Biostrings]{DNAStringSet} instance containing the sequences for all transcripts
  specified in \code{query}. 
} 

\author{Valerie Obenchain <vobencha@fhcrc.org>}

\seealso{
  \link{predictCoding}
  \link[GenomicFeatures]{extractTranscriptSeqs}
  \link[Biostrings]{getSeq}
}

\examples{
## See ?extractTranscriptSeqs in the GenomicFeatures package.
}

\keyword{methods}
\keyword{manip}
