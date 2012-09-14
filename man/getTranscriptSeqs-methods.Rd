\name{getTranscriptSeqs}
\alias{getTranscriptSeqs}
\alias{getTranscriptSeqs,GRangesList,BSgenome-method}
\alias{getTranscriptSeqs,GRangesList,FaFile-method}
\alias{getTranscriptSeqs,GRanges,FaFile-method}

\title{Get transcript sequences}

\description{
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
  \code{extractTranscriptsFromGenome} and \code{getSeq} functions. The 
  purpose is to allow sequence extraction from either a
  \link[BSgenome]{BSgenome} or \link[Rsamtools]{FaFile}. Transcript
  sequences are extracted based on the boundaries of the feature
  provided in the \code{query} (i.e., either exons or cds regions).
}

\value{
  A \link[Biostrings]{DNAStringSet} instance containing the sequences for all transcripts
  specified in \code{query}. 
} 

\author{Valerie Obenchain}

\seealso{
  \link{predictCoding}
  \link[GenomicFeatures]{extractTranscriptsFromGenome}
  \link[Biostrings]{getSeq}
}

\examples{
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  library(BSgenome.Hsapiens.UCSC.hg19)
 
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene 
  cdsByTx <- cdsBy(txdb)
  chr20 <- keepSeqlevels(cdsByTx, "chr20")

  ## BSgenome as sequence source
  bsSource <- getTranscriptSeqs(chr20, Hsapiens)
}

\keyword{methods}
\keyword{manip}
