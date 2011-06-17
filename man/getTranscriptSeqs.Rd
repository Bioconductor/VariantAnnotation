\name{getTranscriptSeqs}
\alias{getTranscriptSeqs}
\alias{getTranscriptSeqs,GRangesList,BSgenome-method}
\alias{getTranscriptSeqs,GRangesList,FaFile-method}
\alias{getTranscriptSeqs,GRanges,FaFile-method}

\title{Get transcript sequences}

\description{
  Get transcript sequences from a \code{\link[BSgenome]{BSgenome}} object or 
  an \link{FaFile}.
}

\usage{
  \S4method{getTranscriptSeqs}{GRangesList,BSgenome}(query, subject, ...)
  \S4method{getTranscriptSeqs}{GRangesList,FaFile}(query, subject, ...)
  \S4method{getTranscriptSeqs}{GRanges,FaFile}(query, subject, ...)
}

\arguments{
  \item{query}{A \link{GRangesList} of coding regions by transcript which
    specifies the sequences to be retrieved from \code{subject}.
  }
  \item{subject}{A \code{\link[BSgenome]{BSgenome}} object or a \link{FaFileList} from which
    the sequences will be taken.
  }
  \item{\dots}{Additional arguments
  }
}

\details{
  Add some ...
}

\value{
  A \link{DNAStringSet} instance containing the sequences for all transcripts
  specified in \code{query}. 
} 

\author{Valerie Obenchain}

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}

\examples{
library(TxDb.Hsapiens.UCSC.hg18.knownGene)
library(BSgenome.Hsapiens.UCSC.hg18)

txdb <- Hsapiens_UCSC_hg18_knownGene_TxDb
cdsByTx <- cdsBy(txdb)
grl <- cdsByTx
grl20 <- keepSeqlevels(grl, "chr20")

## BSgenome as source
bsSource <- getTranscriptSeqs(grl20, Hsapiens)

## fasta files as source
## FIXME : need fasta files
#fafile <- FaFile("/home/vobencha/chr16.fa")
#faSource <- getTranscriptSeqs(grl20, fafile)

# identical(bsSource, faSource)

}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
