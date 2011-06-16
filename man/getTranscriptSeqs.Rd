\name{getTranscriptSeqs}
\alias{getTranscriptSeqs,GRangesList,BSgenome-method}
\alias{getTranscriptSeqs,GRangesList,FaFile-method}
\alias{getTranscriptSeqs,GRanges,FaFile-method}
\alias{getTranscriptSeqs}

\title{Get transcript subject}

\description{
  Get transcript subject from a \link{BSgenome} object or 
  an \link{FaFile}.
}

\usage{
  \S4method{getTranscriptSeqs}{GRangesList, BSgenome}(query, subject, ...)
  \S4method{getTranscriptSeqs}{GRangesList, FaFile}(query, subject, ...)
  \S4method{getTranscriptSeqs}{GRanges, FaFile}(query, subject, ...)
}

\arguments{
  \item{query}{A \link{GRangesList} of coding regions by transcript which
    specifies the sequences to be retrieved from \code{seqSource}.
  }
  \item{subject}{A \link{BSgenome} object or a \link{FaFileList} from which
    the sequences will be taken.
  }
  \item{\dots}{Additional arguments
  }
}

\details{
  Add some ...
}

\value{
  A \code{\link{GRanges}} object containing only the variants from
  \code{query} that overlap with a coding sequence in a transcript. A row is
  returned for each transcript the variant overlapped. Three metadata columns
  are added, 'txID', 'aaRef', and 'aa'. The 'txID' corresponds to the
  transcript id found in the cdsByTx object. 'refAA' and 'obsAA' represent
}

\author{Valerie Obenchain}


\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}

\examples{
library(TxDb.Hsapiens.UCSC.hg18.knownGene)
library(BSgenome.Hsapiens.UCSC.hg18)

## create or find fasta files
#fafile <- FaFile("/home/vobencha/chr16.fa")

txdb <- Hsapiens_UCSC_hg18_knownGene_TxDb
cdsByTx <- cdsBy(txdb)
grl <- cdsByTx
grl20 <- keepSeqlevels(subject, "chr20")

## BSgenome as source
bsSource <- getTranscriptSeqs(glr20, Hsapiens)

## fasta files as source
#faSource <- getTranscriptSeqs(grl20, fafile)

# identical(bsSource, faSource)

}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
