\name{globalToLocal}
\alias{globalToLocal}

\title{globalToLocal}

\description{
  Converts reference-based (aka chromosome-based or genomic)
  locations into transcript-based locations
}

\usage{globalToLocal(global, ranges)}

\arguments{
  \item{global}{A \code{\linkS4class{GRanges}} object}
  \item{ranges}{A \code{\linkS4class{GRangesList}} object}
}

\details{
  This function translates coordinates from a reference location specified
  in \code{global} to the transcript location as defined in \code{ranges}.
}

\value{
  A \link[IRanges]{DataFrame} with three columns of 'global.ind', ranges.ind'
  and 'local'. 'global.ind' and 'ranges.ind' are the indices
  of the \code{global} and \code{ranges} input arguments, respectively. The
  'local' column contains the range coordinates of the elements in
  \code{global} translated with respect to the ranges in 'ranges.ind'. 
}


\author{Michael Lawrence}


\seealso{
\code{getTranscriptSeqs}
\link[GenomicFeatures]{transcriptLocs2refLocs}
\link[GenomicFeatures]{extractTranscriptsFromGenome}
}

\examples{
## under construction

#library(BSgenome.Hsapiens.UCSC.hg19)
#library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#
### snps from dbsnp132
#snvs <- GRanges(seqnames=c("chr6", "chr10", "chr16"), 
#ranges = IRanges(start=c(88375601, 115912481, 69875501), width=c(1,2,2)))
#
#txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene 
#genes <- transcriptsBy(txdb, "gene")
#rgs <- reduce(genes[subjectHits(fo), ])
#fo <- findOverlaps(snps, rgs)
#txseqs <- extractTranscriptsFromGenome(Hsapiens, rgs)
#
#txLocal <- globalToLocal(snvs, rgs[subjectHits(fo), ])
#
### retrieve reference sequences 
### reference coordinates
#rlocs <- transcriptLocs2refLocs(tlocs, start(exbytx), end(exbytx),
#             tx_strand, reorder.exons.on.minus.strand=TRUE)
#
# subject <- subject[unique(subjectHits(fo))]
#        txSeqs <- getTranscriptSeqs(subject, seqSource)
#        xCoding <- query[txLocal$globalInd]
#
# back to reference coordinates

}

\keyword{manip}
