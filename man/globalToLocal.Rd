\name{globalToLocal}
\alias{globalToLocal}

\title{globalToLocal}

\description{
  Converts reference-based (aka chromosome-based or genomic)
  locations into transcript-based locations
}

\usage{globalToLocal(query, subject)}

\arguments{
  \item{query}{A \code{\linkS4class{GRanges}} object}
  \item{subject}{A \code{\linkS4class{GRangesList}} object}
}

\details{
  This function translates coordinates from a reference location specified
  in \code{query} to the transcript location as defined in \code{subject}.

  The \code{subject} must be an object created by the \code{cdsBy} function.
  This should be a method which enforces this ... TBD.
}

\value{
  A \link[IRanges]{DataFrame} with columns of \sQuote{qindex}, \sQuote{sindex},
  \sQuote{cdsid} and \sQuote{local}. The \sQuote{local} column contains the 
  range coordinates of the \code{query} translated with respect to the ranges 
  in the \code{subject}. 
}

\author{Michael Lawrence and Valerie Obenchain <vobencha@fhcrc.org>}

\seealso{
  \code{getTranscriptSeqs}
  \link[GenomicFeatures]{transcriptLocs2refLocs}
  \link[GenomicFeatures]{extractTranscriptsFromGenome}
}

\examples{
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene 
cdsbytx <- cdsBy(txdb)
query <- GRanges("chr22", IRanges(c(1000:2000), width=150))

globalToLocal(query, cdsbytx)
}

\keyword{manip}
