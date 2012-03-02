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
  This function translates reference-based coordinates defined in 
  the \code{query} to transcript-based coordinates.

  The \code{subject} must be an object created by the \code{cdsBy} function.
  This should be a method which enforces this ... TBD.
}

\value{
  A \link[IRanges]{DataFrame} with the following columns of \sQuote{qindex}, \sQuote{sindex},
  \sQuote{cdsid}, \sQuote{strand}, \sQuote{txloc}, and \sQuote{cdsloc}. The 
  \sQuote{txloc} column contains the transcript-based range coordinates of the \code{query}.
  \describe{
    \item{\code{qindex}}{
      Row number of original \code{query}.
    }
    \item{\code{sindex}}{
      Row number of original \code{subject}.
    }
    \item{\code{cdsid}}{
      cds_id taken from the \code{subject}.
    }
    \item{\code{strand}}{
      Strand of the \code{subject} the \code{query} overlapped.
    }
    \item{\code{txloc}}{
      Location of \code{query} in transcript-based coordinates.
    }
    \item{\code{cdsloc}}{
      Location of \code{query} in cds-based coordinates.
    }
  }
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
query <- GRanges("chr22", IRanges(c(16279195, 51208333), width=10))

globalToLocal(query, cdsbytx)
}

\keyword{manip}
