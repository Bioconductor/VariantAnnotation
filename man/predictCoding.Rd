\name{predictCoding}
\alias{predictCoding}
\alias{predictCoding,Ranges,TranscriptDb,ANY,character-method}
\alias{predictCoding,GRanges,TranscriptDb,ANY,character-method}
\alias{predictCoding,Ranges,GRangesList,ANY,character-method}
\alias{predictCoding,GRanges,GRangesList,ANY,character-method}

\title{Predict amino acid coding changes for variants}

\description{
  Predict amino acid coding changes for variants that fall in a coding region.
}

\usage{
  \S4method{predictCoding}{GRanges,GRangesList,ANY,character}(query, subject, 
    seqSource, varAllele, ...)
}

\arguments{
  \item{query}{A \link{Ranges} or \link{GRanges} instance 
    containing the variants to be annotated. If a \link[IRanges]{Ranges}
    instance is provided it will be coerced to a \link{GRanges}.
  }
  \item{subject}{A \link{GRangesList} or a \link{TranscriptDb} instance. 

    When \code{subject} is provided as a \link{GRangesList} 
    it is assumed that it was created with the 'cdsBy' function 
    on a \link{TranscriptDb}. If \code{subject} is a \link{TranscriptDb}, the coding
    regions will be identified using the 'cdsBy' function on this object. 
  } 
  \item{seqSource}{A \code{\link[BSgenome]{BSgenome}} instance or an \link{FaFile}
    to be used for sequence extraction.
  }
  \item{varAllele}{A character representing the name of 
    the column in \code{query} that contains the variant alleles. The
    data in the column should be a \link[Biostrings]{DNAStringSet}. Insertions and
    deletions are represented by a missing value.
  }
  \item{\dots}{Additional arguments
  }
}

\details{
  Reference sequences are extracted from fasta files or a 
  \link[BSgenome]{BSgenome} based on the ranges specified in the
  \code{query}. Variant alleles provided in the \code{varAllele} argument
  are substituted into the reference sequences and transcribed. Variant
  sequences are transcribed only if the substitution, insertion or deletion 
  results in a new sequence length divisible by 3. 
} 

\value{
  A \link[IRanges]{DataFrame} of variants that fall within a coding 
  region. Each row represents a variant-transcript match. 
  If a variant matched multiple transcripts, multiple rows are returned for
  the variant. 

  Columns include \code{queryHits}, \code{txID}, \code{refSeq}, \code{varSeq}, 
  \code{refAA}, \code{varAA}, \code{Consequence}, and any metadata that was 
  present in the \code{subject}. \code{queryHits} provides a map to the variants 
  in the original \code{query}. \code{refSeq} and \code{varSeq} contain the reference 
  and variant-modified DNA sequences. Reference and variant amino acid codes are
  provided in \code{refAA} and \code{varAA}. Variant sequences are transcribed only if the 
  substitution, insertion or deletion results in a new sequence length divisible 
  by 3. When a sequence is not transcribed the \code{varAA} column is empty and the 
  variant is classified as a frameshift in the \code{Consequence} column. Possible 
  values for \code{Consequence} are synonymous, nonsynonymous, or frameshift. See 
  the vignette for more details. 
}

\author{Michael Lawrence and Valerie Obenchain}

\seealso{
  \link{locateVariants}
  \link{getTranscriptSeqs}
}

\examples{
  library(TxDb.Hsapiens.UCSC.hg18.knownGene)
  library(BSgenome.Hsapiens.UCSC.hg18)
 
  data(variants)
  txdb <- TxDb.Hsapiens.UCSC.hg18.knownGene 
  aaCoding <- predictCoding(variants, txdb, seqSource=Hsapiens,
      varAllele="varAllele")

  # TODO : example for fasta
}

\keyword{methods}
