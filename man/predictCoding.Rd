\name{predictCoding}
\alias{predictCoding}
\alias{predictCoding,Ranges,TranscriptDb-method}
\alias{predictCoding,Ranges,character-method}
\alias{predictCoding,GRanges,TranscriptDb-method}
\alias{predictCoding,GRanges,character-method}
\alias{predictCoding,Ranges,GRangesList-method}
\alias{predictCoding,GRanges,GRangesList-method}

\title{Predict amino acid coding changes for variants}

\description{
  Predict amino acid coding changes for variants that fall in a coding region.
}

\usage{
  \S4method{predictCoding}{GRanges,GRangesList}(query, subject, seqSource,
    refAllele, varAllele, ...)
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
  \item{refAllele}{A character representing the name of the column 
    in \code{query} that contains the reference alleles. The
    data in the column should be a DNAStringSet. Insertions and
    deletions are represented by a missing value.
  }
  \item{varAllele}{A character representing the name of the column 
    in \code{query} that contains the variant alleles. The
    data in the column should be a DNAStringSet. Insertions and
    deletions are represented by a missing value.
  }
  \item{\dots}{Additional arguments
  }
}

\details{
  Variants that fall in a coding region are identified and categorized as 
  synonymous, nonsynonymous or a frameshift. Amino acid coding changes 
  are provided for substitutions, insertions and deletions that change
  by a nucleotide multiple of 3.

  Sequences for the variants can be extracted from fasta files or a 
  \code{\link[BSgenome]{BSgenome}}. See the vignette for further details.
} 

\value{
  A \code{\link{DataFrame}} with columns of queryIndex, txID, refCodon,
  varCodon, refAA, varAA, Consequence, and any metadata columns that
  were present in the \code{subject}. 

  The \code{link{DataFrame}} contains data on the variants from
  \code{query} that overlap with a coding sequence in a transcript
  in the \code{subject}. One row is returned for each transcript the variant 
  overlapped. If the variant overlaps multiple transcripts, multiple
  rows are returned for that variant. The variant from \code{query} is
  identified via the 'queryIndex' column. Not all variant codons are
  transcribed. If the substitution, insertion or deletion results in
  a basepair sequence not divisible by 3, the codon(s) are not transcribed.
  In this case, the varAA column will be empty and the variant will be
  classified as a framshift in the Consequence column. Possible values
  for the Consequence column are synonymous, nonsynonymous,
  or frameshift. 
}

\author{Michael Lawrence and Valerie Obenchain}

\seealso{
  \code{\link{locateVariants}}
}

\examples{
  \dontrun{
  library(TxDb.Hsapiens.UCSC.hg18.knownGene)
  library(BSgenome.Hsapiens.UCSC.hg18)
 
  data(variants)
  txdb <- Hsapiens_UCSC_hg18_knownGene_TxDb 
  aaCoding <- predictCoding(variants, txdb, seqSource=Hsapiens,
      refAllele="refAllele", varAllele="varAllele")

  ## FIXME : example for fasta
  }
}

\keyword{}
