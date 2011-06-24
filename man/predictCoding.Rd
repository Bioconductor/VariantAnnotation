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
  \S4method{predictCoding}{GRanges,GRangesList}(query, subject, alleles, 
    seqSource, ...) 
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
  \item{alleles}{A character representing the name of the column 
    in \code{query} that contains the reference and variant alleles. The
    data in the column must be character and have format 'reference/variant'.
    For example, 'A/G' would indicate the reference allele is A and the 
    variant is G. An two nucleotide insertion would be represented by
    '-/AT', and a four nucleotide deletion would be 'GCCT/-'. 
  }
  \item{seqSource}{A \code{\link[BSgenome]{BSgenome}} instance or an \link{FaFile}
    to be used for sequence extraction.
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
  A \code{\link{GRanges}} instance containing only the variants from
  \code{query} that overlap with a coding sequence in a transcript. A row is
  returned for each transcript the variant overlapped. The metadata columns
  from the \code{subject} are included with the results as well as a column 
  for 'tx_ID'. Columns for 'CodonChange', and 'AAChange' contain the coding
  changes in the format 'reference/variant'. The 'Consequence' column
  reports if the variant is synonymous, nonsynonymous,
  or a frameshift. 
}

\author{Michael Lawrence and Valerie Obenchain}

\examples{
  \dontrun{
  library(TxDb.Hsapiens.UCSC.hg18.knownGene)
  library(BSgenome.Hsapiens.UCSC.hg18)
 
  data(variants)
  txdb <- Hsapiens_UCSC_hg18_knownGene_TxDb 
  aaCoding <- predictCoding(variants, txdb, alleles="alleles",
      seqSource=Hsapiens)
 
  ## FIXME : example for fasta
  }
}

\keyword{}
