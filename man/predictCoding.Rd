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
  \item{query}{A \link[IRanges]{Ranges} or \link[GenomicRanges]{GRanges} object
    containing the variants to be annotated. If a \link[IRanges]{Ranges}
    object is provided it will be coerced to a \link{GRanges} object.
  }
  \item{subject}{A character vector, \link[GenomicFeatures]{TranscriptDb} or
    \link[GenomicRanges]{GRangesList} object. 

    When \code{subject} is a character vector it should contain the complete
    path to the bam files to be used for sequence extraction. 
    If \code{subject} is a \link{TranscriptDb}, the coding
    regions will be identified using the 'cdsBy' function on this object. 
    It is assumed that when \code{subject} is provided as a \link{GRangesList} 
    object it was created with the 'cdsBy' function on a \link{TranscriptDb}
    object.
  } 
  \item{alleles}{A character representing the name of the column 
    in \code{query} that contains the reference and variant alleles. The
    data in the column must be character have format 'reference/variant'.
    For example, 'A/G' would indicate the reference allele is A and the 
    variant is G. An two nucleotide insertion would be represented by
    '-/AT', and a four nucleotide deletion would be 'GCCT/-'. 
  }
  \item{seqSource}{A \code{\link[BSgenome]{BSgenome}} object or \link{FaFileList}.
  }
  \item{\dots}{Additional arguments
  }
}

\details{
  Sequences for the variants in \code{query} that fall in a coding region
  can be extracted from fasta files or a \code{\link[BSgenome]{BSgenome}} object. The
  sequences are translated and amino acid coding changes are returned. 
} 

\value{
  A \code{\link{GRanges}} object containing only the variants from
  \code{query} that overlap with a coding sequence in a transcript. A row is
  returned for each transcript the variant overlapped. The metadata columns
  from the \code{subject} are included with the results as well as a column 
  for 'tx_ID'. Columns for 'CodonChange', and 'AAChange' contain the coding
  changes in the format 'reference/variant'. The 'Consequence' column
  indicates if the variant is synonymous or nonsynonymous coding,
  or a frameshift. Amino acid codes are currently computed for snvs and
  substitutions with the appropriate number of variant basepairs but not 
  frameshifts. 
}

\references{
}
\author{Michael Lawrence and Valerie Obenchain}


\examples{
  \dontrun{
  library(TxDb.Hsapiens.UCSC.hg18.knownGene)
  library(BSgenome.Hsapiens.UCSC.hg18)
 
  data(variants)
  strand(variants) <- "+"
  txdb <- Hsapiens_UCSC_hg18_knownGene_TxDb 
  aaCoding <- predictCoding(variants, txdb, alleles="alleles",
      seqSource=Hsapiens)

  ## FIXME : example for fasta
  }
}

\keyword{}
