\name{predictCoding}
\alias{predictCoding}

\title{Predict amino acid coding changes for single nucleotide variants}

\description{
}

\usage{
predictCoding(x, txdb, BSgenomeOrganism = Hsapiens, variantColumn = "read")
}

\arguments{
  \item{x}{A \link[IRanges]{Ranges} or \link[GenomicRanges]{GRanges} object}
  \item{txdb}{A \link[GenomicFeatures]{TranscriptDb} object}
  \item{BSgenomeOrganism}{The name of the BSgenome
    organism to be used for extracting sequences}
  \item{variantColumn}{A character vector with the name of the column
    in \code{x} that contains the observed variant allele}
}

\details{
  A \link[GenomicFeatures]{TranscriptDb} object is created by calling cdsBy on \code{txdb}.
  Overlap locations of the variants in \code{x} and the transcripts in the cds_by_tx object
  are identified and the corresponding transcript sequences are obtained
  from the \link[BSgenome]{BSgenome} object specified in \code{BSgenomeOrganism}. 
  Codons containing the variants are identified and the sequences are translated for the 
  reference and variant alleles.
} 

\value{
  A \code{\linkS4class{GRanges}} object containing only the variants from
  \code{x} that overlap with a coding sequence in a transcript. A row is
  returned for each transcript the variant overlapped. Three metadata columns
  are added, 'tx_id', 'aa', and 'aa.ref'. The 'tx_id' corresponds to the
  transcript id found in the cds_by_tx object. 'aa' and 'aa.ref' represent
  the amino acids codes for the variant and reference alelles, respectively. 
}

\references{
}
\author{M. Lawrence}

\seealso{
}

\examples{
  \dontrun{
  ## using sample data from ANNOVAR
  library(TxDb.Hsapiens.UCSC.hg18.knownGene)
  library(BSgenome.Hsapiens.UCSC.hg18)
  data(variants)

  strand(variants) <- "+"
  txdb <- Hsapiens_UCSC_hg18_knownGene_TxDb 
  BSgenomeOrganism <- Hsapiens
  aaCoding <- predictCoding(variants, txdb, BSgenomeOrganism, variantColumn="obsAllele")
  }
}

\keyword{}
