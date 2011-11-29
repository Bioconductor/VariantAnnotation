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
  \item{query}{A \link{Ranges} or \link{GRanges} instance containing the 
    variants to be annotated. If a \link[IRanges]{Ranges} is provided it will be 
    coerced to a \link{GRanges}.

    A single elementMetadata column containing the variant alleles is required.
    Others are ignored. 
  }
  \item{subject}{A \link{GRangesList} or a \link[GenomicFeatures]{TranscriptDb}. 

    When \code{subject} is provided as a \link{GRangesList} it is assumed that
    the outer list elements define coding regions and the individual elements
    are ranges within that coding region. This object can be created with the 
    'cdsBy' function on a \link[GenomicFeatures]{TranscriptDb}.
  } 
  \item{seqSource}{A \code{\link[BSgenome]{BSgenome}} instance or an \link{FaFile}
    to be used for sequence extraction.
  }
  \item{varAllele}{The character name of the elementMetadata column in 
    \code{query} that contains the variant alleles. These data should be a 
    \link[Biostrings]{DNAStringSet}. Deletions are represented by an empty
    sequence of width 0.
  }
  \item{\dots}{Additional arguments passed to methods.
  }
}

\details{
  This function identifies the variants in coding regions and retrieves
  the reference sequences from a fasta file or \link[BSgenome]{BSgenome}.
  Alleles provided in the \code{varAllele} column are substituted into the 
  reference sequences and transcribed. Variant sequences are transcribed only 
  if the substitution, insertion or deletion results in a new sequence length 
  divisible by 3. 
} 

\value{
  A \link[IRanges]{DataFrame} with a row for each variant-transcript 
  match. This includes only variants that fell in a coding regions.

  Columns include \code{queryHits}, \code{txID}, \code{refSeq}, 
  \code{varSeq}, \code{refAA}, \code{varAA}, \code{Consequence}, and any 
  metadata that was present in the \code{subject}. 

  The \code{queryHits} column provides a map back to the variants in the 
  original \code{query}. \code{refSeq} and \code{varSeq} contain the reference 
  and variant-modified DNA sequences. Reference and variant amino acid codes are
  provided in \code{refAA} and \code{varAA}. 

  Possible values for the \code{Consequence} column are `synonymous',
  `nonsynonymous', or `frameshift'. Variant sequences are transcribed only if they
  are divisible by 3. When a sequence is not transcribed the \code{varAA} column 
  will be empty and the \code{Consequence} will be a `frameshift'. 
}

\author{Michael Lawrence and Valerie Obenchain <vobencha@fhcrc.org>}

\seealso{
  \link{locateVariants}
  \link{getTranscriptSeqs}
}

\examples{
  library(BSgenome.Hsapiens.UCSC.hg19)
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  txdb19 <- TxDb.Hsapiens.UCSC.hg19.knownGene 

  ## Read variants from a VCF file and adjust seqlevels to match the annotation 
  fl <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
  vcf <- readVcf(fl, "hg19")
  rd <- rowData(vcf)
  seqlevels(rd) <- paste("chr", seqlevels(rd), sep="")

  ## Expand the ALT column
  rdflat <- rep(rd, elementLengths(values(rd)[["ALT"]]))
  values(rdflat) <- DataFrame(ALT=unlist(values(rd)[["ALT"]]))

  aaCoding <- predictCoding(
                query=rdflat, 
                subject=txdb19, 
                seqSource=Hsapiens,
                varAllele="ALT")

  ## Consequence for all variant-transcript matches
  table(aaCoding$Consequence)
}

\keyword{methods}
