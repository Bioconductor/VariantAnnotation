\name{predictCoding}
\alias{predictCoding}
\alias{predictCoding,Ranges,TranscriptDb,ANY,DNAStringSet-method}
\alias{predictCoding,GRanges,TranscriptDb,ANY,DNAStringSet-method}
\alias{predictCoding,VCF,TranscriptDb,ANY,missing-method}

\title{Predict amino acid coding changes for variants}

\description{
  Predict amino acid coding changes for variants that fall in a coding region.
}

\usage{
  \S4method{predictCoding}{VCF,TranscriptDb,ANY,missing}(query, subject, 
    seqSource, varAllele, ...)
  \S4method{predictCoding}{Ranges,TranscriptDb,ANY,DNAStringSet}(query, subject, 
    seqSource, varAllele, ...)
  \S4method{predictCoding}{GRanges,TranscriptDb,ANY,DNAStringSet}(query, subject, 
    seqSource, varAllele, ...)
}

\arguments{
  \item{query}{A \linkS4class{VCF}, \link{Ranges} or \link{GRanges} 
    instance containing the variants to be annotated. If a \link[IRanges]{Ranges} 
    is provided it will be coerced to a \link{GRanges}. If a \linkS4class{VCF}
    is provided the \code{GRanges} from the \code{rowData} slot will be used. All 
    elementMetadata columns are ignored.
  }
  \item{subject}{A \link[GenomicFeatures]{TranscriptDb}. 
  } 
  \item{seqSource}{A \code{\link[BSgenome]{BSgenome}} instance or an \link{FaFile}
    to be used for sequence extraction.
  }
  \item{varAllele}{A \link[Biostrings]{DNAStringSet} containing the variant
    (alternate) alleles. The length of \code{varAllele} must equal the length
    of \code{query}. Missing values are represented by a zero width empty 
    element. Ranges with missing values will be ignored. 

    When the \code{query} is a \code{VCF} object the \code{varAllele} will
    be missing. This value is taken internally from the \code{VCF} with
    \code{values(alt(<VCF>))[["ALT"]]}.
  }
  \item{\dots}{Additional arguments passed to methods.
  }
}

\details{
  This function returns the amino acid coding for variants that fall in
  coding regions. The reference sequences are not required as input, rather 
  the ranges in the \code{query} must represent the width of the reference 
  allele. For guidance on how to represent an insertion, deletion or
  substitution see the 1000 Genomes VCF format (references). Reference 
  sequences are retrieved from a fasta file or \link[BSgenome]{BSgenome}. 
  Variant alleles are substituted into the reference sequences and transcribed. 
  Sequences are transcribed only if the substitution, insertion or deletion 
  results in a new sequence length divisible by 3.

  When the \code{query} is a \code{VCF} object the \code{varAllele} will
  be missing. This value is taken internally from the \code{VCF} with
  \code{values(alt(<VCF>))[["ALT"]]}. 
} 

\value{
  A \link[IRanges]{DataFrame} with a row for each variant-transcript 
  match. This includes only variants that fell in a coding regions.

  Columns include \code{queryID}, \code{consequence}, \code{refSeq}, 
  \code{varSeq}, \code{refAA}, \code{varAA}, \code{txID}, \code{geneID},
  and \code{cdsID}.

  The \code{queryHits} column provides a map back to the variants in the 
  original \code{query}. \code{refSeq} and \code{varSeq} contain the reference 
  and variant-modified DNA sequences. Reference and variant amino acid codes are
  provided in \code{refAA} and \code{varAA}. 

  Possible values for the \code{consequence} column are `synonymous',
  `nonsynonymous', or `frameshift'. Variant sequences are transcribed only if they
  are divisible by 3. When a sequence is not transcribed the \code{varAA} column 
  will be empty and the \code{consequence} will be a `frameshift'. 
}

\references{
  \url{http://vcftools.sourceforge.net/specs.html}
}

\author{Michael Lawrence and Valerie Obenchain <vobencha@fhcrc.org>}

\seealso{
  \link{readVcf}
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

  ## rename seqlevels in the VCF to match those in the TxDb
  rd <- renameSeqlevels(rowData(vcf), c("22"="chr22"))
  rowData(vcf) <- rd 
  aacoding <- predictCoding(query=vcf, subject=txdb19, seqSource=Hsapiens)

  ## Consequence for all variant-transcript matches
  table(aacoding$consequence)
}

\keyword{methods}
