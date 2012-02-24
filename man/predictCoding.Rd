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

    When the \code{query} is not a \code{VCF} object a \code{varAllele} must be 
    provided. The \code{varAllele} must be a \code{DNAStringSet} the same length
    as the \code{query}. If there are multiple alternate alleles per variant
    the \code{query} must be expanded to one row per each alternate allele.
    See examples.
  }
  \item{subject}{A \link[GenomicFeatures]{TranscriptDb}. 
  } 
  \item{seqSource}{A \code{\link[BSgenome]{BSgenome}} instance or an \link{FaFile}
    to be used for sequence extraction.
  }
  \item{varAllele}{A \link[Biostrings]{DNAStringSet} containing the variant
    (alternate) alleles. The length of \code{varAllele} must equal the length
    of \code{query}. Missing values are represented by a zero width empty 
    element. Ranges with missing \code{varAllele} values are ignored; those with 
    an \sQuote{N} character are not translated.

    When the \code{query} is a \code{VCF} object the \code{varAllele} argument 
    will be missing. This value is taken internally from the \code{VCF} with
    \code{values(alt(<VCF>))[["ALT"]]}.
  }
  \item{\dots}{Additional arguments passed to methods.
  }
}

\details{
  This function returns the amino acid coding for variants that fall in
  coding regions. The reference sequences are taken from a fasta file or 
  \link[BSgenome]{BSgenome}. The width of the reference is determined
  from the start postion and width of the range in the \code{query}.
  For guidance on how to represent an insertion, deletion or substitution 
  see the 1000 Genomes VCF format (references). 

  Variant alleles are taken from the \code{varAllele} when supplied.
  When the \code{query} is a \code{VCF} object the \code{varAllele} will
  be missing. This value is taken internally from the \code{VCF} with
  \code{values(alt(<VCF>))[["ALT"]]}. The variant allele is substituted 
  into the reference sequences and transcribed. Transcription only 
  occurs if the substitution, insertion or deletion results in a new sequence 
  length divisible by 3.
} 

\value{
  A \link[IRanges]{DataFrame} with a row for each variant-transcript 
  match. This includes only variants that fell in a coding regions.

  Columns include \code{queryID}, \code{consequence}, \code{refSeq}, 
  \code{varSeq}, \code{refAA}, \code{varAA}, \code{txID}, \code{geneID},
  and \code{cdsID}.

  \describe{
    \item{\code{queryID}}{
      The \code{queryID} column provides a map back to the variants in the 
      original \code{query}. If the \code{query} was a \code{VCF} object this
      index corresponds to the row in the \code{GRanges} in the \code{rowData}
      slot. If \code{query} was a \code{GRanges}, \code{RangedData} or
      \code{RangesList} the index corresponds to the row in the object.
    }
    \item{\code{consequence}}{
      Possible values for the \code{consequence} column are `synonymous',
      `nonsynonymous', or `frameshift'. Variant sequences are transcribed only 
      if they are divisible by 3. When a sequence is not transcribed the 
      \code{varAA} column will be empty and the \code{consequence} will be a 
      `frameshift'. 
    }
    \item{\code{refSeq}}{
      This is the reference codon sequence. This range is typically greater
      than the width of the reference sequence because it includes all codons
      involved in the sequence modification. If the reference sequence is of
      width 2 but the alternate allele is of width 4 then at least two codons
      must be included in the \code{refSeq}.
    }
    \item{\code{varSeq}}{
      This sequence is the result of replacing the reference range with the
      alternate allele. If the result of this substitution is not a multiple
      of 3 no translation is performed and the \code{varAA} value will be missing.
    }
    \item{\code{refAA}}{
      The reference amino acid column contains the translated \code{refSeq}.
    }
    \item{\code{varAA}}{
      The variant amino acid column contains the translated \code{varSeq}. When
      translation is not possible this value is missing.
    }
    \item{\code{txID}}{
      The transcript id taken from the \code{TranscriptDb} object.
    }
    \item{\code{geneID}}{
      The gene id taken from the \code{TranscriptDb} object.
    }
    \item{\code{cdsID}}{
      The coding sequence id taken from the \code{TranscriptDb} object.
    }
  }
}

\references{
  \url{http://vcftools.sourceforge.net/specs.html}
}

\author{Michael Lawrence and Valerie Obenchain <vobencha@fhcrc.org>}

\seealso{
  \link{readVcf},
  \link{locateVariants},
  \link{getTranscriptSeqs}
}

\examples{
  library(BSgenome.Hsapiens.UCSC.hg19)
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  txdb19 <- TxDb.Hsapiens.UCSC.hg19.knownGene 

  ## Read variants from a VCF file 
  fl <- system.file("extdata", "ex1.vcf", package="VariantAnnotation")
  vcf <- readVcf(fl, "hg19")

  ## Rename seqlevels in the VCF object to match those in the TxDb
  vcf_adj <- renameSeqlevels(vcf, c("16"="chr16", "21"="chr21"))
  ## Confirm common seqlevels
  intersect(seqlevels(vcf_adj), seqlevels(txdb19))

 
  ## Note that when 'query' is a VCF object the varAllele argument is missing.
  coding <- predictCoding(vcf_adj, txdb19, Hsapiens)

  ## The queryID column in the result corresponds to the row numbers of the
  ## of the original 'query' which was the VCF. To map the result from 
  ## predictCoding() back to the VCF we use the rowData.
  var <- rowData(vcf_adj)[coding$queryID]

  ## Consequence for all variant-transcript matches
  table(coding$consequence)

  ## Alternatively, we can use a GRanges as the 'query' to predictCoding().
  ## The seqlevels have been adjusted in the VCF object so the GRanges we
  ## extract has the correct seqlevels.
  rd <- rowData(vcf_adj)
 
  ## The GRanges provided as the 'query' must have one row for each alternate
  ## allele. We see we that variants 1, 2 and 10 have two alternate alleles.
  altallele <- values(alt(vcf_adj))[["ALT"]]
  eltlen <- elementLengths(altallele)
  eltlen

  ## Expand the GRanges by the element lengths
  rd_exp <- rep(rd, eltlen)
 
  ## Call predictCoding() with the expanded GRanges as the 'query'
  ## and the unlisted alternate allele as the 'varAllele'
  ## Note that the queryID values in coding2 correspond to the rows
  ## in the expanded GRanges. 
  coding2 <- predictCoding(query=rd_exp, subject=txdb19, seqSource=Hsapiens, 
                           varAllele=unlist(altallele)) 
}

\keyword{methods}
