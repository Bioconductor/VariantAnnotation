\name{predictCoding}
\alias{predictCoding}
\alias{predictCoding,Ranges,TxDb,ANY,DNAStringSet-method}
\alias{predictCoding,GRanges,TxDb,ANY,DNAStringSet-method}
\alias{predictCoding,CollapsedVCF,TxDb,ANY,missing-method}
\alias{predictCoding,ExpandedVCF,TxDb,ANY,missing-method}
\alias{predictCoding,VRanges,TxDb,ANY,missing-method}

\title{Predict amino acid coding changes for variants}

\description{
  Predict amino acid coding changes for variants a coding regions
}

\usage{
\S4method{predictCoding}{CollapsedVCF,TxDb,ANY,missing}(query, subject, seqSource, varAllele, ..., ignore.strand=FALSE)
\S4method{predictCoding}{ExpandedVCF,TxDb,ANY,missing}(query, subject, seqSource, varAllele, ..., ignore.strand=FALSE)
\S4method{predictCoding}{Ranges,TxDb,ANY,DNAStringSet}(query, subject, seqSource, varAllele, ..., ignore.strand=FALSE)
\S4method{predictCoding}{GRanges,TxDb,ANY,DNAStringSet}(query, subject, seqSource, varAllele, ..., ignore.strand=FALSE)
\S4method{predictCoding}{VRanges,TxDb,ANY,missing}(query, subject, seqSource, varAllele, ..., ignore.strand=FALSE)
}

\arguments{
  \item{query}{A \linkS4class{VCF}, \link{Ranges}, \link{GRanges} or
    \code{VRanges} instance containing the variants to be annotated. If
    \code{query} is a \link[IRanges]{Ranges} or \code{VRanges} it is coerced to
    a \link{GRanges}.  If a \linkS4class{VCF} is provided the \code{GRanges}
    returned by the \code{rowRanges()} accessor will be used. All metadata
    columns are ignored.

    When \code{query} is not a \code{VCF} object a \code{varAllele} must be
    provided. The \code{varAllele} must be a \code{DNAStringSet} the same length
    as the \code{query}. If there are multiple alternate alleles per variant the
    \code{query} must be expanded to one row per each alternate allele.  See
    examples.

    NOTE: Variants are expected to conform to the VCF specs as described
    on the 1000 Genomes page (see references). Indels must include
    the reference allele; zero-width ranges are not valid and return no matches.
  }
  \item{subject}{A \link[GenomicFeatures]{TxDb} object that serves 
    as the annotation. GFF files can be converted to 
    \link[GenomicFeatures]{TxDb} objects with 
    \code{makeTxDbFromGFF()} in the \code{GenomicFeatures} package.
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
    \code{alt(<VCF>)}.
  }
  \item{\dots}{Additional arguments passed to methods. Arguments
    \code{genetic.code} and \code{if.fuzzy.codon} are supported for the
    \code{translate} function.
  }
  \item{ignore.strand}{A \code{logical} indicating if strand should be ignored
    when performing overlaps.

    When \code{ignore.strand == TRUE} the \code{query} strand is set to
    '*' and can overlap with any strand of the \code{subject}. The
    return \code{GRanges} reflects the strand of the \code{subject} hit,
    however, the positions and alleles reported are computed as if 
    both were from the '+' strand.

    \code{ignore.strand == FALSE} requires the \code{query} and \code{subject}
    to have compatible strand. Again the return \code{GRanges} reports
    the strand of the \code{subject} hit but in this case positions and
    alleles are computed according to the strand of the \code{subject}.
  }
}

\details{
  This function returns the amino acid coding for variants that fall 
  completely `within' a coding region. The reference sequences are 
  taken from a fasta file or \link[BSgenome]{BSgenome}. The width of 
  the reference is determined from the start position and width of the 
  range in the \code{query}. For guidance on how to represent an insertion, 
  deletion or substitution see the 1000 Genomes VCF format (references). 

  Variant alleles are taken from the \code{varAllele} when supplied.
  When the \code{query} is a \code{VCF} object the \code{varAllele} will
  be missing. This value is taken internally from the \code{VCF} with
  \code{alt(<VCF>)}. The variant allele is substituted 
  into the reference sequences and transcribed. Transcription only 
  occurs if the substitution, insertion or deletion results in a new sequence 
  length divisible by 3.

  When the \code{query} is an unstranded (*) \code{GRanges} \code{predictCoding} 
  will attempt to find overlaps on both the positive and negative strands of the
  \code{subject}. When the subject hit is on the negative strand the return 
  \code{varAllele} is reverse complemented. The strand of the returned 
  \code{GRanges} represents the strand of the subject hit.
} 

\value{
  A \link[GenomicRanges]{GRanges} with a row for each variant-transcript 
  match. The result includes only variants that fell within coding regions.
  The strand of the output \code{GRanges} represents the strand of the 
  \code{subject} hit.

  At a minimum, the metadata columns (accessible with \code{mcols}) include,
  \describe{
    \item{\code{varAllele}}{
      Variant allele. This value is reverse complemented for an unstranded 
      \code{query} that overlaps a \code{subject} on the negative strand.
    }
    \item{\code{QUERYID}}{
      Map back to the row in the original query
    }
    \item{\code{TXID}}{
      Internal transcript id from the annotation
    }
    \item{\code{CDSID}}{
      Internal coding region id from the annotation
    }
    \item{\code{GENEID}}{
      Internal gene id from the annotation
    }
    \item{\code{CDSLOC}}{
      Variant location in coding region-based coordinates. This position is 
      relative to the start of the coding (cds) region defined in the 
      \code{subject} annotation.
    }
    \item{\code{PROTEINLOC}}{
      Variant codon triplet location in coding region-based coordinates.
      This position is relative to the start of the coding (cds) region
      defined in the \code{subject} annotation. 
    }
    \item{\code{CONSEQUENCE}}{
      Possible values are `synonymous', `nonsynonymous', `frameshift', 
      `nonsense' and `not translated'. Variant sequences are translated only 
      when the codon sequence is a multiple of 3. The value will be `frameshift' 
      when a sequence is of incompatible length. `not translated' is used
      when \code{varAllele} is missing or there is an \sQuote{N} in the 
      sequence. `nonsense' is used for premature stop codons.
    }
    \item{\code{REFCODON}}{
      The reference codon sequence. This range is typically greater
      than the width of the range in the \code{GRanges} because it includes 
      all codons involved in the sequence modification. If the reference 
      sequence is of width 2 but the alternate allele is of width 4 then at 
      least two codons must be included in the \code{REFSEQ}.
    }
    \item{\code{VARCODON}}{
      This sequence is the result of inserting, deleting or replacing the 
      position(s) in the reference sequence alternate allele. If the result 
      of this substitution is not a multiple of 3 it will not be translated.
    }
    \item{\code{REFAA}}{
      The reference amino acid column contains the translated \code{REFSEQ}.
    }
    \item{\code{VARAA}}{
      The variant amino acid column contains the translated \code{VARSEQ}. When
      translation is not possible this value is missing.
    }
  }
}

\references{
  \url{http://www.1000genomes.org/wiki/analysis/variant-call-format/}
  \url{http://vcftools.sourceforge.net/specs.html}
}

\author{Michael Lawrence and Valerie Obenchain <vobencha@fhcrc.org>}

\seealso{
  \link{readVcf},
  \link{locateVariants},
  \link{refLocsToLocalLocs}
  \link{getTranscriptSeqs}
}

\examples{
  library(BSgenome.Hsapiens.UCSC.hg19)
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene 

  ## ----------------------------
  ## VCF object as query 
  ## ----------------------------
  ## Read variants from a VCF file 
  fl <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
  vcf <- readVcf(fl, "hg19")

  ## Rename seqlevels in the VCF object to match those in the TxDb.
  vcf <- renameSeqlevels(vcf, "chr22")
  ## Confirm common seqlevels
  intersect(seqlevels(vcf), seqlevels(txdb))
 
  ## When 'query' is a VCF object the varAllele argument is missing.
  coding1 <- predictCoding(vcf, txdb, Hsapiens)
  head(coding1, 3)

  ## Exon-centric or cDNA locations:
  exonsbytx <- exonsBy(txdb, "tx")
  cDNA <- mapToTranscripts(coding1[!duplicated(ranges(coding1))], exonsbytx)
  coding1$cDNA <- ranges(cDNA)[togroup(coding1$QUERYID)]

  ## ----------------------------
  ## GRanges object as query 
  ## ----------------------------
  ## A GRanges can also be used as the 'query'. The seqlevels in the VCF
  ## were adjusted in previous example so the GRanges extracted with
  ## has the correct seqlevels.
  rd <- rowRanges(vcf)
 
  ## The GRanges must be expanded to have one row per alternate allele. 
  ## Variants 1, 2 and 10 have two alternate alleles.
  altallele <- alt(vcf)
  eltlen <- elementLengths(altallele)
  rd_exp <- rep(rd, eltlen)
 
  ## Call predictCoding() with the expanded GRanges and the unlisted
  ## alternate allele as the 'varAllele'.
  coding2 <- predictCoding(rd_exp, txdb, Hsapiens, unlist(altallele)) 
}

\keyword{methods}
