\name{summarizeVariants}
\alias{summarizeVariants}

\alias{summarizeVariants,TranscriptDb,VCF,CodingVariants-method}
\alias{summarizeVariants,TranscriptDb,VCF,FiveUTRVariants-method}
\alias{summarizeVariants,TranscriptDb,VCF,ThreeUTRVariants-method}
\alias{summarizeVariants,TranscriptDb,VCF,SpliceSiteVariants-method}
\alias{summarizeVariants,TranscriptDb,VCF,IntronVariants-method}
\alias{summarizeVariants,TranscriptDb,VCF,PromoterVariants-method}
\alias{summarizeVariants,GRangesList,VCF,ANY-method}

\title{Summarize variants}

\description{Summarize variant counts within samples grouped by gene or transcript}

\usage{
\S4method{summarizeVariants}{TranscriptDb,VCF,CodingVariants}(subject, query, mode, ...)
\S4method{summarizeVariants}{TranscriptDb,VCF,FiveUTRVariants}(subject, query, mode, ...)
\S4method{summarizeVariants}{TranscriptDb,VCF,ThreeUTRVariants}(subject, query, mode, ...)
\S4method{summarizeVariants}{TranscriptDb,VCF,SpliceSiteVariants}(subject, query, mode, ...)
\S4method{summarizeVariants}{TranscriptDb,VCF,IntronVariants}(subject, query, mode, ...)
\S4method{summarizeVariants}{TranscriptDb,VCF,PromoterVariants}(subject, query, mode, ...)
\S4method{summarizeVariants}{GRangesList,VCF,ANY}(subject, query, mode, ..., subjectFactor=factor(seq_len(length(subject))))
}

\arguments{
  \item{query}{A \linkS4class{VCF} object containing the variants.
  }
  \item{subject}{A \link[GenomicFeatures]{TranscriptDb} or \code{GRangesList}
    object that serves as the annotation. GFF files can be converted to
    \link[GenomicFeatures]{TranscriptDb} objects with 
    \code{makeTranscriptDbFromGFF()} in the \code{GenomicFeatures} package.
  }
  \item{mode}{When \code{mode} is an instance of \code{CodingVariants},
    \code{IntronVariants}, \code{FiveUTRVariants}, \code{ThreeUTRVariants},
    \code{SpliceSiteVariants} or \code{PromoterVariants}, the 
    \code{locateVariants} function will be used for counting. 
    \code{AllVariants()} and \code{IntergenicVariants} types are not 
    supported. See ?\code{locateVariants} for more detail on the variant 
    classes.

    Other counting methods can be supplied as the \code{mode} argument
    with a \code{subjectFactor} indicating how the counting should be
    summarized. Currently only \code{mode=findOverlaps} is suppported.
  }
  \item{subjectFactor}{A \code{factor} specifying how the annotation 
    should be partitioned for counting. The factor can represent genes, 
    transcripts or any other desired grouping. Must be the same length as
    \code{subject}. If NAs are present in \code{subjectFactor} they
    will be dropped along with the corresponding elements of the 
    annotation.
  }
  \item{\dots}{Additional arguments passed to methods such as
    \describe{
      \item{ignore.strand}{A \code{logical} indicating if strand should be
        igored when performing overlaps.
      }
    }
  }
}

\details{
  \code{summarizeVariants} uses genotype information in a VCF
  file to summarize variants by a factor for each sample. Factor 
  levels can represent genes, transcripts, or any desired grouping.
  Variants are counted against the annotation and summarized for 
  samples that are positive for the variant genotype.

  Currently two counting methods are available, \code{locateVariants}
  and \code{findOverlaps}. 

  \itemize{
    \item{Counting with locateVariants() :}{

      \code{subject} can be a \code{TranscriptDb} or \code{GRangesList}.
      In the case of the \code{TranscriptDb}, a region-by-transcript
      \code{GRangesList} is extracted internally and used as the annotation. 
      By default, counts are summarized by the outer list elements of the
      \code{GRangesList}. If provided, counts will be summarized by a 
      \code{subjectFactor}. \code{mode} can be any of the VariantType classes 
      defined for \code{locateVariants} except \code{IntergenicVariants} and 
      \code{AllVariants}. 

      When \code{subject} is a \code{GRangesList}, counts are summarized by
      the outer list elements of the \code{GRangesList} unless a different
      \code{subjectFactor} is provided. The \code{GRangesList} must be 
      compatible with the \code{mode} argument as follows,
      \describe{
        \item{CodingVariants :}{coding (CDS) by transcript}
        \item{IntronVariants :}{introns by transcript}
        \item{FiveUTRVariants :}{five prime UTR by transcript}
        \item{ThreeUTRVariants :}{three prime UTR by transcript}
        \item{SpliceSiteVariants :}{introns by transcript}
        \item{PromoterVariants :}{list of transcripts}
      }
    }
 
    \item{Counting with findOverlaps() :}{

      \code{subject} must be a \code{GRangesList} and \code{mode} must
      be findOverlaps (not a character string). Counts are summarized
      by the outer list elements of the \code{GRangesList} unless a
      \code{subjectFactor} is provided.
    }
  }
}

\value{
  A \code{SummarizedExperiment} object with count summaries in the 
  \code{assays} slot. The \code{rowData} contains the annotation
  corresponding to the \code{subjectFactor} provided. By default
  the \code{subjectFactor} is the length of the \code{GRangesList}
  and the origional annotation would is returned. If the \code{subjectFactor}
  provided indicates a different grouping, the annotation is relisted
  to reflect this grouping. Information in \code{colData} and 
  \code{exptData} are copied from the VCF file.
}

\author{Valerie Obenchain <vobencha@fhcrc.org>}

\seealso{
  \code{\link{readVcf}},
  \code{\link{predictCoding}}
  \code{\link{locateVariants}}
  \code{\link{summarizeVariants}},
}

\examples{
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene 

  ## Read variants from VCF and adjust seqlevels
  fl <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
  vcf <- readVcf(fl, "hg19")
  newnames <- paste("chr", seqlevels(vcf), sep="")
  names(newnames) <- seqlevels(vcf)
  vcf <- renameSeqlevels(vcf, newnames) 
  ## Confirm
  intersect(seqlevels(vcf), seqlevels(txdb))

  ## --------------------------------
  ## Counting with locateVariants()
  ## --------------------------------
  ## When 'subject' is a TranscriptDb a region-by-transcript
  ## GRangesList is extracted internally as the annotation.
  ## By default counts are summarized by the outer list elements 
  ## of the GRangesList.
  coding1 <- summarizeVariants(txdb, vcf, CodingVariants())
  colSums(assays(coding1)$counts)
  ## Dimention of the annotation used in counting
  length(rowData(coding1))

  ## When a 'subjectFactor' of gene ids is provided counting
  ## is summarized by gene. Transcripts with no gene id 
  ## (i.e, factor level is NA) are dropped and a warning is issued.
  cdsbytx <- cdsBy(txdb, "tx")
  geneid <- select(txdb, names(cdsbytx), "GENEID", "TXID")$GENEID
  coding2 <- summarizeVariants(txdb, vcf, CodingVariants(), 
                              subjectFactor=geneid)
  ## The dimention of the annotation used in counting matches the 
  ## number of non-NA gene levels.
  length(na.omit(unique(geneid))) 
  length(rowData(coding2))

  ## GRangesList can be used instead of the TranscriptDb
  coding3 <- summarizeVariants(cdsbytx, vcf, CodingVariants(), 
                               subjectFactor=geneid)
  stopifnot(identical(assays(coding2)$counts, assays(coding3)$counts))

  ## Promoter region variants summarized by transcript
  tx <- transcripts(txdb, columns="gene_id")
  txlst <- splitAsList(tx, seq_len(length(tx)))
  promoter <- summarizeVariants(txlst, vcf, 
                                PromoterVariants(upstream=100, downstream=10))
  colSums(assays(promoter)$counts)

  ## --------------------------------
  ## Counting with findOverlaps()
  ## --------------------------------
  ## 'subject' is a GRangesList. Counts are summarized
  ## by the outer list elements if no 'subjectFactor' is 
  ## provided.

  ## Summarize by transcript
  txolap <- summarizeVariants(txlst, vcf, findOverlaps)
  colSums(assays(txolap)$counts)

  ## Summarize by gene
  myfactor <- factor(as.character(mcols(tx)$gene_id)) 
  geneolap <- summarizeVariants(txlst, vcf, findOverlaps, subjectFactor=myfactor)
  colSums(assays(geneolap)$counts)
}

\keyword{methods}
