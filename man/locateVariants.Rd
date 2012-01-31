\name{locateVariants}
\alias{locateVariants}
\alias{locateVariants,Ranges,TranscriptDb-method}
\alias{locateVariants,GRanges,TranscriptDb-method}
\alias{locateVariants,VCF,TranscriptDb-method}

\title{Locate variants}

\description{Variant location with respect to gene functionality}

\usage{
  locateVariants(query, subject, ...)
  \S4method{locateVariants}{Ranges,TranscriptDb}(query, subject, ...)
  \S4method{locateVariants}{GRanges,TranscriptDb}(query, subject, ...)
  \S4method{locateVariants}{VCF,TranscriptDb}(query, subject, ...)
}

\arguments{
  \item{query}{A \link[IRanges]{Ranges}, \link[GenomicRanges]{GRanges} or 
    \linkS4class{VCF} object containing the variants. Metadata columns are 
    allowed but are ignored. 
  }
  \item{subject}{A \link[GenomicFeatures]{TranscriptDb} object that serves
    as the annotation reference.
  }
  \item{\dots}{Additional arguments passed to methods.
  }
}

\details{
  \describe{
    \item{Range representation :}{
      \code{locateVariants} accepts a \link[GenomicRanges]{GRanges} object that 
      contains the variants of interest. The ranges should reflect the
      position(s) of the reference allele. For snps the range will be of width 1. 
      For range insertions or deletions the reference allele could be a sequence 
      such as GGTG in which case the width of the range should be 4. 
    }
    \item{Location classification :}{
      Variants that fall within a transcript are classified as 'coding',
      'intron', '3'UTR' or '5'UTR'. All overlap operations require variants to 
      fall completely 'within' the defined region to be classified as such. When 
      a variant falls within a transcript but not completely within one of the 
      defined regions it is designated 'transcript_region'.

      Variants that do not fall in a transcript are classified as 'intergenic' 
      if flanking genes can be identified. A variant location will be 'NA' if 
      it meets one of the following conditions :

      (1) the variant does not hit a transcript and no flanking genes can be
      identified

      (2) the variant was called on a chromosome present in the \code{query} but
      not in the \code{subject}
    }
    \item{Values of txID, geneID and cdsID :}{
      Variants that fall in a trancscript ('coding', '5'UTR', '3'UTR' and 
      'intron') will have values for all of \code{txID}, \code{geneID} and 
      \code{cdsID} if they are present in the annotation. 'intergenic' variants 
      will have \code{geneID} values for preceding and following genes and NA
      for \code{txID} and \code{cdsID}. All ID values will be NA for variants 
      with a location of 'transcript_region' or NA.
    }
  }
}

\value{
  A \link[IRanges]{DataFrame} with a row for each variant-transcript match.
  Columns are \code{queryID}, \code{location}, \code{txID}, \code{geneID} and 
  \code{cdsID}. 

 \describe{
    \item{\code{queryHits}}{
      The \code{queryHits} column provides a map back to the variants in the
      original \code{query}. If the \code{query} was a \code{VCF} object this
      index corresponds to the row in the \code{GRanges} in the \code{rowData}
      slot. 
    }
    \item{\code{location}}{
      Possible values for the \code{consequence} column are `synonymous',
      `nonsynonymous', or `frameshift'. Variant sequences are transcribed only
      if they are divisible by 3. When a sequence is not transcribed the
      \code{varAA} column will be empty and the \code{consequence} will be a
      `frameshift'.
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

\author{Valerie Obenchain <vobencha@fhcrc.org>}

\seealso{
  \code{\link{readVcf}}
  \code{\link{predictCoding}}
}

\examples{
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  txdb19 <- TxDb.Hsapiens.UCSC.hg19.knownGene 

  ## Read variants from a VCF file 
  fl <- system.file("extdata", "ex1.vcf", package="VariantAnnotation")
  vcf <- readVcf(fl, "hg19")

  ## Note the variant seqlevels do not match those in the TxDb
  rd <- rowData(vcf)
  head(seqlevels(rd))
  head(seqlevels(txdb19))
  intersect(seqlevels(rd), seqlevels(txdb19))

  ## Use the renameSeqlevels() helper to rename seqlevels
  rd_adj <- renameSeqlevels(rd, c("16"="chr16", "21"="chr21"))
  intersect(seqlevels(rd_adj), seqlevels(txdb19))

  ## Call locateVariants() with GRanges as the query 
  loc <- locateVariants(rd_adj, txdb19)

  ## Location for all variant-transcript matches
  table(loc$Location)

  ## Alternatively, we can provide a VCF object as the 'query'.
  ## First, the seqlevels in the VCF must be modified. 
  ## Adjusting the seqlevels at the VCF-level will ensure all future data 
  ## extractions from this VCF object have the modified seqlevels.
  ##
  ## vcf_adj <- renameSeqlevels(vcf, c("16"="chr16", "21"="chr21")) 
  ## loc <- locateVariants(vcf_adj, txdb19)
}

\keyword{methods}
