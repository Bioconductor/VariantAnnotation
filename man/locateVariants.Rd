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
      Variants that fall within a transcript are classified as \sQuote{coding},
      \sQuote{intron}, \sQuote{3'UTR} or \sQuote{5'UTR}. All overlap operations 
      require variants to fall completely \sQuote{within} the defined region to 
      be classified as such. When a variant falls within a transcript but not 
      completely within one of the defined regions it is designated 
      \sQuote{transcript_region}. Variants that do not fall in a transcript are 
      classified as \sQuote{intergenic}. 

      A variant location will be \sQuote{NA} if the variant was called on a 
      chromosome present in the \code{query} but not in the \code{subject}
    }
    \item{Values of txID, geneID and cdsID :}{
      Variants that fall in a trancscript regions will have values for all of 
      \code{txID}, \code{geneID} and \code{cdsID} if they are present in the 
      annotation. \sQuote{intergenic} variants will have \code{geneID} values 
      for preceding and following genes and \code{NA} for \code{txID} and 
      \code{cdsID}. All ID values will be NA for variants with a location of
      \code{transcript_region} or \code{NA}.
    }
  }
}

\value{
  A \link[IRanges]{DataFrame} with a row for each variant-transcript match.
  Columns are \code{queryID}, \code{location}, \code{txID}, \code{geneID} and 
  \code{cdsID}. 

 \describe{
    \item{\code{queryHits}}{
      The \code{queryHits} column provides a map back to the row in the
      original \code{query}. If the \code{query} was a \code{VCF} object this
      index corresponds to the row in the \code{GRanges} in the \code{rowData}
      slot. 
    }
    \item{\code{location}}{
      Possible values for the \code{consequence} column are \sQuote{synonymous},
      \sQuote{nonsynonymous}, or \sQuote{frameshift}. Variant sequences are 
      transcribed only if they are divisible by 3. When a sequence is not 
      transcribed the \code{varAA} column will be empty and the 
      \code{consequence} will be a \sQuote{frameshift}.
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
  \code{\link{readVcf}},
  \code{\link{predictCoding}}
}

\examples{
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  txdb19 <- TxDb.Hsapiens.UCSC.hg19.knownGene 

  ## Read variants from a VCF file 
  fl <- system.file("extdata", "ex1.vcf", package="VariantAnnotation")
  vcf <- readVcf(fl, genome="hg19")

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

  ## Alternatively, a VCF object can be provided as the query.
  ## First, the seqlevels in the VCF must be modified. 
  ## Adjusting the seqlevels at the VCF-level ensures all future data 
  ## extractions from this VCF object have the modified seqlevels.
  ##
  ## vcf_adj <- renameSeqlevels(vcf, c("16"="chr16", "21"="chr21")) 
  ## loc <- locateVariants(vcf_adj, txdb19)
}

\keyword{methods}
