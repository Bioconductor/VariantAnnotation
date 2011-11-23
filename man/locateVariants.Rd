\name{locateVariants}
\alias{locateVariants}
\alias{locateVariants,GRanges,TranscriptDb-method}

\title{Locate variants}

\description{Variant location with respect to gene functionality}

\usage{
  locateVariants(query, subject, ...)
  \S4method{locateVariants}{GRanges,TranscriptDb}(query, subject, ...)
}

\arguments{
  \item{query}{A \link[GenomicRanges]{GRanges} object containing the variants.
    Metadata columns can be present but are ignored. Only the \code{seqnames}
    \code{ranges} and \code{strand} are pertinent to the method. 
  }
  \item{subject}{A \link[GenomicFeatures]{TranscriptDb} object that serves
    as the annotation reference.
  }
  \item{\dots}{Additional arguments passed to methods.
  }
}

\details{
  Range representation :

  \code{locateVariants} accepts a \link[GenomicRanges]{GRanges} object that 
  contains the variants of interest. The ranges should reflect the position(s) 
  of the reference allele. For snps the range will be of width 1. For range
  insertions or deletions the reference allele could be a sequences such as
  GGTG in which case the width of the range should be 4. 

  Location classification :

  Variants that fall within a transcript are classified as 'coding', 'intron', 
  '3'UTR' or '5'UTR'. All overlap operations require variants to fall completely 
  within one of the defined regions to be classified. When a variant falls 
  within a transcript but not completely within one of the defined regions it 
  is designated 'transcript_region'.

  Variants that do not fall in a transcript are classified as 'intergenic' if 
  flanking genes can be identified. A variant is classified as NA if it
  meets one of the following conditions :

  (1) the variant does not hit a transcript and no flanking genes can be
  identified

  (2) the variant was called on a chromosome present in the \code{query} but
  not in the \code{subject}
}

\value{
  A \link[IRanges]{DataFrame} with a row for each variant-transcript match.
  Only variants that fall within coding regions are reported. Columns are 
  \code{queryHits}, which is a map back to the original variants in the 
  \code{query}, \code{txID}, \code{geneID} and \code{Location}. 

  Variants that fall in a trancscript will have both a \code{txID} and 
  \code{geneID} if present in the annotation. The 'intergenic' variants will 
  have no \code{txID} and the \code{geneID} will contain identifiers for the 
  preceding and following genes. Both \code{txID} and \code{geneID} will be 
  NA for variants with Location classification of 'transcript_region' and 
  NA.
}

\author{Valerie Obenchain <vobencha@fhcrc.org>}

\seealso{
  \code{\link{predictCoding}}
}

\examples{
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  txdb19 <- TxDb.Hsapiens.UCSC.hg19.knownGene 

  ## Read variants from a VCF file and adjust seqlevels to match the annotation 
  fl <- system.file("extdata", "ex1.vcf", package="VariantAnnotation")
  vcf <- readVcf(fl, "hg19")
  rd <- rowData(vcf)
  seqlevels(rd) <- paste("chr", seqlevels(rd), sep="")

  loc <- locateVariants(
           query=rd, 
           subject=txdb19)

  ## Location for all variant-transcript matches
  table(loc$Location)
}

\keyword{methods}
