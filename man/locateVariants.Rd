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
  Range representation of variants :

  \code{locateVariants} accepts a \link[GenomicRanges]{GRanges} object that 
  contains the variants of interest. The ranges should reflect the position(s) 
  of the reference allele. For snps the range will be of width 1. For range
  insertions or deletions the reference allele could be a sequences such as
  `GGTG' in which case the width of the range should be 4. 

  Location classification :

  Variants that fall in a transcript are classified as `coding', `intron', 
  `3UTR' or `5UTR'. All overlap operations require that variants fall
  completely "within" one of the defined regions to be classified. Those not 
  meeting this criteria are labeled `unknown'.

  Variants that do not overlap a transcript are classified as `intergenic' if 
  flanking genes can be identified and `noGeneMatch' if not. A case where 
  flanking genes may not be available is a \code{query} that contains variants 
  for unassembled chromosomes.
}

\value{
  A \link[IRanges]{DataFrame} with a row for each variant-transcript match.
  Columns include \code{queryHits}, which is a map back to the original variants
  in the \code{query}, \code{txID}, \code{geneID} and \code{Location}. 

  Variants that fall in a trancscript will have both a \code{txID} and 
  \code{geneID}. The `intergenic' variants will have no \code{txID} and 
  the \code{geneID} will contain identifiers for the preceding and following
  genes. Both \code{txID} and \code{geneID} will be `NA' for `noGeneMatch'.
}

\author{Valerie Obenchain <vobencha@fhcrc.org>}

\seealso{
  \code{\link{predictCoding}}
}

\examples{
  library(TxDb.Hsapiens.UCSC.hg18.knownGene)
  data(variants)
 
  txdb <- TxDb.Hsapiens.UCSC.hg18.knownGene 
  loc <- locateVariants(variants, txdb)
}

\keyword{methods}
