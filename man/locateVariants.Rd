\name{locateVariants}
\alias{locateVariants}
\alias{locateVariants,GRanges,TranscriptDb-method}

\title{Locate variants}

\description{Identify variant location with respect to gene functionality}

\usage{
  locateVariants(query, subject, ...)
  \S4method{locateVariants}{GRanges,TranscriptDb}(query, subject, ...)
}

\arguments{
  \item{query}{A \link[GenomicRanges]{GRanges} object containing the variants
  }
  \item{subject}{A \link[GenomicFeatures]{TranscriptDb} object that will serve
    as the annotation reference
  }
  \item{\dots}{Additional arguments passed to methods
  }
}

\details{
  Range representation :

  \code{locateVariants} accepts a \link[GenomicRanges]{GRanges} object that 
  contains the variants of interest. The range should reflect the positions of
  the reference basepairs being substituted or deleted. The width will be the
  same as the number of basepairs being substituted or deleted. For insertions,
  the width is negative and the range is represented as \code{start} =
  \code{end} + 1. For examples of how different variant types should
  be represented see the \code{variants} sample data in the package.
 

  Location :

  Variant positions are overlaped with feature ranges
  of the \link[GenomicFeatures]{TranscriptDb} object. Possible locations include
  coding, intron, 3UTR, 5UTR, intergenic or unknown. Intergenic
  variants are those that do not match (i.e., overlap with) a transcript. Currently all 
  \code{findOverlaps} operations in \code{locateVariants} are performed with
  \code{type} = `within'. This requires that the variant fall
   completely within one of the defined regions or else it is classified as 
   `unknown'.
}


\value{
  A \link[IRanges]{DataFrame} is returned with columns of \code{queryHits}, \code{txID}, 
  \code{geneID} and \code{Location}. \code{queryHits} is a map back to the original variants
  from the \code{query}. The \code{txID} and \code{geneID} values come from the
  \link[GenomicFeatures]{TranscriptDb} object. For intergenic variants, the \code{geneID} 
  contains identifiers for genes on either side of the variant (preceding and
  following). Each row in the result represents a transcript hit by a variant. When a variant hits 
  multiple transcripts there will be multiple rows for that variant.
}

\author{Valerie Obenchain}

\seealso{
  \code{\link{predictCoding}}
}

\examples{
  #library(TxDb.Hsapiens.UCSC.hg18.knownGene)
  #library(BSgenome.Hsapiens.UCSC.hg18)
  #data(variants)
 
  #txdb <- Hsapiens_UCSC_hg18_knownGene_TxDb
  #loc <- locateVariants(variants, txdb)
}

\keyword{methods}
