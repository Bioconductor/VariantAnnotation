\name{locateVariants}
\alias{locateVariants}

\title{Locate variants}

\description{Identify variant location with respect to gene functionality}

\usage{
  locateVariants(query, subject, ...)
}

\arguments{
  \item{query}{A \link[GenomicRanges]{GRanges} object containing the variants
  }
  \item{subject}{A \link[GenomicFeatures]{TranscriptDb} object that will serve
    as the annotation reference
  }
  \item{\dots}{Additional arguments
  }
}

\details{
  Add some ..
}


\value{
  The \code{query} object is returned with three columns added to the metadata,
  Location, txID and geneID. The Location column describes the variant location 
  as one of exon, intron, intergenic, 3UTR, or 5UTR. The txID and geneID
  columns contain the transcript and gene IDs from the organism package specified. 
  If the variant is intergenic, the geneID column contains the IDs for genes
  on either side of the variant.
}

\author{Valerie Obenchain}

\seealso{
  \code{\link{predictCoding}}
}

\examples{
  \dontrun{
  library(TxDb.Hsapiens.UCSC.hg18.knownGene)
  library(BSgenome.Hsapiens.UCSC.hg18)
  data(variants)
 
  txdb <- Hsapiens_UCSC_hg18_knownGene_TxDb
  loc <- locateVariants(variants, txdb)
  }
}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
