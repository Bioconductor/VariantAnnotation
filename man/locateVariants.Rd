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
  The \code{query} object is returned with two columns added to the metadata.
  The 'location' column describes the variant location as being in an exon, 
  intron, intergenic, 3UTR, or 5UTR. The 'genes' column contains the gene ID 
  from the organism package specified. If the variant is intergenic, 
  the 'genes' column contains the names of the genes on either side.
}

\references{
%% ~put references to the literature/web site here ~
}

\author{Valerie Obenchain}

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
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
