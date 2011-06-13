\name{locateVariants}
\alias{locateVariants}

\title{Locate variants}

\description{Identify variant location with respect to genes}

\usage{
  locateVariants(query, subject, ...)
}

\arguments{
  \item{query}{A \link[GenomicRanges]{GRanges} object}
  \item{subject}{A \link[GenomicFeatures]{TranscriptDb} object}
  \item{\dots}{Additional arguments}
}

\details{
%%  ~~ If necessary, more details than the description above ~~
}

\value{
  The \code{query} object is returned with two columns added to the metadata.
  The 'location' column describes the variant location as exonic, intronic,
  intergenic, 3UTR, or 5UTR. The 'genes' column contains the gene symbol
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
  ## using sample data from ANNOVAR
  library(TxDb.Hsapiens.UCSC.hg18.knownGene)
  library(BSgenome.Hsapiens.UCSC.hg18)
  data(variants)
 
  txdb <- Hsapiens_UCSC_hg18_knownGene_TxDb
  location <- locateVariants(variants, txdb)
  }
}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
