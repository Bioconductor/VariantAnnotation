\name{annotateVariants}
\alias{annotateVariants}

\title{Annotate variants}

\description{
  Annotate variants with respect to their functional consequence on genes
}

\usage{
  annotateVariants(query, subject, orgPackage, BSgenomeOrganism, 
                   variantColumn, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{query}{
%%     ~~Describe \code{x} here~~
}
  \item{subject}{
%%     ~~Describe \code{cds_by_tx} here~~
}
  \item{orgPackage}{
%%     ~~Describe \code{tx_seqs} here~~
}
  \item{BSgenomeOrganism}{
%%     ~~Describe \code{tx_seqs} here~~
}
  \item{variantColumn}{
%%     ~~Describe \code{variant_column} here~~
}
  \item{\dots}{
%%     ~~Describe \code{variant_column} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}


\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
  \dontrun{
  ## using sample data from ANNOVAR
  library(TxDb.Hsapiens.UCSC.hg18.knownGene)
  library(BSgenome.Hsapiens.UCSC.hg18)
  library(org.Hs.eg.db)
  data(variants)

  txdb <- Hsapiens_UCSC_hg18_knownGene_TxDb
  strand(variants) <- "+"
  res <- annotateVariants(variants, txdb, orgPackage="org.Hs.eg.db", 
      BSgenomeOrganism=Hsapiens, variantColumn="obsAllele")
  }
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
