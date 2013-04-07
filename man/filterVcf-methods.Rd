\name{filterVcf}
\alias{filterVcf}
\alias{filterVcf,character-method}
\alias{filterVcf,TabixFile-method}
\alias{dbSNPFilter}
\alias{regionFilter}

\title{Filter VCF files}

\description{Filter Variant Call Format (VCF) files from one file to another}

\usage{
\S4method{filterVcf}{character}(file, genome, destination, ..., verbose = TRUE,
    index = FALSE, prefilters = FilterRules(), filters = FilterRules(),
    param = ScanVcfParam())

\S4method{filterVcf}{TabixFile}(file, genome, destination, ..., verbose = TRUE,
    index = FALSE, prefilters = FilterRules(), filters = FilterRules(),
    param = ScanVcfParam())
}

\arguments{

  \item{file}{A \code{character(1)} file path or \code{\link{TabixFile}}
    specifying the VCF file to be filtered.}

  \item{genome}{A \code{character(1)} identifier}

  \item{destination}{A \code{character(1)} path to the location where
    the filtered VCF file will be written.}

  \item{...}{Additional arguments, possibly used by future methods.}

  \item{verbose}{A \code{logical(1)} indicating whether progress
    messages should be printed.}

  \item{index}{A \code{logical(1)} indicating whether the filtered file
    should be compressed and indexed (using \code{\link{bgzip}} and
    \code{indexTabix}).}

  \item{prefilters}{A \code{\link{FilterRules}} instance contains rules for
    filtering un-parsed lines of the VCF file.}

  \item{filters}{A \code{\link{FilterRules}} instance contains rules for
    filtering fully parsed VCF objects.}

  \item{param}{A \code{\link{ScanVcfParam}} instance restricting input
    of particular \code{info} or \code{geno} fileds, or genomic
    locations.}

}

\details{

  This function transfers content of one VCF file to another, removing
  records that fail to satisfy \code{prefilters} and
  \code{filters}. Filtering is done in a memory efficient manner,
  iterating over the input VCF file in chunks of default size 100,000
  (when invoked with \code{character(1)} for \code{file}) or as
  specified by the \code{yieldSize} argument of \code{TabixFile} (when
  invoked with \code{TabixFile}).

  There are up to two passes. In the first pass, unparsed lines are
  passed to \code{prefilters} for filtering, e.g., searching for a fixed
  character string. In the second pass lines successfully passing
  \code{prefilters} are parsed into \code{VCF} instances and made
  available for further filtering. One or both of \code{prefilter} and
  \code{filter} can be present.

}

\value{The destination file path as a \code{character(1)}.}

\author{
  Martin Morgan \url{mailto:mtmorgan@fhcrc.org} and Paul Shannon
  \url{mailto:pshannon@fhcrc.org}.
}

\seealso{
  \code{\link{readVcf}}, \code{\link{writeVcf}}.
}

\examples{
fl <- system.file(package="VariantAnnotation", "extdata",
                  "chr22.vcf.gz")
destination <- tempfile()
pre <- FilterRules(list(isLowCoverageExomeSnp = function(x) {
    grepl("LOWCOV,EXOME", x, fixed=TRUE)
}))
filt <- FilterRules(list(isSNP = function(x) info(x)$VT == "SNP"))
filtered <- filterVcf(fl, "hg19", destination, prefilters=pre, filters=filt)
vcf <- readVcf(filtered, "hg19")
}

\keyword{manip}
