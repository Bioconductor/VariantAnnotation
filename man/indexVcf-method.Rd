\name{indexVcf}
\alias{indexVcf}
\alias{indexVcf,character-method}
\alias{indexVcf,VcfFile-method}
\alias{indexVcf,VcfFileList-method}

\title{Create Index files for VCF file}

\description{Create Index file for VCF if it does not exist}

\usage{
\S4method{indexVcf}{character}(x, ...)
\S4method{indexVcf}{VcfFile}(x, ...)
\S4method{indexVcf}{VcfFileList}(x, ...)
}

\arguments{
  \item{x}{Either character(), \code{VcfFile}, or \code{VcfFileList}}
  \item{...}{Additional arguments to \code{\link[Rsamtools]{indexTabix}}}
}

\details{
  If a character vector is given, assumes they are the path(s) to the
  VCF file. If the index file for that VCF file does not exist, one is
  created. A \code{VcfFile} or \code{VcfFileList} is returned.

  If a \code{VcfFile} or \code{VcfFileList} is given, the index file is
  checked, if it does not exist it will crete one. If the index file was
  NA or missing, the path of the associated VCF file is used as the
  index file path. An updated \code{VcfFile} or \code{VcfFileList} is
  returned. 
}

\value{VcfFile or VcfFileList}

\author{Lori Shepherd}

\seealso{
  \code{\link{VcfFile}}
}

\examples{
fl <- system.file("extdata", "chr7-sub.vcf.gz", package="VariantAnnotation",
                  mustWork=TRUE)
vcf1 <- indexVcf(fl)
vcf1
}

\keyword{manip}
