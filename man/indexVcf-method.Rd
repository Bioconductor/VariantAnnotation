\name{indexVcf}
\alias{indexVcf}
\alias{indexVcf,character-method}
\alias{indexVcf,VcfFile-method}
\alias{indexVcf,VcfFileList-method}

\title{Create index files for VCF files}

\description{
  \code{indexVcf()} creates an index file for compressed VCF files, if
  the index file does not exist.
}

\usage{
\S4method{indexVcf}{character}(x, ...)
\S4method{indexVcf}{VcfFile}(x, ...)
\S4method{indexVcf}{VcfFileList}(x, ...)
}

\arguments{
  \item{x}{
    \code{character()}, \code{VcfFile}, or \code{VcfFileList} pointing
    to bgzf-compressed VCF files.
  }

  \item{...}{
    Additional arguments to \code{\link[Rsamtools]{indexTabix}}
  }
}

\details{
  If \code{x} is a character vector, assumes they are the path(s) to
  bgzf-compressed VCF file(s). If an index does not exist, one is
  created. VCF files can be compreseed using \code{\link{bgzip}}. A
  \code{VcfFile} or \code{VcfFileList} is returned.

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
fl <- system.file(
    "extdata", "chr7-sub.vcf.gz", package="VariantAnnotation",
    mustWork=TRUE
)
vcf1 <- indexVcf(fl)
vcf1
}

\keyword{manip}
