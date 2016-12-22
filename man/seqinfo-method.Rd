\name{seqinfo}
\alias{seqinfo}
\alias{seqinfo,VcfFile-method}
\alias{seqinfo,VcfFileList-method}

\title{Get seqinfo for VCF file}

\description{Get seqinfo for VCF file}

\usage{
\S4method{seqinfo}{VcfFile}(x)
\S4method{seqinfo}{VcfFileList}(x)
}

\arguments{
  \item{x}{Either character(), \code{VcfFile}, or \code{VcfFileList}}
}

\details{

  If a \code{VcfFile}The file header is scanned an appropriate seqinfo
  object in given. If a \code{VcfFileList} is given, all file headers are
  scanned, and appropriate combined seqinfo object is given.
}

\value{Seqinfo object}

\author{Lori Shepherd}

\seealso{
  \code{\link{VcfFile}}, \code{\link[GenomeInfoDb]{Seqinfo}}
}

\examples{
fl <- system.file("extdata", "chr7-sub.vcf.gz", package="VariantAnnotation",
                  mustWork=TRUE)
vcf <- VcfFile(fl)
seqinfo(vcf)
}

\keyword{manip}
