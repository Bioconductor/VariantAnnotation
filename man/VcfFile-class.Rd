\name{VcfFile}
\Rdversion{1.1}
\docType{class}
\alias{VcfFile-class}
\alias{VcfFileList-class}
% con/destructors
\alias{VcfFile}
\alias{VcfFileList}

\title{Manipulate Variant Call Format (Vcf) files.}

\description{

  Use \code{VcfFile()} to create a reference to a Vcf file (and its
  index). Once opened, the reference remains open across calls to
  methods, avoiding costly index re-loading.

  \code{VcfFileList()} provides a convenient way of managing a list of
  \code{VcfFile} instances.

}

\usage{

## Constructors

VcfFile(file, index = paste(file, "tbi", sep="."), ...,
        yieldSize=NA_integer_)
VcfFileList(..., yieldSize=NA_integer_)

## Accessors
index(object)
path(object, ...)
isOpen(con, rw="")
yieldSize(object, ...)
yieldSize(object, ...) <- value
show(object)

## Opening / closing
open(con, ...)
close(con, ...)

}

\arguments{

  \item{con}{An instance of \code{VcfFile}.}

  \item{file}{A character(1) vector to the Vcf file
    path; can be remote (http://, ftp://).}

  \item{index}{A character(1) vector of the Vcf file index (.tbi file).}

  \item{yieldSize}{Number of records to yield each time the file is read
    from using \code{scanVcf} or \code{readVcf}.}

  \item{...}{Additional arguments. For \code{VcfFileList}, this can
    either be a single character vector of paths to Vcf files, or
    several instances of \code{VcfFile} objects.}

  \item{rw}{character() indicating mode of file.}

}

\section{Objects from the Class}{

  Objects are created by calls of the form \code{VcfFile()}.

}

\section{Fields}{

  \code{VcfFile} and \code{VcfFileList} classes inherit fields from the
  \code{\linkS4class{TabixFile}} and \code{\linkS4class{TabixFileList}} 
  classes.
}

\section{Functions and methods}{

  \code{VcfFile} and \code{VcfFileList} classes inherit methods from the
  \code{\linkS4class{TabixFile}} and \code{\linkS4class{TabixFileList}} 
  classes.

  Opening / closing:
  \describe{
    \item{open}{Opens the (local or remote) \code{path} and
      \code{index}. Returns a \code{VcfFile} instance. 
      \code{yieldSize} determines the number of records parsed during
      each call to \code{scanVcf} or \code{readVcf}; \code{NA} indicates 
      that all records are to be parsed.}

    \item{close}{Closes the \code{VcfFile} \code{con}; returning
      (invisibly) the updated \code{VcfFile}. The instance may be
      re-opened with \code{open.VcfFile}.}
  }

  Accessors: 
  \describe{
    \item{path}{Returns a character(1) vector of the Vcf path name.}

    \item{index}{Returns a character(1) vector of Vcf index (tabix file) name.}

    \item{yieldSize, yieldSize<-}{Return or set an integer(1) vector
      indicating yield size.}
  }
}
\author{Valerie Obenchain}

\examples{
fl <- system.file("extdata", "chr7-sub.vcf.gz", package="VariantAnnotation",
                  mustWork=TRUE)
vcffile <- VcfFile(fl)
vcffile

param <- GRanges("7", IRanges(c(55000000,  55900000), width=10000))
vcf <- readVcf(vcffile, "hg19", param=param)
dim(vcf)
}

\keyword{classes}
