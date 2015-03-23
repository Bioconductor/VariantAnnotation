\name{scanVcf}
\Rdversion{1.1}

\alias{scanVcfHeader}
\alias{scanVcfHeader,character-method}
\alias{scanVcf}
\alias{scanVcf,character,ScanVcfParam-method}
\alias{scanVcf,character,missing-method}
\alias{scanVcf,connection,missing-method}

\alias{scanVcfHeader,TabixFile-method}
\alias{scanVcf,TabixFile,GRanges-method}
\alias{scanVcf,TabixFile,RangedData-method}
\alias{scanVcf,TabixFile,RangesList-method}
\alias{scanVcf,TabixFile,ScanVcfParam-method}
\alias{scanVcf,TabixFile,missing-method}

\alias{scanVcfHeader,TabixFile-method}
\alias{scanVcf,TabixFile,GRanges-method}
\alias{scanVcf,TabixFile,RangedData-method}
\alias{scanVcf,TabixFile,RangesList-method}
\alias{scanVcf,TabixFile,ScanVcfParam-method}
\alias{scanVcf,TabixFile,missing-method}


\title{
  Import VCF files
}

\description{
  Import Variant Call Format (VCF) files in text or binary format
}

\usage{
scanVcfHeader(file, ...)
\S4method{scanVcfHeader}{character}(file, ...)

scanVcf(file, ..., param)
\S4method{scanVcf}{character,ScanVcfParam}(file, ..., param)
\S4method{scanVcf}{character,missing}(file, ..., param)
\S4method{scanVcf}{connection,missing}(file, ..., param)

\S4method{scanVcfHeader}{TabixFile}(file, ...)
\S4method{scanVcf}{TabixFile,missing}(file, ..., param)
\S4method{scanVcf}{TabixFile,ScanVcfParam}(file, ..., param)
\S4method{scanVcf}{TabixFile,GRanges}(file, ..., param)
\S4method{scanVcf}{TabixFile,RangedData}(file, ..., param)
\S4method{scanVcf}{TabixFile,RangesList}(file, ..., param)

\S4method{scanVcfHeader}{TabixFile}(file, ...)
\S4method{scanVcf}{TabixFile,missing}(file, ..., param)
\S4method{scanVcf}{TabixFile,ScanVcfParam}(file, ..., param)
\S4method{scanVcf}{TabixFile,GRanges}(file, ..., param)
\S4method{scanVcf}{TabixFile,RangedData}(file, ..., param)
\S4method{scanVcf}{TabixFile,RangesList}(file, ..., param)
}

\arguments{
  \item{file}{For \code{scanVcf} and \code{scanVcfHeader}, the character() 
    file name, \code{\link{TabixFile}}, or class \code{connection} 
    (\code{file()} or \code{bgzip()}) of the \sQuote{VCF} file to be 
    processed.
  }
  \item{param}{A instance of \code{\linkS4class{ScanVcfParam}} influencing 
    which records are parsed and the \sQuote{INFO} and \sQuote{GENO} information
    returned.
  }
  \item{...}{Additional arguments for methods
  }
}

\details{
  The argument \code{param} allows portions of the file to be input, but 
  requires that the file be bgzip'd and indexed as a 
  \code{\linkS4class{TabixFile}}.

  \code{scanVcf} with \code{param="missing"} and \code{file="character"}
  or \code{file="connection"} scan the entire file. With
  \code{file="connection"}, an argument \code{n} indicates the number of
  lines of the VCF file to input; a connection open at the beginning of
  the call is open and incremented by \code{n} lines at the end of the
  call, providing a convenient way to stream through large VCF files.

  The INFO field of the scanned VCF file is returned as a single
  \sQuote{packed} vector, as in the VCF file. The GENO field is a list of 
  matricies, each matrix corresponds to a field as defined in the FORMAT 
  field of the VCF header. Each matrix has as many rows as
  scanned in the VCF file, and as many columns as there are samples. As
  with the INFO field, the elements of the matrix are
  \sQuote{packed}. The reason that INFO and GENO are returned packed is
  to facilitate manipulation, e.g., selecting particular rows or
  samples in a consistent manner across elements.
}

\value{
  \code{scanVcfHeader} returns a \code{VCFHeader} object with
  header information parsed into five categories, \code{samples}, 
  \code{meta}, \code{fixed}, \code{info} and \code{geno}. Each
  can be accessed with a `getter' of the same name
  (e.g., info(<VCFHeader>)).

  \code{scanVcf} returns a list, with one element per range. Each list
  has 7 elements, obtained from the columns of the VCF specification:
    \describe{
      \item{rowRanges}{
        \code{GRanges} instance derived from \code{CHROM}, \code{POS}, 
        \code{ID}, and the width of \code{REF}
      }
      \item{REF}{
        reference allele
      }
      \item{ALT}{
        alternate allele
      }
      \item{QUAL}{
        phred-scaled quality score for the assertion made in ALT
      }
      \item{FILTER}{
        indicator of whether or not the position passed all filters applied
      }
      \item{INFO}{
        additional information
      }
      \item{GENO}{
        genotype information immediately following the FORMAT field in the VCF
      }
    }
  The \code{GENO} element is itself a list, with elements corresponding
  to those defined in the VCF file header. For \code{scanVcf}, elements
  of GENO are returned as a matrix of records x samples; if the
  description of the element in the file header indicated multiplicity
  other than 1 (e.g., variable number for \dQuote{A}, \dQuote{G}, or
  \dQuote{.}), then each entry in the matrix is a character string with
  sub-entries comma-delimited.
}

\references{
  \url{http://vcftools.sourceforge.net/specs.html} outlines the VCF
  specification.

  \url{http://samtools.sourceforge.net/mpileup.shtml} contains
  information on the portion of the specification implemented by
  \code{bcftools}.

  \url{http://samtools.sourceforge.net/} provides information on
  \code{samtools}.
}

\seealso{
  \code{\link{readVcf}} 
  \code{\link{BcfFile}} 
  \code{\link{TabixFile}}
}

\author{
  Martin Morgan and Valerie Obenchain <vobencha@fhcrc.org>
}

\examples{
  fl <- system.file("extdata", "ex2.vcf", package="VariantAnnotation")
  scanVcfHeader(fl)
  vcf <- scanVcf(fl)
  ## value: list-of-lists
  str(vcf)
  names(vcf[[1]][["GENO"]])
  vcf[[1]][["GENO"]][["GT"]]
}

\keyword{ manip }
