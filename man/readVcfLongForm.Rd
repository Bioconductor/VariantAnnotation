\name{readVcfLongForm}
\alias{readVcfLongForm}
\alias{readVcfLongForm,character,character,ScanVcfParam-method}
\alias{readVcfLongForm,character,character,missing-method}
\alias{readVcfLongForm,character,missing,missing-method}
\alias{readVcfLongForm,TabixFile,character,ScanVcfParam-method}
\alias{readVcfLongForm,TabixFile,character,RangedData-method}
\alias{readVcfLongForm,TabixFile,character,RangesList-method}
\alias{readVcfLongForm,TabixFile,character,GRanges-method}
\alias{readVcfLongForm,TabixFile,character,missing-method}


\title{Read VCF files into a long form GRanges}

\description{Read Variant Call Format (VCF) files into a long form GRanges}

\usage{
  \S4method{readVcfLongForm}{TabixFile,character,ScanVcfParam}(file, genome, param, ...)
  \S4method{readVcfLongForm}{TabixFile,character,RangedData}(file, genome, param, ...)
  \S4method{readVcfLongForm}{TabixFile,character,RangesList}(file, genome, param, ...)
  \S4method{readVcfLongForm}{TabixFile,character,GRanges}(file, genome, param, ...)
  \S4method{readVcfLongForm}{TabixFile,character,missing}(file, genome, param, ...)
  \S4method{readVcfLongForm}{character,character,ScanVcfParam}(file, genome, param, ...)
  \S4method{readVcfLongForm}{character,character,missing}(file, genome, param, ...)
  \S4method{readVcfLongForm}{character,missing,missing}(file, genome, param, ...)
}

\arguments{
  \item{file}{A \code{\link{TabixFile}} instance or character() name of the VCF 
    file to be processed. When ranges are specified in \code{param},
    \code{file} must be a \code{\link{TabixFile}}. 

    Use of the \code{\link{TabixFile}} methods are encouraged as they are
    more efficient than the character() methods. See ?\code{TabixFile} and
    ?\code{indexTabix} for help creating a \code{\link{TabixFile}}.
  }
  \item{genome}{The character() name of the genome the variants are mapped 
    to. This information is stored in the genome slot of the \code{seqinfo}
    associated with the \code{\linkS4class{GRanges}} object in \code{rowData}.
  }
  \item{param}{A instance of \code{\linkS4class{ScanVcfParam}}, \code{GRanges},
    \code{RangedData} or \code{RangesList}. VCF files can be subset on genomic
    coordinates (ranges) or elements in the VCF fields. Both genomic coordinates
    and VCF elements can be specified in a \code{\linkS4class{ScanVcfParam}}.
    See ?\code{ScanVcfParam} for details.
  }
  \item{\dots}{Additional arguments, passed to methods.
  }
}

\details{
  \describe{
     \item{Long Form GRanges object :}{
       \code{readVcfLongForm} reads data from a VCF file in the same manner 
       as \code{readVcf}. Input arguments and the ability to subset the
       data on ranges or VCF fields are the same. The return object is a 
       long form \code{GRanges} expanded to the length of the \sQuote{unlisted}
       alternate allele (\code{ALT}). The \code{fixed} and \code{info} data
       are also expanded as elementMetadata columns. Currently no \code{geno}
       information is included.
     }
  }
}

\value{
  A \code{GRanges} of variant locations. When the alternate allele column
  (\code{ALT}) column in the VCF file has a single value per variant the 
  \code{GRanges} will have a single row per variant (i.e., not expanded).
  When variants have multiple alternate allele values the \code{GRanges}
  will have repeated rows for the variants with multiple values.
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

\author{
  Valerie Obenchain <vobencha@fhcrc.org>
}

\seealso{
  \code{\link{readVcf}}
}

\examples{
  ## All data
  fl <- system.file("extdata", "ex2.vcf", package="VariantAnnotation")
  gr1 <- readVcfLongForm(fl, "hg19")

  ## Subset on ranges :
  rngs <- GRanges("20", IRanges(c(14370, 1234567), c(17330, 1234567)))
  names(rngs) <- c("geneA", "geneB")
  param <- ScanVcfParam(which=rngs) 
  compressVcf <- bgzip(fl, tempfile())
  idx <- indexTabix(compressVcf, "vcf")
  tab <- TabixFile(compressVcf, idx)
  gr2 <- readVcfLongForm(tab, "hg19", param)

  ## 'paramRangeID' associates the records to the ranges in the param
  gr2 

  ## Subset on 'fixed' or 'info' VCF elements :
  param <- ScanVcfParam(which=rngs, fixed="ALT", info=c("NS", "AF"))
  gr3 <- readVcfLongForm(tab, "hg19", param)
}

\keyword{manip}
