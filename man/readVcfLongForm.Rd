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
  fl <- system.file("extdata", "ex1.vcf", package="VariantAnnotation") 
  vcf <- readVcfLongForm(fl, "hg19")

  ## 

  ## ---------------------------------
  ## Reading in subsets
  ## ---------------------------------
  fl <- system.file("extdata", "ex2.vcf", package="VariantAnnotation")

  ## Subset on genome coordinates :
  ## Range subsets are specified in a GRanges, RangedData or RangesList
  ## and supplied as the 'which' argument to ScanVcfParam().
  ## When range subsets are specified the 'file' must have a Tabix index.
  rngs <- GRanges("20", IRanges(1110000, 1234600))
  param <- ScanVcfParam(which=rngs) 
  compressVcf <- bgzip(fl, tempfile())
  idx <- indexTabix(compressVcf, "vcf")
  tab <- TabixFile(compressVcf, idx)
  vcf <- readVcfLongForm(tab, "hg19", param)

  ## Subset on 'fixed', 'info' or 'geno' VCF elements :
  ## An index file is not required if genome coordinates (ranges)
  ## are not specified.
  param <- ScanVcfParam(fixed="ALT", geno=c("GT", "HQ"), info=c("NS", "AF"))
  vcf_tab <- readVcfLongForm(tab, "hg19", param)
  vcf_fname <- readVcfLongForm(fl, "hg19", param)
  identical(vcf_tab, vcf_fname) 

  ## Subsetting on both genome coordinates and VCF fields : 
  param <- ScanVcfParam(geno="HQ", info="AF", which=rngs)
  vcf <- readVcfLongForm(tab, "hg19", param)

  ## When any of 'fixed', 'info' or 'geno' are omitted (i.e., no
  ## elements specified) all records are retrieved. Use NA to indicate
  ## that no records should be retrieved. This param specifies
  ## all 'fixed fields, the "GT" 'geno' field and none of 'info'.
  ScanVcfParam(geno="GT", info=NA)

  ## ---------------------------------
  ## Long form GRanges : 
  ## ---------------------------------
  ## A 'long form' GRanges is one that has been expanded to accomodate
  ## the 'unlisted' length of the 'fixed', 'info' or 'geno' fields
  ## requested. This results in a GRanges with multiple rows per record 
  ## where each row represents a unique combination of the elementMetadata 
  ## columns. 

  ## To specify a long form GRanges, set asGRanges=TRUE in the ScanVcfParam
  ## object and specify 'fixed', 'info' and/or 'geno' elements.
  ## Possible values for 'fixed' are ALT, QUAL or FILTER. To see available
  ## values for INFO and FORMAT use scanVcfHeader. 
  fl <- system.file("extdata", "ex2.vcf", package="VariantAnnotation")
  hdr <- scanVcfHeader(fl)
  info(hdr)
  geno(hdr) 

  ## Request all 'fixed' fields, "GT" from 'geno' and "H2" from 'info'.
  ## See ?ScanVcfParam for more details.
  param <- ScanVcfParam(asGRanges=TRUE, info="H2", geno="GT")
  vcf <- readVcfLongForm(fl, "hg19", param)
  head(vcf)
}

\keyword{manip}
