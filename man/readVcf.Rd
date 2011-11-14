\name{readVcf}
\alias{readVcf}
\alias{readVcf,character,missing-method}
\alias{readVcf,character,ANY-method}
\alias{readVcf,TabixFile,GRanges-method}
\alias{readVcf,TabixFile,RangedData-method}
\alias{readVcf,TabixFile,RangesList-method}
\alias{readVcf,TabixFile,ScanVcfParam-method}
\alias{readVcf,TabixFile,missing-method}


\title{Read VCF (variant call format) files into a SummarizedExperiment object}

\description{Read VCF files into a SummarizedExperiment object}

\usage{
\S4method{readVcf}{character,ANY}(file, ..., param)
\S4method{readVcf}{character,missing}(file, ..., param)
\S4method{readVcf}{TabixFile,GRanges}(file, ..., param)
\S4method{readVcf}{TabixFile,RangedData}(file, ..., param)
\S4method{readVcf}{TabixFile,RangesList}(file, ..., param)
\S4method{readVcf}{TabixFile,ScanVcfParam}(file, ..., param)
\S4method{readVcf}{TabixFile,missing}(file, ..., param)
}

\arguments{
  \item{file}{The character() file name(s) of the VCF file to be
    processed, or an instance of class \code{\link{TabixFile}}.
  }
  \item{param}{A instance of \code{GRanges}, \code{RangedData},
    \code{RangesList}, or \code{\link[Rsamtools]{ScanVcfParam}} 
    providing the sequence names and regions to be
    parsed.
  }
  \item{\dots}{Additional arguments, passed to methods.
  }
}


\details{
  \code{readVcf} imports records from a VCF file into a 
  \code{\linkS4class{SummarizedExperiment}} object. 
  The \code{param} argument can be used to select a subset of
  records when reading from compressed (bgzip) VCF files.

  The function calls \code{Rsamtools} \code{\link{scanVcf}}, the 
  details of which can be found with \code{?scanVcf}. 
}

\value{
  The object returned is a \code{\linkS4class{SummarizedExperiment}}.
  See the help page for \code{\linkS4class{SummarizedExperiment}} for
  complete details of the class structure.

  The \code{rowData} slot contains a \code{\linkS4class{GRanges}} instance 
  with metadata columns of QUAL (phred-scaled quality for the assertion 
  made in ALT), FILTER (pass if position passed all filters), INFO 
  (additional information), REF (reference base(s)), and ALT 
  (alternate non-reference alleles called on at least one of the
  samples) fields in the VCF file. The row names of the 
  \code{\linkS4class{GRanges}} object are the ID values
  from the VCF file if they exist. See references for complete details
  of the VCF file format. The ranges in the \code{rowData} reflect the
  range of the REF alleles.

  The \code{assays} slot is a list the same length as the number of fields
  in the FORMAT column of the VCF file. Each list element is a matrix
  of data from the GENO field where there is a row for each variant
  and a column for each sample.

  The \code{colData} slot contains a \code{DataFrame} describing the samples.
  If present, the sample names (i.e., the names following GENO in the VCF file), 
  become the row names.

  The \code{exptData} slot contains a \code{SimpleList} with the VCF file header
  information. 
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

\author{Valerie Obenchain <vobencha@fhcrc.org>}

\seealso{
  \code{\link{indexTabix}}
  \code{\link{TabixFile}}
  \code{\link{scanTabix}}
  \code{\link{scanBcf}}
}

\examples{

## read all records of VCF into a SummarizedExperiment
vcfFile <- system.file("extdata", "ex1.vcf", package="VariantAnnotation") 
vcf <- readVcf(vcfFile)
names(assays(vcf))

## compress on the fly and read a subset of data
rngs <- GRanges(seqnames="1", ranges=IRanges(start=0, end=100000))
compressVcf <- bgzip(vcfFile, tempfile())
idx <- indexTabix(compressVcf, "vcf")
tab <- TabixFile(compressVcf, idx)
vcf <- readVcf(tab, param=rngs)

## vcf file header information
exptData(vcf)$HEADER
}


\keyword{manip}
