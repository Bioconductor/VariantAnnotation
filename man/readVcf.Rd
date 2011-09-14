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
\S4method{readVcf}{character,ANY}(file, ..., param, vcfRanges=FALSE)
\S4method{readVcf}{character,missing}(file, ..., param, vcfRanges=FALSE)
\S4method{readVcf}{TabixFile,GRanges}(file, ..., param, vcfRanges=FALSE)
\S4method{readVcf}{TabixFile,RangedData}(file, ..., param, vcfRanges=FALSE)
\S4method{readVcf}{TabixFile,RangesList}(file, ..., param, vcfRanges=FALSE)
\S4method{readVcf}{TabixFile,ScanVcfParam}(file, ..., param, vcfRanges=FALSE)
\S4method{readVcf}{TabixFile,missing}(file, ..., param, vcfRanges=FALSE)
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
  \item{vcfRanges}{A logical indicating if ranges and alleles should appear 
    exactly as in VCF file or if they should be modified for compatibility 
    with the Bioconductor standards for representing ranges. 

    In a VCF file, the value in the POS field refers to the position of the 
    first base in the REF field. For indels, the REF string includes the 
    base before the event which is reflected in the POS field.
    Reguardless of the value of `vcfRanges', the ranges in the \code{rowData} will 
    be of width 1 and refect the exact position in the VCF POS field.

    If \code{vcfRanges}=TRUE, the only ranges available will be those in the 
    \code{GRanges} object in the \code{rowData}. If \code{vcfRanges}=FALSE
    modified ranges will be supplied in the alt metadata field. These
    ranges will reflect the \code{Bioconductor} standard of insertions
    having a width of length 0 and the range reflecting only those 
    basepairs involved in the substitution, insertion or deletion.

    We use the string ACTG with positions 1,2,3,4 as an example.
    In a VCF file a deletion of the T at position 3 would be represented as
    POS = 2, REF = CT, ALT = ".". If \code{vcfRanges}=TRUE
    the only ranges present for this variant will be in the \code{GRanges}
    instance in the \code{rowData} and will have start = end = 2.
    If \code{vcfRanges} = FALSE, the range in the \code{rowData} will be the
    same but an additional range of start = end = 3 will be supplied in 
    the \code{alt} metatdata slot of the \code{rowData}. The alleles
    will also be modified such that the reference allele is T instead
    of CT. The alternate allele will be the same missing value of ".".
  }
  \item{\dots}{Additional arguments, passed to methods.
  }
}


\details{
  \code{readVcf} imports records from a VCF file into a 
  \code{\linkS4class{SummarizedExperiment}} object. 
  The \code{param} argument can be used to select a subset of
  records when reading from compressed (bgzip) VCF files.

  The function calls \code{Rsamtools} \code{link{scanVcf}}, the 
  details of which can be found with \code{?scanVcf}. 
}

\value{
  The object returned is a \code{\linkS4class{SummarizedExperiment}}.
  See the help page for \code{\linkS4class{SummarizedExperiment}} for
  complete details of the class structure.

  The `rowData' slot contains a \code{\linkS4class{GRanges}} instance 
  with metadata columns of `pos', `qual', `filter', `info', `ref' and `alt'. 
  These correspond to the `POS' (reference position), `QUAL' (phred-scaled
  quality for assertion made in ALT), `FILTER' (pass if position passed
  all filters), `INFO' (additional information), `REF' (reference base(s)),
  `ALT' (alternate non-reference alleles called on at least one of the
  samples) fields in the VCF file. The row names of the 
  \code{\linkS4class{GRanges}} object are the `ID' values
  from the VCF file if they exist. See references for complete details
  of the VCF file format.

  The `assays' slot is a list the same length as the number of fields
  in the `FORMAT' column of the VCF file. Each list element is a matrix
  of data from the `GENO' field where there is a row for each variant
  and a column for each sample.

  The `colData' slot contains a \code{DataFrame} describing the samples.
  If present, the sample names (i.e., the names following `GENO' in the VCF file), 
  become the row names.
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

## These ranges and alleles have been adjusted to reflect the exact variant
## and are compatible with the Bioconductor standards for representing ranges. 
vcf_adj <- readVcf(tab, param=rngs)
values(rowData(vcf_adj))[["REF"]]
values(rowData(vcf_adj))[["ALT"]]

## Here vcfRanges=TRUE and the alleles and ranges appear exactly 
## as represented in the VCF file. Care must be taken when
## using these ranges for matching or overlap methods as the
## basepairs and range positions reflect more than the exact variant.
vcf_noadj <- readVcf(tab, param=rngs, vcfRanges=TRUE)
values(rowData(vcf_noadj))[["REF"]]
values(rowData(vcf_noadj))[["ALT"]]

}


\keyword{manip}
