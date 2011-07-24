\name{readVcf}
\alias{readVcf}
\alias{readVcf,character,missing-method}
\alias{readVcf,character,ANY-method}
\alias{readVcf,TabixFile,GRanges-method}
\alias{readVcf,TabixFile,RangedData-method}
\alias{readVcf,TabixFile,RangesList-method}


\title{Read VCF (variant call format) files into a SummarizedExperiment object}

\description{Read VCF files into a SummarizedExperiment object}

\usage{
\S4method{readVcf}{character,ANY}(file, index=paste(file, "tbi", sep = "."),
..., param=NULL, raw=FALSE)
\S4method{readVcf}{character,missing}(file, ..., param, raw=FALSE)
\S4method{readVcf}{TabixFile,GRanges}(file, ..., param, raw=FALSE)
\S4method{readVcf}{TabixFile,RangedData}(file, ..., param, raw=FALSE)
\S4method{readVcf}{TabixFile,RangesList}(file, ..., param, raw=FALSE)
}

\arguments{
  \item{file}{The character() file name(s) of the VCF file to be
    processed, or more flexibly an instance of class
    \code{\link{TabixFile}}.
  }
  \item{param}{A instance of \code{GRanges}, \code{RangedData},  or
    \code{RangesList} provide the sequence names and regions to be
    parsed.
  }
  \item{raw}{A logical indicating if modified ranges should appear in 
    the alt metadata field of \code{rowData}. 

    In a VCF file, the value in the POS field refers to the position of the 
    first base in the REF field. For indels, the REF string includes the 
    base before the event which is reflected in the POS field.
    Reguardless of the value of `raw', the ranges in the \code{rowData} will 
    be of width 1 and refect the exact position in the VCF POS field.

    If \code{raw}=TRUE, the only ranges available will be those in the 
    \code{GRanges} object in the \code{rowData}. If \code{raw}=FALSE
    modified ranges will be supplied in the alt metadata field. These
    ranges will reflect the \code{Bioconductor} standard of insertions
    having a width of length 0 and the range reflecting only those 
    basepairs involved in the substitution, insertion or deletion.

    We use the string ACTG with positions 1,2,3,4 as an example.
    In a VCF file a deletion of the T at position 3 would be represented as
    POS = 2, REF = CT, ALT = ".". If \code{raw}=TRUE
    the only ranges present for this variant will be in the \code{GRanges}
    instance in the \code{rowData} and will have start = end = 2.
    If \code{raw} = FALSE, the range in the \code{rowData} will be the
    same but an additional range of start = end = 3 will be supplied in 
    the \code{alt} metatdata slot of the \code{rowData}. The alleles
    will also be modified such that the reference allele is T instead
    of CT. The alternate allele will be the same missing value of ".".
  }
  \item{index}{The character() file names of the index file of the VCF file
    being processed; this is given \emph{without} the `.tbi' extension. 
    The default value is the name of the VCF file. 
  }
  \item{\dots}{Additional arguments, passed to methods.
  }
}


\details{
  The \code{readVcf} function parses VCF and compressed VCF files into a
  \code{\linkS4class{SummarizedExperiment}} object.  \code{readVcf}
  calls \code{Rsamtools} \code{link{scanBcf}} and \code{link{scanTabix}}
  functions, the details of which can be found on their respective help pages. 
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
  \code{\link{scanBcf}}
  \code{\link{scanBcfHeader}}
  \code{\link{ScanBcfParam}}
}

\examples{

## read entire VCF file
vcfFile <- system.file("extdata", "ex1.vcf", package="VariantAnnotation") 
vcf <- readVcf(vcfFile, index=character(0))
values(rowData(vcf))[["ref"]][2:5]
values(rowData(vcf))[["alt"]][2:5]


## read a portion of a compressed VCF file
gr <- GRanges(seqnames="1",
              ranges=IRanges(start=0, end=100000))
from <- vcfFile 
to <- tempfile()
compressVcf <- bgzip(from, to)
idx <- indexTabix(compressVcf, "vcf")
tab <- TabixFile(compressVcf, idx)
vcf_chr1 <- readVcf(tab, param=gr)
values(rowData(vcf_chr1))[["ref"]]
values(rowData(vcf_chr1))[["alt"]]


## return alleles and ranges exactly as represented in the VCF file
vcf_chr1_raw <- readVcf(tab, param=gr, raw=TRUE)
values(rowData(vcf_chr1_raw))[["ref"]]
values(rowData(vcf_chr1_raw))[["alt"]]

}


\keyword{manip}
