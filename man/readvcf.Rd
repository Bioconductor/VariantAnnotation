\name{readVcf}
\alias{readVcf}

\title{Read VCF (variant call) files into a SummarizedExperiment object}

\description{Read VCF files into a SummarizedExperiment object}

\usage{
  readVcf(file, index=file, param=ScanBcfParam, raw=FALSE, ...)
}

\arguments{
  \item{file}{The character() file names(s) of the VCF file to be processed.
  }
  \item{index}{The character() file names of the index file of the VCF file
    being processed; this is given \emph{without} the `.tbi' extension. 
    The default value is the name of the VCF file. 
  }
  \item{param}{An instance of \code{\linkS4class{ScanBcfParam}}. This 
    influences what fields and which records are imported. If a 
    \code{param} is provided, an index for the VCF file must exist.
  }
  \item{raw}{A logical indicating if the ranges in \code{rowData} 
    should appear exactly as represented in the VCF file or if
    the ranges should reflect only the basepairs involved in
    the variation.

    In VCF format, the value in the POS field refers to the position of the 
    first base in the REF field. For InDels, the REF string includes the 
    base before the event which is reflected in the POS field.

    As an example, we use the string ACTG which has positions 1,2,3,4.
    An deletion of the T at position 3 would be represented as
    POS = 2, REF = CT, ALT = "." in VCF format. If \code{raw}=TRUE
    the range for this variant would have start = end = 2.
    If \code{raw} = FALSE, the range for this variant would be
    start = end = 3.
  }
  \item{\dots}{Additional arguments, passed to methods.
  }
}

\details{
  The \code{readVcf} function parses VCF files into a
  \code{\linkS4class{SummarizedExperiment}}.  \code{readVcf}
  calls \code{Rsamtools} \code{link{scanBcf}} function
  the details of which can be found on its help page. 
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
vcfFile <- system.file("extdata", "ex1.vcf", package="VariantAnnotation") 
res <- readVcf(vcfFile, index=character(0))
head(rowData(res))

## reference and alternate alleles with ranges
values(rowData(res))[["ref"]][2:5]
values(rowData(res))[["alt"]][2:5]
}


\keyword{manip}
