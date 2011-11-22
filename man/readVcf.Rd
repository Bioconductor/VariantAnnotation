\name{readVcf}
\alias{readVcf}
\alias{readVcf,character,character,ANY-method}
\alias{readVcf,character,character,missing-method}
\alias{readVcf,TabixFile,character,GRanges-method}
\alias{readVcf,TabixFile,character,RangedData-method}
\alias{readVcf,TabixFile,character,RangesList-method}
\alias{readVcf,TabixFile,character,ScanVcfParam-method}
\alias{readVcf,TabixFile,character,missing-method}


\title{Read VCF (variant call format) files into a SummarizedExperiment object}

\description{Read VCF files into a SummarizedExperiment object}

\usage{
  \S4method{readVcf}{character,character,ANY}(file, genome, ..., param)
  \S4method{readVcf}{character,character,missing}(file, genome, ..., param)
  \S4method{readVcf}{TabixFile,character,GRanges}(file, genome, ..., param)
  \S4method{readVcf}{TabixFile,character,RangedData}(file, genome, ..., param)
  \S4method{readVcf}{TabixFile,character,RangesList}(file, genome, ..., param)
  \S4method{readVcf}{TabixFile,character,ScanVcfParam}(file, genome, ..., param)
  \S4method{readVcf}{TabixFile,character,missing}(file, genome, ..., param)
}

\arguments{
  \item{file}{The character() file name of the VCF file to be
    processed, or an instance of class \code{\link{TabixFile}}.
  }
  \item{genome}{The character() name of the genome the variants are mapped 
    to. This information is stored in the genome slot of the \code{seqinfo}
    associated with the \code{\linkS4class{GRanges}} object in \code{rowData}.
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
  \code{readVcf} imports records from bzip compressed or uncompressed VCF 
  files. Data are parsed into a \code{\linkS4class{SummarizedExperiment}} 
  object using the VCF file header information if available. To import a
  subset of records the VCF file must have a Tabix index file. An index
  file can be created with \code{bzip} and \code{indexTabix} functions.
  See examples below.

  The function calls \code{Rsamtools} \code{\link{scanVcf}}, and 
  \code{\link{unpackVcf}} the details of which can be found with 
  \code{?scanVcf} and \code{?unpackVcf}.
}

\value{
  The object returned is a \code{\linkS4class{SummarizedExperiment}}.
  See ?\code{SummarizedExperiment} for complete details of the class structure.

  \describe{
    \item{rowData:}{
      The 8 required fields of the VCF file are parsed into a \code{GRanges} in 
      the \code{rowData} slot. CHROM, POS and ID fields are in the seqnames, 
      ranges and rownames. QUAL, FILTER, INFO, REF and ALT are included as 
      metadata columns. The ranges are created using the position(POS) and width 
      of the reference allele (REF). Genome information from the \code{genome}
      argument is stored in the seqinfo associated with the \code{GRanges} and
      can be accessed with \code{genome(rowData(readVcf(file)))}.
    }
    \item{assays:}{
      If present, genomic data from the samples specified in the FORMAT field
      are parsed into the \code{assays} slot. This is a list of matrices with
      one matrix for each of fields in the FORMAT column. The rows are the
      variants, the columns the samples. 
    }
    \item{colData:}{
      This slot contains a \code{DataFrame} describing the samples. If present, 
      the sample names following FORMAT in the VCF file become the row names.
    }
    \item{exptData:}{ 
      Header information present in the file is put into a \code{SimpleList}
      in \code{exptData}.
    }
  }
  See references for complete details of the VCF file format. 
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
  vcfFile <- system.file("extdata", "ex1.vcf", package="VariantAnnotation") 
  vcf <- readVcf(vcfFile, "hg19")
  names(assays(vcf))

  ## Explore the header information 
  vcf
  exptData(vcf)[["HEADER"]]
  exptData(vcf)[["HEADER"]][["META"]]

  ## rowData contains the 8 required VCF fields
  exptData(vcf)[["HEADER"]][["FILTER"]]
  exptData(vcf)[["HEADER"]][["INFO"]]
  unique(genome(rowData(vcf)))
  head(rowData(vcf), 5)

  ## Genotype data described in FORMAT are parsed into the assays slot 
  exptData(vcf)[["HEADER"]][["FORMAT"]]
  assays(vcf)
  head(assays(vcf)$GT)

  ## Reading in a subset of ranges :
  ## To create a Tabix index the file must be compressed
  compressVcf <- bgzip(vcfFile, tempfile())
  idx <- indexTabix(compressVcf, "vcf")
  tab <- TabixFile(compressVcf, idx)
  ## Define a subset of ranges and read them in
  rngs <- GRanges(seqnames="16", ranges=IRanges(start=0, end=100000))
  vcf <- readVcf(tab, "hg19", param=rngs)
}

\keyword{manip}
