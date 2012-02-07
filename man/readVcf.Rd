\name{readVcf}
\alias{readVcf}
\alias{readVcf,character,character,ScanVcfParam-method}
\alias{readVcf,character,character,RangedData-method}
\alias{readVcf,character,character,RangesList-method}
\alias{readVcf,character,character,GRanges-method}
\alias{readVcf,character,character,missing-method}
\alias{readVcf,character,missing,ANY-method}
\alias{readVcf,TabixFile,character,ScanVcfParam-method}
\alias{readVcf,TabixFile,character,RangedData-method}
\alias{readVcf,TabixFile,character,RangesList-method}
\alias{readVcf,TabixFile,character,GRanges-method}
\alias{readVcf,TabixFile,character,missing-method}
\alias{readVcf,TabixFile,missing,ANY-method}


\title{Read VCF files}

\description{Read Variant Call Format (VCF) files}

\usage{
  \S4method{readVcf}{character,character,ScanVcfParam}(file, genome, param, ...)
  \S4method{readVcf}{character,character,RangedData}(file, genome, param, ...)
  \S4method{readVcf}{character,character,RangesList}(file, genome, param, ...)
  \S4method{readVcf}{character,character,GRanges}(file, genome, param, ...)
  \S4method{readVcf}{character,character,missing}(file, genome, param, ...)
  \S4method{readVcf}{TabixFile,character,ScanVcfParam}(file, genome, param, ...)
  \S4method{readVcf}{TabixFile,character,RangedData}(file, genome, param, ...)
  \S4method{readVcf}{TabixFile,character,RangesList}(file, genome, param, ...)
  \S4method{readVcf}{TabixFile,character,GRanges}(file, genome, param, ...)
  \S4method{readVcf}{TabixFile,character,missing}(file, genome, param, ...)
}

\arguments{
  \item{file}{The character() file name of the VCF file to be
    processed, a \code{GRanges}, \code{RangedData}, \code{RangesList} or an 
    instance of class \code{\link{TabixFile}}.
  }
  \item{genome}{The character() name of the genome the variants are mapped 
    to. This information is stored in the genome slot of the \code{seqinfo}
    associated with the \code{\linkS4class{GRanges}} object in \code{rowData}.
  }
  \item{param}{A instance of \code{\linkS4class{ScanVcfParam}} 
    providing the desired ranges, info fields or geno fields to be parsed.
  }
  \item{\dots}{Additional arguments, passed to methods.
  }
}

\details{
  \describe{
    \item{Data Import}{
      \code{readVcf} imports records from bzip compressed or uncompressed VCF 
      files. Data are parsed into a \code{\linkS4class{VCF}} 
      object using the VCF file header information if available. To import a
      subset of records the VCF file must have a Tabix index file. An index
      file can be created with \code{bzip} and \code{indexTabix} functions.
      See examples below.

      The function calls \code{Rsamtools} \code{\link{scanVcf}}, the details of 
      which can be found with \code{?scanVcf}.
    }
    \item{Data Representation}{
      Data from the INFO and FORMAT (GENO) fields in the VCF file are parsed
      into a \code{vector}, \code{matrix}, \code{arrays} or 
      \code{CompressedList} based on the 'Number' field in the header.

      If 'Number' is 1, 'info' data are parsed into a \code{vector}. 'geno' 
      data are parsed into a \code{matrix} where the columns are the samples.

      If 'Number' is an integer >1, \sQuote{info} data are parsed into a
      \code{DataFrame} with the same number of columns. \sQuote{geno} data are
      parsed into an \code{array} with the same dimentions as 'Number'. The 
      columns of the \sQuote{geno} matrices are the samples.

      If 'Number' is '.', 'A' or 'G', a \code{CompressedList} is used for both
      'info' and 'geno' data.

      When the VCF header does not contain data type information for the
      sQuote{INFO} or \sQuote{GENO} elements the data are returned as a single
      unparsed column named as \sQuote{INFO} or \sQuote{GENO}.
    }
    \item{Flat Format}{
      When the \code{asGRanges} argument to \code{\linkS4class{ScanVcfParam}} 
      is specified as \sQuote{info} or \sQuote{geno}, the return value is a 
      'long form' \code{GRanges} object instead of a \code{VCF} object. In 
      this form, each row represents a variant with the \sQuote{info} or
      \sQuote{geno} data expanded in the metadata columns.
    }
  }
}

\value{
  The object returned is a \code{\linkS4class{VCF}}.
  See ?\code{VCF} for complete details of the class structure.

  \describe{
    \item{rowData:}{
      The first 7 required fields of the VCF file are parsed into a 
      \code{GRanges}. CHROM, POS and ID fields are the seqnames, 
      ranges and rownames. If the ID is missing (i.e., '.') a colon separated
      string of chromosome and start position will be used as the rowname.
      REF, ALT, QUAL and FILTER are included as elementMetadata.

      The ranges are created using the position(POS) and width of the reference 
      allele (REF). The \code{genome} argument is stored in the seqinfo 
      of the \code{GRanges} and can be accessed with 
      \code{genome(rowData(readVcf(file)))}.
    }
    \item{info:}{
      Data from the INFO field of the VCF is parsed into a \code{DataFrame}.
      Individual elements can be a \code{vector}, \code{CompressedList},
      or an \code{array}.
    }
    \item{geno:}{
      If present, genomic data from the samples specified in the FORMAT field
      are parsed into a list of \code{matrices} or \code{arrays} with a matrix 
      for each field in the FORMAT column. The rows are the variants, the 
      columns the samples. 
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

\author{
  Valerie Obenchain <vobencha@fhcrc.org>
}

\seealso{
  \code{\link{indexTabix}}
  \code{\link{TabixFile}}
  \code{\link{scanTabix}}
  \code{\link{scanBcf}}
}

\examples{
  fl <- system.file("extdata", "ex1.vcf", package="VariantAnnotation") 
  vcf <- readVcf(fl, "hg19")

  ## ---------------------------------
  ## Header and genome information 
  ## ---------------------------------
  vcf
  exptData(vcf)[["HEADER"]]
  unique(genome(rowData(vcf)))

  ## ---------------------------------
  ## Accessors
  ## ---------------------------------
  ## Fixed fields together
  head(fixed(vcf), 5)

  ## Fixed fields separately 
  filt(vcf)
  ref(vcf) 

  ## info data 
  exptData(vcf)[["HEADER"]][["INFO"]]
  info(vcf)
  values(info(vcf))[["DP"]]

  ## geno data 
  exptData(vcf)[["HEADER"]][["FORMAT"]]
  geno(vcf)
  head(geno(vcf)$GT)

  identical(geno(vcf)[["GQ"]], geno(vcf)$GQ)

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
  vcf <- readVcf(tab, "hg19", param)

  ## Subset on 'info' or 'geno' elements :
  ## An index file is not required if no ranges (genome coordinates) 
  ## are specified.
  param <- ScanVcfParam(geno=c("GT", "HQ"), info=c("NS", "AF"))
  vcf_tab <- readVcf(tab, "hg19", param)
  vcf_fname <- readVcf(fl, "hg19", param)
  identical(vcf_tab, vcf_fname) 

  ## Subsetting on both genome coordinates and 'info' or 'geno' elements : 
  param <- ScanVcfParam(geno="HQ", info="AF", which=rngs)
  vcf <- readVcf(tab, "hg19", param)

  ## When either of 'fixed', 'info' or 'geno' are omitted (i.e., no
  ## elements specified) all records are retrieved. Use NA to indicate
  ## that no records should be retrieved. This param object specifies
  ## the "GT" 'geno' field and none of 'info' fields.
  ScanVcfParam(geno="GT", info=NA)

  ## ---------------------------------
  ## Long form GRanges : 
  ## ---------------------------------
  ##
  ## -  under construction - 

  ## A 'long form' GRanges has been expanded to match the expanded
  ## length of elements in the 'info' or 'geno' fields that have more
  ## than one value per record. This results in a GRanges with
  ## potentially multiple rows per record but each row represents a
  ## unique combination of the elementMetadata columns. The
  ## elementMetadata columns are the values supplied in 'info' or 'geno'
  ## in the ScanVcfParam object.
  ## To specify a long form GRanges, set asGRanges=TRUE and specify one
  ## or more elements in the 'info' filed OR specify one element in the
  ## 'geno' field.  In this usage only one of 'info' or 'geno' can be
  ## specified and 'info' can accept multiple elements but 'geno' can
  ## take only one. See ?ScanVcfParam for more details.
  ##
  ## First retrieve available INFO and FORMAT fields from the file
  fl <- system.file("extdata", "ex2.vcf", package="VariantAnnotation")
  hdr <- scanVcfHeader(fl)
  hdr[[1]]$Header$INFO
  hdr[[1]]$Header$FORMAT

  ## Specify one or more fields to 'info'. The "AF" element has two
  ## values for variant rs6040355. In the result there are now two rows
  ## for this variant with the same ranges but unique elementMetadata.

  param <- ScanVcfParam(asGRanges=TRUE, info=c("NS", "AF", "H2"))
  #readVcf(fl, "hg19", param)

  ## When creating a long form GRanges only one 'geno' element can be
  ## specified at a time.
  param <- ScanVcfParam(asGRanges=TRUE, geno="GT")
  #readVcf(fl, "hg19", param)
}

\keyword{manip}
