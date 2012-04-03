\name{readVcf}
\alias{readVcf}
\alias{readVcf,character,character,ScanVcfParam-method}
\alias{readVcf,character,character,missing-method}
\alias{readVcf,character,missing,missing-method}
\alias{readVcf,TabixFile,character,ScanVcfParam-method}
\alias{readVcf,TabixFile,character,RangedData-method}
\alias{readVcf,TabixFile,character,RangesList-method}
\alias{readVcf,TabixFile,character,GRanges-method}
\alias{readVcf,TabixFile,character,missing-method}


\title{Read VCF files}

\description{Read Variant Call Format (VCF) files}

\usage{
  \S4method{readVcf}{TabixFile,character,ScanVcfParam}(file, genome, param, ...)
  \S4method{readVcf}{TabixFile,character,RangedData}(file, genome, param, ...)
  \S4method{readVcf}{TabixFile,character,RangesList}(file, genome, param, ...)
  \S4method{readVcf}{TabixFile,character,GRanges}(file, genome, param, ...)
  \S4method{readVcf}{TabixFile,character,missing}(file, genome, param, ...)
  \S4method{readVcf}{character,character,ScanVcfParam}(file, genome, param, ...)
  \S4method{readVcf}{character,character,missing}(file, genome, param, ...)
  \S4method{readVcf}{character,missing,missing}(file, genome, param, ...)
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
    \item{Data Import :}{
      \describe{
        \item{VCF object :}{
          \code{readVcf} imports records from bzip compressed or uncompressed 
          VCF files. Data are parsed into a \code{\linkS4class{VCF}} object 
          using the file header information if available. To import a subset 
          of ranges the VCF must have a Tabix index file. An index file can be 
          created with \code{bzip} and \code{indexTabix} functions. See
          examples. 
          }
      } 
      \code{readVcf} calls \code{\link{scanVcf}}, the details of which can be 
      found with \code{?scanVcf}.
    }
    \item{Data type :}{
      CHROM, POS, ID and REF fields are used to create the \code{GRanges}
      stored in the \code{rowData} slot of the \code{VCF} object. Access
      with \code{rowData} accessor.

      REF, ALT, QUAL and FILTER are parsed into the \code{DataFrame} in the 
      \code{fixed} slot. Because ALT can have more than one value per variant 
      it is represented as a \code{DNAStringSetList}. REF is a \code{DNAStringSet},
      QUAL is \code{numeric} and FILTER is a \code{character}. Accessors include
      \code{fixed}, \code{ref}, \code{alt}, \code{qual}, and \code{filt}.

      Data from the INFO field can be accessed with the \code{info} accessor.
      Genotype data (i.e., data immediately following the FORMAT field in the 
      VCF) can be accessed with the \code{geno} accessor. INFO and genotype data 
      types are determined according to the \sQuote{Number} and \sQuote{Type} 
      information in the file header as follows :

      If \sQuote{Number} is 1, \sQuote{info} data are parsed into a 
      \code{vector}. \sQuote{geno} data are parsed into a \code{matrix} where 
      the columns are the samples.

      If \sQuote{Number} is an integer >1, \sQuote{info} data are parsed into a
      \code{DataFrame} with the indicated number of columns. \sQuote{geno} are
      parsed into an \code{array} with the same dimentions as \sQuote{Number}. 
      Columns of the \sQuote{geno} matrices are the samples.

      If \sQuote{Number} is \sQuote{.}, \sQuote{A} or \sQuote{G}, a \code{matrix} 
      is used for both \sQuote{info} and \sQuote{geno} data.

      When the VCF header does not contain data type information, the data are 
      returned as a single unparsed column named \sQuote{INFO} or \sQuote{GENO}.
    }
    \item{Missing data  :}{
      Missing data in VCF files are represented by a dot ("."). \code{readVcf} 
      retains the dot as a character string for data type character and converts 
      it to \code{NA} for data types numeric or double. 

      Because the data are stored in rectangluar data structures there is a
      value for each \code{info} and \code{geno} field element in the \code{VCF} 
      class. If the element was missing or was not collected for a particular 
      variant the value will be \code{NA}.
    }
  }
}

\value{
  The object returned is a \code{\linkS4class{VCF}} class instance.
  See ?\code{VCF} for complete details of the class structure.

  \describe{
    \item{rowData:}{
      The CHROM, POS, ID and REF fields are used to create a \code{GRanges}
      object. The ranges are created using POS as the start value and width of 
      the reference allele (REF). The IDs become the rownames. If they are 
      missing (i.e., \sQuote{.}) a colon separated string of chromosome and 
      start position will be used instead. The \code{genome} argument is stored 
      in the seqinfo of the \code{GRanges} and can be accessed with 
      \code{genome(<VCF>)}.

      One elementMetadata column, \code{paramRangeID}, is included with the
      \code{rowData}. This ID is meaningful when multiple ranges are specified 
      in the \code{ScanVcfParam}. This ID distinguishes which records match
      each range.
    }
    \item{fixed:}{
      REF, ALT, QUAL and FILTER fields of the VCF are parsed into a 
      \code{DataFrame}.
    }
    \item{info:}{
      Data from the INFO field of the VCF is parsed into a \code{DataFrame}.
    }
    \item{geno:}{
      If present, the genotype data are parsed into a list of \code{matrices} 
      or \code{arrays}. Each list element represents a field in the FORMAT 
      column of the VCF file. Rows are the variants, columns are the samples. 
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
  \code{\link{indexTabix}},
  \code{\link{TabixFile}},
  \code{\link{scanTabix}},
  \code{\link{scanBcf}},
  \code{readVcfLongForm}
}

\examples{
  fl <- system.file("extdata", "ex2.vcf", package="VariantAnnotation") 
  vcf <- readVcf(fl, "hg19")

  ## ---------------------------------
  ## Header and genome information 
  ## ---------------------------------
  vcf

  ## extract the VCFHeader object from exptData
  hdr <- exptData(vcf)[["header"]]
  ## use accessors to extract the data
  info(hdr)
  fixed(hdr)

  ## genome information is stored in the GRanges
  unique(genome(rowData(vcf)))

  ## ---------------------------------
  ## Accessors
  ## ---------------------------------
  ## fixed fields together
  head(fixed(vcf), 5)

  ## fixed fields separately 
  filt(vcf)
  ref(vcf) 

  ## info data 
  info(hdr)
  info(vcf)
  values(info(vcf))[["DP"]]

  ## geno data 
  geno(hdr)
  geno(vcf)
  head(geno(vcf)$GT)

  identical(geno(vcf)[["GQ"]], geno(vcf)$GQ)

  ## ---------------------------------
  ## Data subsets 
  ## ---------------------------------

  ## Subset on genome coordinates :
  ## 'file' must have a Tabix index
  rngs <- GRanges("20", IRanges(c(14370, 1110000), c(17330, 1234600)))
  names(rngs) <- c("geneA", "geneB")
  param <- ScanVcfParam(which=rngs) 
  compressVcf <- bgzip(fl, tempfile())
  idx <- indexTabix(compressVcf, "vcf")
  tab <- TabixFile(compressVcf, idx)
  vcf <- readVcf(tab, "hg19", param)

  ## 'paramRangeID' associates the records to the ranges in the param
  rowData(vcf)

  ## Subset on 'fixed', 'info' or 'geno' VCF elements :
  param <- ScanVcfParam(fixed="ALT", geno=c("GT", "HQ"), info=c("NS", "AF"))
  vcf_tab <- readVcf(tab, "hg19", param)
  vcf_fname <- readVcf(fl, "hg19", param)
  identical(vcf_tab, vcf_fname) 

  ## Subset on both genome coordinates and VCF elements : 
  param <- ScanVcfParam(geno="HQ", info="AF", which=rngs)
  vcf <- readVcf(tab, "hg19", param)

  ## When any of 'fixed', 'info' or 'geno' are omitted (i.e., no
  ## elements specified) all records are retrieved. Use NA to indicate
  ## that no records should be retrieved. This param specifies
  ## all 'fixed fields, the "GT" 'geno' field and none of 'info'.
  ScanVcfParam(geno="GT", info=NA)
}

\keyword{manip}
