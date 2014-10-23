\name{readVcf}
\alias{readVcf}
\alias{readVcf,character,ANY,ScanVcfParam-method}
\alias{readVcf,character,ANY,missing-method}
\alias{readVcf,character,missing,missing-method}
\alias{readVcf,TabixFile,ANY,ScanVcfParam-method}
\alias{readVcf,TabixFile,ANY,RangedData-method}
\alias{readVcf,TabixFile,ANY,RangesList-method}
\alias{readVcf,TabixFile,ANY,GRanges-method}
\alias{readVcf,TabixFile,ANY,GRangesList-method}
\alias{readVcf,TabixFile,ANY,missing-method}

%% lightweight read* functions
\alias{readInfo}
\alias{readGeno}
\alias{readGT}


\title{Read VCF files}

\description{Read Variant Call Format (VCF) files}

\usage{

  \S4method{readVcf}{TabixFile,ANY,ScanVcfParam}(file, genome, param, 
      ..., row.names=TRUE)
  \S4method{readVcf}{TabixFile,ANY,RangedData}(file, genome, param, 
      ..., row.names=TRUE)
  \S4method{readVcf}{TabixFile,ANY,RangesList}(file, genome, param, 
      ..., row.names=TRUE)
  \S4method{readVcf}{TabixFile,ANY,GRanges}(file, genome, param, 
      ..., row.names=TRUE)
  \S4method{readVcf}{TabixFile,ANY,GRangesList}(file, genome, param, 
      ..., row.names=TRUE)
  \S4method{readVcf}{TabixFile,ANY,missing}(file, genome, param, 
      ..., row.names=TRUE)
  \S4method{readVcf}{character,ANY,ScanVcfParam}(file, genome, param, 
      ..., row.names=TRUE)
  \S4method{readVcf}{character,ANY,missing}(file, genome, param, 
      ..., row.names=TRUE)
  \S4method{readVcf}{character,missing,missing}(file, genome, param, 
      ..., row.names=TRUE)

## Lightweight functions to read a single variable
readInfo(file, x, param=ScanVcfParam(), ..., row.names=TRUE)
readGeno(file, x, param=ScanVcfParam(), ..., row.names=TRUE)
readGT(file, nucleotides=FALSE, param=ScanVcfParam(), ..., row.names=TRUE)
}

\arguments{
  \item{file}{A \code{\link{TabixFile}} instance or character() name of the VCF 
    file to be processed. When ranges are specified in \code{param},
    \code{file} must be a \code{\link{TabixFile}}. 

    Use of the \code{\link{TabixFile}} methods are encouraged as they are
    more efficient than the character() methods. See ?\code{TabixFile} and
    ?\code{indexTabix} for help creating a \code{\link{TabixFile}}.
  }
  \item{genome}{A \code{character} or \code{Seqinfo} object.
    \itemize{
      \item{\code{character}:}{ Genome identifier as a single string or named 
        character vector. Names of the character vector correspond to 
        chromosome names in the file. This identifier replaces the genome 
        information in the VCF \code{Seqinfo} (i.e., \code{seqinfo(vcf)}).
      }
      \item{\code{Seqinfo}:}{ When \code{genome} is provided as a \code{Seqinfo}
        it is propagated to the VCF output. If seqinfo information can be 
        obtained from the file, 
        (i.e., seqinfo(scanVcfHeader(fl)) is not empty), the output
        \code{Seqinfo} is a product of merging the two.

        If a param (i.e., ScanVcfParam) is used in the call to \code{readVcf}, 
        the seqlevels  of the param ranges must be present in \code{genome}.
      }
    }
  }
  \item{param}{A instance of \code{\linkS4class{ScanVcfParam}}, \code{GRanges},
    \code{GRangesList}, \code{RangedData} or \code{RangesList}. VCF files can be 
    subset on genomic coordinates (ranges) or elements in the VCF fields. Both 
    genomic coordinates and VCF elements can be specified in a 
    \code{\linkS4class{ScanVcfParam}}. See ?\code{ScanVcfParam} for details.
  }
  \item{x}{\code{character} name of single \code{info} or \code{geno} field to 
    import. Applicable to \code{readInfo} and \code{readGeno} only.
  }
  \item{row.names}{A \code{logical} specifying if rownames should be returned. 
    In the case of \code{readVcf}, rownames appear on the \code{GRanges} in the
    \code{rowData} slot.
  }
  \item{nucleotides}{A \code{logical} indicating if genotypes should be returned
    as nucleotides instead of the numeric representation.
    Applicable to \code{readGT} only.
  }
  \item{\dots}{Additional arguments, passed to methods.
  }
}

\details{
  \describe{
    \item{Data Import: }{
      \describe{
        \item{VCF object: }{
          \code{readVcf} imports records from bzip compressed or uncompressed 
          VCF files. Data are parsed into a \code{\linkS4class{VCF}} object 
          using the file header information if available. To import a subset 
          of ranges the VCF must have a Tabix index file. An index file can be 
          created with \code{bzip} and \code{indexTabix} functions.

          The \code{readInfo}, \code{readGeno} and \code{readGT} functions
          are lightweight versions of \code{readVcf} and import a single
          variable. The return object is a vector, matrix or CompressedList
          instead of a VCF class. 
          }
      } 
      \code{readVcf} calls \code{\link{scanVcf}}, the details of which can be 
      found with \code{?scanVcf}.
    }
    \item{Data type: }{
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
      information in the file header as follows:

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
    \item{Missing data: }{
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
  \code{readVcf} returns a \code{\linkS4class{VCF}} object. See ?\code{VCF} for 
  complete details of the class structure. \code{readGT}, \code{readInfo} and
  \code{readGeno} return a \code{matrix}.

  \describe{
    \item{rowData: }{
      The CHROM, POS, ID and REF fields are used to create a \code{GRanges}
      object. The ranges are created using POS as the start value and width of 
      the reference allele (REF). The IDs become the rownames. If they are 
      missing (i.e., \sQuote{.}) a string of CHROM:POS_REF/ALT is used instead. 
      The \code{genome} argument is stored in the seqinfo of the \code{GRanges} 
      and can be accessed with \code{genome(<VCF>)}.

      One elementMetadata column, \code{paramRangeID}, is included with the
      \code{rowData}. This ID is meaningful when multiple ranges are specified 
      in the \code{ScanVcfParam} and distinguishes which records match
      each range.
    }
    \item{fixed: }{
      REF, ALT, QUAL and FILTER fields of the VCF are parsed into a 
      \code{DataFrame}.
    }
    \item{info: }{
      Data from the INFO field of the VCF is parsed into a \code{DataFrame}.
    }
    \item{geno: }{
      If present, the genotype data are parsed into a list of \code{matrices} 
      or \code{arrays}. Each list element represents a field in the FORMAT 
      column of the VCF file. Rows are the variants, columns are the samples. 
    }
    \item{colData: }{
      This slot contains a \code{DataFrame} describing the samples. If present, 
      the sample names following FORMAT in the VCF file become the row names.
    }
    \item{exptData: }{ 
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
  \code{\link{expand,CollapsedVCF-method}}
}

\examples{
  fl <- system.file("extdata", "ex2.vcf", package="VariantAnnotation") 
  vcf <- readVcf(fl, "hg19")
  ## vcf <- readVcf(fl, c("20"="hg19"))  ## 'genome' as named vector

  ## ---------------------------------------------------------------------
  ## Header and genome information 
  ## ---------------------------------------------------------------------
  vcf

  ## all header information
  hdr <- header(vcf)

  ## header information for 'info' and 'fixed' fields
  info(hdr)
  fixed(hdr)

  ## ---------------------------------------------------------------------
  ## Accessors
  ## ---------------------------------------------------------------------
  ## fixed fields together
  head(fixed(vcf), 5)

  ## fixed fields separately 
  filt(vcf)
  ref(vcf) 

  ## info data 
  info(hdr)
  info(vcf)
  info(vcf)$DP

  ## geno data 
  geno(hdr)
  geno(vcf)
  head(geno(vcf)$GT)

  ## genome
  unique(genome(rowData(vcf)))

  ## ---------------------------------------------------------------------
  ## Data subsets with lightweight read* functions 
  ## ---------------------------------------------------------------------

  ## Import a single 'info' or 'geno' variable
  DP <- readInfo(fl, "DP")
  HQ <- readGeno(fl, "HQ")

  ## Import GT as numeric representation 
  GT <- readGT(fl)
  ## Import GT as nucleotides 
  GT <- readGT(fl, nucleotides=TRUE)

  ## ---------------------------------------------------------------------
  ## Data subsets with ScanVcfParam
  ## ---------------------------------------------------------------------

  ## Subset on genome coordinates:
  ## 'file' must have a Tabix index
  rngs <- GRanges("20", IRanges(c(14370, 1110000), c(17330, 1234600)))
  names(rngs) <- c("geneA", "geneB")
  param <- ScanVcfParam(which=rngs) 
  compressVcf <- bgzip(fl, tempfile())
  idx <- indexTabix(compressVcf, "vcf")
  tab <- TabixFile(compressVcf, idx)
  vcf <- readVcf(tab, "hg19", param)

  ## When data are subset by range ('which' argument in ScanVcfParam),
  ## the 'paramRangeID' column provides a map back to the original 
  ## range in 'param'.
  rowData(vcf)[,"paramRangeID"]
  vcfWhich(param)

  ## Subset on samples:
  ## Consult the header for the sample names.
  samples(hdr) 
  ## Specify one or more names in 'samples' in a ScanVcfParam.
  param <- ScanVcfParam(samples="NA00002")
  vcf <- readVcf(tab, "hg19", param)
  geno(vcf)$GT

  ## Subset on 'fixed', 'info' or 'geno' fields:
  param <- ScanVcfParam(fixed="ALT", geno=c("GT", "HQ"), info=c("NS", "AF"))
  vcf_tab <- readVcf(tab, "hg19", param)
  info(vcf_tab)
  geno(vcf_tab)

  ## No ranges are specified in the 'param' so tabix file is not
  ## required. Instead, the uncompressed VCF can be used as 'file'.
  vcf_fname <- readVcf(fl, "hg19", param)

  ## The header will always contain information for all variables
  ## in the original file reguardless of how the data were subset.
  ## For example, all 'geno' fields are listed in the header 
  geno(header(vcf_fname))

  ## but only 'GT' and 'HQ' are present in the VCF object.
  geno(vcf_fname)

  ## Subset on both genome coordinates and 'info', 'geno' fields: 
  param <- ScanVcfParam(geno="HQ", info="AF", which=rngs)
  vcf <- readVcf(tab, "hg19", param)

  ## When any of 'fixed', 'info' or 'geno' are omitted (i.e., no
  ## elements specified) all records are retrieved. Use NA to indicate
  ## that no records should be retrieved. This param specifies
  ## all 'fixed fields, the "GT" 'geno' field and none of 'info'.
  ScanVcfParam(geno="GT", info=NA)

  ## ---------------------------------------------------------------------
  ## Iterate through VCF with 'yieldSize' 
  ## ---------------------------------------------------------------------
  fl <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
  param <- ScanVcfParam(fixed="ALT", geno=c("GT", "GL"), info=c("LDAF"))
  tab <- TabixFile(fl, yieldSize=4000)
  open(tab)
  while (nrow(vcf_yield <- readVcf(tab, "hg19", param=param)))
      cat("vcf dim:", dim(vcf), "\n")
  close(tab)

}

\keyword{manip}
