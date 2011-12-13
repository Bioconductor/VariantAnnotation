\name{MatrixToSnpMatrix}
\alias{MatrixToSnpMatrix}
\alias{MatrixToSnpMatrix,matrix,GRanges-method}


\title{Convert genotype calls from a VCF file to a SnpMatrix object}

\description{
  Convert a matrix of genotype calls from the "GT" FORMAT field of a VCF
  file to a \linkS4class{SnpMatrix}.
}

\usage{MatrixToSnpMatrix(callMatrix, mapGRanges, ...)}

\arguments{
  \item{callMatrix}{A \code{matrix} of genotype data from the "GT" FORMAT field 
    of a VCF file.  This \code{matrix} is located in the \code{assays} slot of 
    the \code{SummarizedExperiment} produced from a call to \code{readVcf}.
  }
  \item{mapGRanges}{An optional \code{GRanges} with REF and ALT allele information
    in the elementMetadata. This \code{GRanges} is located in the \code{rowData} 
    slot of the \code{SummarizedExperiment} produced by a call to \code{readVcf}. 
  }
  \item{\dots}{Additional arguments, passed to methods.
  }
}

\details{
  \code{MatrixToSnpMatrix} converts a matrix of genotype calls from the 
  "GT" FORMAT field of a VCF file into a \linkS4class{SnpMatrix}. The
  following caveats apply, 

  - no distinction is made between phased and unphased genotypes

  - only diploid calls are included; others are set to NA

  - only single nucleotide variants are included; others are set to NA

  - variants with >1 ALT allele are set to NA

  In VCF files, 0 represents the reference allele and integers greater than 0
  represent the alternate alleles (i.e., 2, 3, 4 would indicate the 2nd, 3rd
  or 4th allele in the ALT field for a particular variant). This function only 
  supports variants with a single alternate allele and therefore the alternate 
  values will always be 1. Genotypes are stored in the \linkS4class{SnpMatrix} 
  as 0, 1, 2 or 3 where 0 = missing, 1 = "0/0", 2 = "0/1" or "1/0" and 
  3 = "1/1". In \linkS4class{SnpMatrix} terminology, "A" is the reference allele 
  and "B" is the risk allele. Equivalent statements to those made with 0 and 1 
  allele values would be 0 = missing, 1 = "A/A", 2 = "A/B" or "B/A" and 3 = "B/B".
}

\value{
  A list with the following elements,

  \item{genotypes}{The output genotype data as an object of class 
    \code{"SnpMatrix"}. The columns are snps and the rows are the samples. See 
    the help page for \linkS4class{SnpMatrix} for complete details of the class 
    structure.}

  \item{map}{A \code{DataFrame} giving the snp names and alleles at each locus.
    The \code{setNA} column indicates which variants were set to NA because they 
    met one or more of the caveats stated above.} 
}

\author{Valerie Obenchain <vobencha@fhcrc.org>}

\seealso{
  \linkS4class{SnpMatrix}
}

\examples{
  ## Read a vcf file into a SummarizedExperiment object
  vcfFile <- system.file("extdata", "ex2.vcf", package="VariantAnnotation") 
  se <- readVcf(vcfFile, "hg19")
 
  ## Genotype calls are stored in the assays slot and the
  ## GRanges with the map information is in the rowData slot.
  assays(se)
  calls <- assays(se)$GT
  mapgr <- rowData(se)
  MatrixToSnpMatrix(calls, mapgr)
}

\keyword{manip}
