\name{MatrixToSnpMatrix}
\alias{MatrixToSnpMatrix}
\alias{MatrixToSnpMatrix,matrix,DNAStringSet,DNAStringSetList-method}


\title{Convert genotype calls from a VCF file to a SnpMatrix object}

\description{
  Convert a matrix of genotype calls from the "GT" FORMAT field of a VCF
  file to a \link[snpStats:SnpMatrix-class]{SnpMatrix}.
}

\usage{MatrixToSnpMatrix(callMatrix, ref, alt, ...)}

\arguments{
  \item{callMatrix}{A \code{matrix} of genotype data from the "GT" FORMAT field 
    of a VCF file.  This \code{matrix} is created with a call to \code{readVcf}
    and can be accessed with \code{geno<VCF>}.
  }
  \item{ref}{A \code{DNAStringSet} of reference alleles.
  }
  \item{alt}{A \code{DNAStringSetList} of alternate alleles. 
  }
  \item{\dots}{Additional arguments, passed to methods.
  }
}

\details{
  \code{MatrixToSnpMatrix} converts a matrix of genotype calls from the 
  "GT" FORMAT field of a VCF file into a 
  \link[snpStats:SnpMatrix-class]{SnpMatrix}. The following caveats apply, 

  - no distinction is made between phased and unphased genotypes

  - only diploid calls are included; others are set to NA

  - only single nucleotide variants are included; others are set to NA

  - variants with >1 ALT allele are set to NA

  In VCF files, 0 represents the reference allele and integers greater than 0
  represent the alternate alleles (i.e., 2, 3, 4 would indicate the 2nd, 3rd
  or 4th allele in the ALT field for a particular variant). This function only 
  supports variants with a single alternate allele and therefore the alternate 
  values will always be 1. Genotypes are stored in the SnpMatrix 
  as 0, 1, 2 or 3 where 0 = missing, 1 = "0/0", 2 = "0/1" or "1/0" and 
  3 = "1/1". In SnpMatrix terminology, "A" is the reference allele 
  and "B" is the risk allele. Equivalent statements to those made with 0 and 1 
  allele values would be 0 = missing, 1 = "A/A", 2 = "A/B" or "B/A" and 3 = "B/B".
}

\value{
  A list with the following elements,

  \item{genotypes}{The output genotype data as an object of class 
    \code{"SnpMatrix"}. The columns are snps and the rows are the samples. See 
    the help page for \link[snpStats:SnpMatrix-class]{SnpMatrix} for complete 
    details of the class structure.}

  \item{map}{A \code{DataFrame} giving the snp names and alleles at each locus.
    The \code{ignore} column indicates which variants were set to NA because they 
    met one or more of the caveats stated above.} 
}

\author{
  Valerie Obenchain <vobencha@fhcrc.org>
}

\seealso{
  \link{readVcf},
  \linkS4class{VCF},
  \link[snpStats:SnpMatrix-class]{SnpMatrix}
}

\examples{
  ## Read a vcf file into a VCF object.
  vcfFile <- system.file("extdata", "ex2.vcf", package="VariantAnnotation") 
  vcf <- readVcf(vcfFile, "hg19")
 
  ## Genotype calls are stored in the assays slot.
  ## Reference and alternate (variant) alleles are in ref and
  ## alt slots.
  geno(vcf)
  calls <- geno(vcf)$GT
  a0 <- values(ref(vcf))[["REF"]]
  a1 <- values(alt(vcf))[["ALT"]]
  mat <- MatrixToSnpMatrix(calls, a0, a1)

  ## The result is a list of length 2.
  names(mat)

  ## Compare original coding in the VCF file to the SnpMatrix coding.
  geno(vcf)$GT
  t(as(mat$genotype, "character"))

  ## The 'ignore' column of the map element indicates which variants 
  ## were set to missing as per the criteria stated in the details section. 
  mat$map

  ## Variant rs6040355 was ignored because it has more than one alternate
  ## allele. The microsat1 variant has a reference allele with length
  ## greater than 1 and more than 1 alternate allele.
  ## Variant chr20:1230237 was ignored because the alternate allele is not
  ## of length 1. 

  ## View the reference and alternate alleles. 
  fixed <- fixed(vcf)
  DataFrame(reference = values(fixed)[["REF"]], 
            alternate = CharacterList(as.list(values(fixed)[["ALT"]])),
            row.names = rownames(vcf))
}

\keyword{manip}
