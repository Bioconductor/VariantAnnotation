\name{genotypeToSnpMatrix}

\alias{genotypeToSnpMatrix}
\alias{genotypeToSnpMatrix,CollapsedVCF-method}
\alias{genotypeToSnpMatrix,array-method}


\title{Convert genotype calls from a VCF file to a SnpMatrix object}

\description{
  Convert an array of genotype calls from the "GT", "GP", or "GL" 
  FORMAT field of a VCF file to a \link[snpStats:SnpMatrix-class]{SnpMatrix}.
}

\usage{
\S4method{genotypeToSnpMatrix}{CollapsedVCF}(x, uncertain=FALSE, ...)
\S4method{genotypeToSnpMatrix}{array}(x, ref, alt, ...)
}

\arguments{
  \item{x}{
    A \code{CollapsedVCF} object or a \code{array} of genotype data 
    from the "GT" FORMAT field of a VCF file.  This \code{array} is created 
    with a call to \code{readVcf} and can be accessed with \code{geno(<VCF>)}.
  }
  \item{uncertain}{
    A logical indicating whether the genotypes to convert
    should come from the "GT" field (\code{uncertain=FALSE}) or the "GP" or 
    "GL" field (\code{uncertain=TRUE}).
  }
  \item{ref}{
    A \code{DNAStringSet} of reference alleles.
  }
  \item{alt}{
    A \code{DNAStringSetList} of alternate alleles. 
  }
  \item{\dots}{
    Additional arguments, passed to methods.
  }
}

\details{
  \code{genotypeToSnpMatrix} converts an array of genotype calls from the 
  "GT", "GP", or "GL" FORMAT field of a VCF file into a 
  \link[snpStats:SnpMatrix-class]{SnpMatrix}. The following caveats apply, 
  \itemize{
    \item{no distinction is made between phased and unphased genotypes}
    \item{variants with >1 ALT allele are set to NA}
    \item{only single nucleotide variants are included; others are set to NA}
    \item{only diploid calls are included; others are set to NA}
  }

  In VCF files, 0 represents the reference allele and integers greater than 0
  represent the alternate alleles (i.e., 2, 3, 4 would indicate the 2nd, 3rd
  or 4th allele in the ALT field for a particular variant). This function only 
  supports variants with a single alternate allele and therefore the alternate 
  values will always be 1. Genotypes are stored in the SnpMatrix 
  as 0, 1, 2 or 3 where 0 = missing, 1 = "0/0", 2 = "0/1" or "1/0" and 
  3 = "1/1". In SnpMatrix terminology, "A" is the reference allele 
  and "B" is the risk allele. Equivalent statements to those made with 0 and 1 
  allele values would be 0 = missing, 1 = "A/A", 2 = "A/B" or "B/A" and
  3 = "B/B".

  If \code{uncertain=TRUE}, the posterior probabilities of the three
  genotypes ("A/A", "A/B", "B/B") are encoded (approximately) as byte
  values.  See the \link[snpStats]{snpStats} documentation for more
  details.
}

\value{
  A list with the following elements,

  \item{genotypes}{
    The output genotype data as an object of class 
    \code{"SnpMatrix"}. The columns are snps and the rows are the samples. 
    See ?\code{SnpMatrix} details of the class structure.
  }
  \item{map}{ 
    A \code{DataFrame} giving the snp names and alleles at each locus.
    The \code{ignore} column indicates which variants were set to \code{NA}
    (see \code{NA} criteria in 'details' section).
  } 
}

\author{
  Stephanie Gogarten, Valerie Obenchain <vobencha@fhcrc.org>
}

\seealso{
  \link{readVcf},
  \linkS4class{VCF},
  \link[snpStats:SnpMatrix-class]{SnpMatrix}
}

\examples{
  ## ----------------------------------------------------------------
  ## Non-probability based snp encoding using "GT"
  ## ----------------------------------------------------------------
  fl <- system.file("extdata", "ex2.vcf", package="VariantAnnotation") 
  vcf <- readVcf(fl, "hg19")

  ## This file has no "GL" or "GP" field so we use "GT".
  geno(vcf)

  ## Convert the "GT" FORMAT field to a SnpMatrix.
  mat <- genotypeToSnpMatrix(vcf)

  ## The result is a list of length 2.
  names(mat)

  ## Compare coding in the VCF file to the SnpMatrix.
  geno(vcf)$GT
  t(as(mat$genotype, "character"))

  ## The 'ignore' column in 'map' indicates which variants 
  ## were set to NA. Variant rs6040355 was ignored because 
  ## it has multiple alternate alleles, microsat1 is not a 
  ## snp, and chr20:1230237 has no alternate allele.
  mat$map

  ## ----------------------------------------------------------------
  ## Probability-based encoding using "GL" or "GP"
  ## ----------------------------------------------------------------
  ## Read a vcf file with a "GL" field.
  fl <- system.file("extdata", "gl_chr1.vcf", package="VariantAnnotation") 
  vcf <- readVcf(fl, "hg19")
  geno(vcf)

  ## Convert the "GL" FORMAT field to a SnpMatrix
  mat <- genotypeToSnpMatrix(vcf, uncertain=TRUE)

  ## Only 3 of the 9 variants passed the filters.
  mat$map

  ## Compare original vs SnpMatrix coding for a subset of
  ## samples in variant rs180734498. 
  geno(vcf)$GL["rs180734498", 14:16]
  t(as(mat$genotype, "character"))["rs180734498", 14:16]
}

\keyword{manip}
\keyword{methods}
