\name{snpSummary}

\alias{snpSummary}
\alias{snpSummary,CollapsedVCF-method}


\title{Counts and distribution statistics for SNPs in a VCF object}

\description{
  Counts and distribution statistics for SNPs in a VCF object
}

\usage{
  \S4method{snpSummary}{CollapsedVCF}(x, ...)
}

\arguments{
  \item{x}{
    A \link{CollapsedVCF} object.
  }
  \item{\dots}{
    Additional arguments to methods.
  }
}

\details{
  Genotype counts, allele counts and Hardy Weinberg equilibrium
  (HWE) statistics are calculated for single nucleotide variants
  in a \link{CollapsedVCF} object. HWE has been established as a 
  useful quality filter on genotype data. This equilibrium should 
  be attained in a single generation of random mating. Departures
  from HWE are indicated by small p values and are almost invariably 
  indicative of a problem with genotype calls.

  The following caveats apply:
  \itemize{
    \item No distinction is made between phased and unphased genotypes. 
    \item Only diploid calls are included.
    \item Only `valid' SNPs are included. A `valid' SNP is defined
          as having a reference allele of length 1 and a single 
          alternate allele of length 1.
  }
  Variants that do not meet these criteria are set to NA. 
}

\value{
  The object returned is a \code{data.frame} with seven columns.
  \describe{
    \item{g00}{
      Counts for genotype 00 (homozygous reference).
    }
    \item{g01}{
      Counts for genotype 01 or 10 (heterozygous).
    }
    \item{g11}{
      Counts for genotype 11 (homozygous alternate).
    }
    \item{a0Freq}{
      Frequency of the reference allele.
    }
    \item{a1Freq}{
      Frequency of the alternate allele.
    }
    \item{HWEzscore}{
      Z-score for departure from a null hypothesis of Hardy Weinberg equilibrium.
    }
    \item{HWEpvalue}{
      p-value for departure from a null hypothesis of Hardy Weinberg equilibrium.
    }
  }
}

\author{
  Chris Wallace <cew54@cam.ac.uk>
}

\seealso{
  \link{genotypeToSnpMatrix},
  \link{probabilityToSnpMatrix}
}

\examples{
  fl <- system.file("extdata", "ex2.vcf", package="VariantAnnotation")
  vcf <- readVcf(fl, "hg19")

  ## The return value is a data.frame with genotype counts
  ## and allele frequencies.
  df <- snpSummary(vcf)
  df

  ## Compare to ranges in the VCF object:
  rowRanges(vcf)

  ## No statistics were computed for the variants in rows 3, 4 
  ## and 5. They were omitted because row 3 has two alternate 
  ## alleles, row 4 has none and row 5 is not a SNP.
}

\keyword{manip}
