\name{isSNV}

\alias{isSNV}
\alias{isInsertion}
\alias{isDeletion}
\alias{isIndel}
\alias{isDelins}
\alias{isTransition}
\alias{isSubstitution}

% VCF
\alias{isSNV,ExpandedVCF-method}
\alias{isSNV,CollapsedVCF-method}
\alias{isInsertion,ExpandedVCF-method}
\alias{isInsertion,CollapsedVCF-method}
\alias{isDeletion,ExpandedVCF-method}
\alias{isDeletion,CollapsedVCF-method}
\alias{isIndel,ExpandedVCF-method}
\alias{isIndel,CollapsedVCF-method}
\alias{isDelins,ExpandedVCF-method}
\alias{isDelins,CollapsedVCF-method}
\alias{isSubstitution,ExpandedVCF-method}
\alias{isSubstitution,CollapsedVCF-method}
\alias{isTransition,ExpandedVCF-method}
\alias{isTransition,CollapsedVCF-method}

% VRanges 
\alias{isSNV,VRanges-method}
\alias{isInsertion,VRanges-method}
\alias{isDeletion,VRanges-method}
\alias{isIndel,VRanges-method}
\alias{isDelins,VRanges-method}
\alias{isSubstitution,VRanges-method}
\alias{isTransition,VRanges-method}

\title{
  Identification of genomic variant types.
}

\description{
  Functions for identifying variant types such as SNVs, insertions, 
  deletions, transitions, and structural rearrangements. 
}

\usage{
\S4method{isSNV}{VRanges}(x, ...)
\S4method{isSNV}{ExpandedVCF}(x, ...)
\S4method{isSNV}{CollapsedVCF}(x, ..., singleAltOnly = TRUE)

\S4method{isInsertion}{VRanges}(x, ...)
\S4method{isInsertion}{ExpandedVCF}(x, ...)
\S4method{isInsertion}{CollapsedVCF}(x, ..., singleAltOnly = TRUE)

\S4method{isDeletion}{VRanges}(x, ...)
\S4method{isDeletion}{ExpandedVCF}(x, ...)
\S4method{isDeletion}{CollapsedVCF}(x, ..., singleAltOnly = TRUE)

\S4method{isIndel}{VRanges}(x, ...)
\S4method{isIndel}{ExpandedVCF}(x, ...)
\S4method{isIndel}{CollapsedVCF}(x, ..., singleAltOnly = TRUE)

\S4method{isDelins}{VRanges}(x, ...)
\S4method{isDelins}{ExpandedVCF}(x, ...)
\S4method{isDelins}{CollapsedVCF}(x, ..., singleAltOnly = TRUE)

\S4method{isTransition}{VRanges}(x, ...)
\S4method{isTransition}{ExpandedVCF}(x, ...)
\S4method{isTransition}{CollapsedVCF}(x, ..., singleAltOnly = TRUE)

\S4method{isSubstitution}{VRanges}(x, ...)
\S4method{isSubstitution}{ExpandedVCF}(x, ...)
\S4method{isSubstitution}{CollapsedVCF}(x, ..., singleAltOnly = TRUE)
}

\arguments{
  \item{x}{A \linkS4class{VCF} or \linkS4class{VRanges} object.
  }
  \item{singleAltOnly}{A \code{logical} only applicable when \code{x}
    is a \link{CollapsedVCF} class.

    When \code{TRUE} (default) only variants with a single alternate 
    allele are evaluated. When \code{FALSE} all ref / alt pairs for 
    each variant are evaluated. If any ref / alt pairs meet the 
    test criteria a value of TRUE is returned for the variant; this
    may result in a value of TRUE for a variant with a mixture of
    alternate alleles, some that pass the criteria and some that do not.
    To retain single ref / alt pairs that pass the critera use
    \code{expand} on the \code{CollapsedVCF} and then apply the test.
  }
  \item{\dots}{Arguments passed to other methods.
  }
}

\details{
  All functions return a logical vector the length of \code{x}.
  \itemize{
    \item{isSNV: }{
      Reference and alternate alleles are both a single nucleotide long.
    }
    \item{isInsertion: }{
      Reference allele is a single nucleotide and the alternate allele 
      is greater (longer) than a single nucleotide and the first
      nucleotide of the alternate allele matches the reference. 
    }
    \item{isDeletion: }{
      Alternate allele is a single nucleotide and the reference allele 
      is greater (longer) than a single nucleotide and the first
      nucleotide of the reference allele matches the alternate. 
    }
    \item{isIndel: }{
      The variant is either a deletion or insertion as determined
      by \code{isDeletion} and \code{isInsertion}.
    }
    \item{isDelins: }{
      The variant is a deletion followed by an insertion, either of them
      involving two or more nucleotides.
    }
    \item{isSubstition: }{
      Reference and alternate alleles are the same length (1 or
      more nucleotides long).
    }
    \item{isTransition: }{
      Reference and alternate alleles are both a single nucleotide long.
      The reference-alternate pair interchange is of either two-ring 
      purines (A <-> G) or one-ring pyrimidines (C <-> T).
    }
  }
}

\value{
  A \code{logical} vector the same length as \code{x}.
}

\author{Michael Lawrence, Valerie Obenchain and Robert Castelo}

\examples{
  fl <- system.file("extdata", "ex2.vcf", package="VariantAnnotation")
  ## ---------------------------------------------------------------------
  ## VCF objects 
  ## ---------------------------------------------------------------------
  vcf <- readVcf(fl, "hg19")
  DataFrame(ref(vcf), alt(vcf))

  ## This vcf has transitions in row 2 and 3. When 'singleAltOnly=TRUE' 
  ## only the row 2 variant is identified:
  isTransition(vcf)

  ## Both row 2 and 3 are identified when 'singleAltOnly=FALSE':
  isTransition(vcf, singleAltOnly=FALSE)

  ## Expand the CollapsedVCF to ExpandedVCF
  evcf <- expand(vcf)
 
  ## All ref / alt pairs are now expanded and there is no need to 
  ## use 'singleAltOnly'. The return length is now 7 instead of 5:
  transition <- isTransition(evcf)
  transition
  DataFrame(ref(evcf)[transition], alt(evcf)[transition])
 
  ## ---------------------------------------------------------------------
  ## VRanges objects 
  ## ---------------------------------------------------------------------
  ## A VRanges object holds data from a VCF class in a completely
  ## 'flat' fashion. INFO and FORMAT variables for all subjects are
  ## 'repped' out such that each row is a unique combination of data.
  vr <- as(vcf, "VRanges")
  isSNV(vr, singleAltOnly=FALSE)
}

\keyword{methods}
