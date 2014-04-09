\name{VRanges-class}
\docType{class}

% Class:
\alias{class:VRanges}
\alias{VRanges-class}
\alias{VRanges}

% Constructors:
\alias{VRanges}

% Coercion:
\alias{asVCF}
\alias{asVCF,VRanges-method}
\alias{coerce,VRanges,VCF-method}
\alias{coerce,VCF,VRanges-method}

% Accessors:
\alias{alt,VRanges-method}
\alias{alt<-,VRanges,ANY-method}
\alias{ref,VRanges-method}
\alias{ref<-,VRanges,ANY-method}
\alias{altDepth}
\alias{altDepth,VRanges-method}
\alias{altDepth<-}
\alias{altDepth<-,VRanges-method}
\alias{refDepth}
\alias{refDepth,VRanges-method}
\alias{refDepth<-}
\alias{refDepth<-,VRanges-method}
\alias{totalDepth}
\alias{totalDepth,VRanges-method}
\alias{totalDepth<-}
\alias{totalDepth<-,VRanges-method}
\alias{altFraction}
\alias{altFraction,VRanges-method}
\alias{called}
\alias{called,VRanges-method}
\alias{hardFilters<-}
\alias{hardFilters<-,VRanges-method}
\alias{hardFilters}
\alias{hardFilters,VRanges-method}
\alias{sampleNames,VRanges-method}
\alias{sampleNames<-,VRanges,ANY-method}
\alias{softFilterMatrix}
\alias{softFilterMatrix,VRanges-method}
\alias{softFilterMatrix<-}
\alias{softFilterMatrix<-,VRanges-method}
\alias{resetFilter}

% Aggregation:
\alias{tabulate}
\alias{tabulate,VRanges-method}

% VCF reading/writing:
\alias{writeVcf,VRanges,ANY-method}
\alias{VRangesScanVcfParam}
\alias{readVcfAsVRanges}

% Utilities:
\alias{match,VRanges,VRanges-method}
\alias{softFilter}

% Typed Rle classes (at least for now)
\alias{characterRle-class}
\alias{characterOrRle-class}
\alias{complexRle-class}
\alias{factorRle-class}
\alias{factorOrRle-class}
\alias{integerRle-class}
\alias{integerOrRle-class}
\alias{logicalRle-class}
\alias{numericRle-class}
\alias{rawRle-class}

\title{VRanges objects}

\description{
  The VRanges class is a container for variant calls, including SNVs and
  indels. It extends \code{\link[GenomicRanges]{GRanges}} to provide
  special semantics on top of a simple vector of genomic locations. While
  it is not as expressive as the \code{\linkS4class{VCF}} object, it is
  a simpler alternative that may be convenient for variant
  calling/filtering and similar exercises.
}

\details{
  VRanges extends GRanges to store the following components. Except
  where noted, the components are considered columns in the dataset,
  i.e., their lengths match the number of variants. Many columns can be
  stored as either an atomic vector or an Rle. 
  \describe{
    \item{\code{ref}}{(\code{character}), the reference
      allele. The range (start/end/width) should always correspond to
      this sequence.}
    \item{\code{alt}}{(\code{character/Rle}),
      the alternative allele (NA allowed). By convention there is only
      a single alt allele per element (row) of the VRanges. Many methods,
      like \code{match}, make this assumption. 
    }
    \item{\code{refCount}}{(\code{integer/Rle}), read count for the
      reference allele (NA allowed)}
    \item{\code{altCount}}{(\code{integer/Rle}), read count for the
      alternative allele (NA allowed)}
    \item{\code{totalCount}}{(\code{integer/Rle}), total read count at the
      position, must be at least \code{refCount+altCount} (NA allowed)}
    \item{\code{sampleNames}}{(\code{factor/Rle}), name of the sample -
      results from multiple samplse can be combined into the same object
      (NA allowed)}
    \item{\code{softFilterMatrix}}{(\code{matrix/FilterMatrix}),
      variant by filter matrix, \code{TRUE} where variant passed the
      filter; use a \code{\link[IRanges]{FilterMatrix}} to store the
      actual \code{FilterRules} object that was applied}
    \item{\code{hardFilters}}{(\code{FilterRules}) record of hard
      filters applied, i.e., only the variants that passed the filters
      are present in this object; this is the only component that is not
      a column, i.e., its length does not match the number of variants}
  }
  Except in the special circumstances described here, a \code{VRanges}
  may be treated like a \code{GRanges}. The range should span the
  sequence in \code{ref}. Indels are typically represented by the VCF
  convention, i.e., the start position is one upstream of the event. The
  strand is always constrained to be positive (+).

  Indels, by convention, should be encoded VCF-style, with the upstream
  reference base prepended to the indel sequence. The ref/alt for a
  deletion of GCGT before A might be AGCGT/A and for an insertion might
  be A/AGCGT. Since the range always matches the \code{ref} sequence,
  this means a deletion will be the width of the deletion + 1, and an
  insertion is always of width 1.
}

\section{Constructor}{
  \describe{
    \item{}{
      \code{VRanges(seqnames = Rle(), ranges = IRanges(), ref = character(), 
        alt = NA_character_, totalDepth = NA_integer_, refDepth = NA_integer_, 
        altDepth = NA_integer_, ..., sampleNames = NA_character_, 
        softFilterMatrix = FilterMatrix(matrix(nrow = length(gr), 
        ncol = 0L), FilterRules()), hardFilters = FilterRules())}:
      Creates a VRanges object.
      \describe{
        \item{\code{seqnames}}{Rle object, character vector, or factor
          containing the sequence names.}
        \item{\code{ranges}}{IRanges object containing the ranges.}
        \item{\code{ref}}{character vector, containing the reference allele.}
        \item{\code{alt}}{character vector/Rle,
          containing the alternative allele (NA allowed).}
        \item{\code{totalDepth}}{integer vector/Rle, containing the
          total read depth (NA allowed).}
        \item{\code{refDepth}}{integer vector/Rle, containing the
          reference read depth (NA allowed).}
        \item{\code{altDepth}}{integer vector/Rle, containing the
          reference read depth (NA allowed).}
        \item{\code{\ldots}}{Arguments passed to the \code{GRanges}
          constructor.}
        \item{\code{sampleNames}}{character/factor vector/Rle, containing the
          sample names (NA allowed).}
        \item{\code{softFilterMatrix}}{a matrix (typically
          a \code{\link[IRanges]{FilterMatrix}}) of dimension variant by
          filter, with logical values indicating whether a variant
          passed the filter.}
        \item{\code{hardFilters}}{a \code{\link[IRanges]{FilterRules}},
          containing the filters that have already been applied to
          subset the object to its current state.}
      }
    }
  }
}

\section{Coercion}{
  These functions/methods coerce objects to and from \code{VRanges}:
  
  \describe{
    \item{}{
      \code{asVCF(x, info = character(), filter = character(), meta =
        character())}: Creates a VCF object from a VRanges object. The
      following gives the mapping from VRanges components to VCF:
      \describe{
        \item{seqnames(x)}{CHROM column}
        \item{start(x)}{POS column}
        \item{names(x)}{ID column}
        \item{ref(x)}{REF column}
        \item{alt(x)}{ALT column}
        \item{totalDepth(x)}{DP in FORMAT column}
        \item{altDepth(x), refDepth(x)}{AD in FORMAT column}
        \item{sampleNames(x)}{Names the sample columns}
        \item{softFilterMatrix(x)}{FT in FORMAT column, except filters
          named in \code{filter} argument, which are considered
          per-position and placed in the FILTER column}
        \item{hardFilters(x)}{Not yet exported}
        \item{mcols(x)}{Become fields in the FORMAT column; unless they
          are named in the \code{info} argument, in which case they
          are considered per-position and placed in the INFO column}
        \item{metadata(x)}{If named in the \code{meta} argument, output
          in the VCF header; a component is required to be coercible to
          a character vector of length one.}
      }

      Note that \code{identical(x, as(as(x, "VCF"), "VRanges"))}
      generally return \code{FALSE}.  During coercion to VCF, the "geno"
      components are reshaped into matrix form, with NAs filling the
      empty cells. The reverse coercion will not drop the NA values, so
      rows are added to the new VRanges. All logical values will become
      integers in VCF, and there is no automatic way of regenerating the
      logical column with the reverse coercion. There are many other
      cases of irreversibility.
    }
    \item{}{
      \code{as(from, "VCF")}: Like calling \code{asVCF(from)}.
    }
    \item{}{
      \code{as(from, "VRanges")}:
      Currently supported when \code{from} is a \code{VCF}. Essentially
      the inverse of \code{asVCF}. Information missing in the VCF
      is imputed as NA.
    }
  }
}

\section{Accessors}{
  In addition to all of the \code{GRanges} accessors, \code{VRanges}
  provides the following, where \code{x} is a VRanges object.

  \describe{
    \item{}{
      \code{alt(x), alt(x) <- value}: Get or set the alt allele (character).
    }
    \item{}{
      \code{ref(x), ref(x) <- value}: Get or set the ref allele (character).
    }
    \item{}{
      \code{altDepth(x), altDepth(x) <- value}: Get or set the alt allele
      read depth (integer).
    }
    \item{}{
      \code{refDepth(x), refDepth(x) <- value}: Get or set the ref
      allele read depth (integer).
    }
    \item{}{
      \code{totalDepth(x), totalDepth(x) <- value}: Get or set the total
      read depth (integer).
    }
    \item{}{
      \code{altFraction(x)}: Returns \code{altDepth(x)/totalDepth(x)} (numeric).
    }
    \item{}{
      \code{sampleNames(x), sampleNames(x) <- value}: Get or set the
      sample names (character/factor).
    }
    \item{}{
      \code{softFilterMatrix(x), softFilterMatrix(x) <- value}: Gets or
      sets the soft filter matrix (any matrix, but ideally a
      \code{FilterMatrix}).
    }
    \item{}{
      \code{resetFilter(x)}: Removes all columns from \code{softFilterMatrix}.
    }
    \item{}{
      \code{called(x)}: Returns whether all filter results in
      \code{softFilterMatrix(x)} are \code{TRUE} for each variant.
    }
    \item{}{
      \code{hardFilters(x), hardFilters(x) <- value}: Gets or
      sets the hard filters (those applied to yield the current subset).
    }
  }
}

\section{Utilities and Conveniences}{
  \describe{
    \item{}{
      \code{match(x)}: Like GRanges \code{match}, except matches on the
      combination of chromosome, start, width, and \strong{alt}.
    }
    \item{}{
      \code{tabulate(bin)}: Finds \code{unique(bin)} and counts how many
      times each unique element occurs in \code{bin}. The result is
      stored in \code{mcols(bin)$sample.count}.
    }
    \item{}{
      \code{softFilter(x, filters, ...)}: applies the \code{FilterRules}
      in \code{filters} to \code{x}, storing the results in
      \code{softFilterMatrix}.
    }
  }  
}

\section{Input/Output to/from VCF}{
  \describe{
    \item{}{
      \code{writeVcf(obj, filename, ...)}: coerces to a VCF object and
      writes it to a file; see \code{\link{writeVcf}}.
    }
    \item{}{
      \code{readVcfAsVRanges(x, genome, param = VRangesScanVcfParam(), ...)}:
      Reads a VCF \code{x} directly into a \code{VRanges};
      see \code{\link{readVcf}} for details on the arguments.
    }
    \item{}{
      \code{VRangesScanVcfParam(fixed = "ALT", info = NA, geno = "AD", ...)}:
      Convenience constructor for a \code{\linkS4class{ScanVcfParam}}
      object that is well suited for import to a \code{VRanges}. The
      \code{fixed}, \code{info} and \code{geno} parameters, and anything
      in \code{\dots} are passed to \code{\link{ScanVcfParam}}.
    }
  }
}

\section{Variant Type}{
  Functions to identify variant type include \link{isSNV}, 
  \link{isInsertion}, \link{isDeletion}, \link{isIndel}, 
  \link{isSubstitution} and \link{isTransition}. See the ?\code{isSNV} 
  man page for details.
}

\author{Michael Lawrence}

\seealso{
  \link{VRangesList}, a list of \code{VRanges}; \code{bam_tally} in the
  gmapR package, which generates a \code{VRanges}.
}

\examples{
## construction
vr <- VRanges(seqnames = c("chr1", "chr2"),
              ranges = IRanges(c(1, 10), c(5, 20)),
              ref = c("T", "A"), alt = c("C", "T"),
              refDepth = c(5, 10), altDepth = c(7, 6),
              totalDepth = c(12, 17), sampleNames = letters[1:2],
              hardFilters =
                FilterRules(list(coverage = function(x) totalDepth > 10)),
              softFilterMatrix =
                FilterMatrix(matrix = cbind(depth = c(TRUE, FALSE)),
                             FilterRules(depth = function(x) altDepth(x) > 6)),
              tumorSpecific = c(FALSE, TRUE))

## simple accessors
ref(vr)
alt(vr)
altDepth(vr)
vr$tumorSpecific
called(vr)

## coerce to VCF and write
vcf <- as(vr, "VCF")
## writeVcf(vcf, "example.vcf")
## or just
## writeVcf(vr, "example.vcf")

## other utilities
match(vr, vr[2:1])
}
