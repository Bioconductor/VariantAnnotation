\name{VCF-class}
\docType{class}

% Class
\alias{VCF-class}

% Constructor-like function:
\alias{VCF}

% Coercion methods:

% Accessor methods:
\alias{ref}
\alias{ref,VCF-method}
\alias{ref<-}
\alias{ref<-,VCF,DNAStringSet-method}
\alias{alt}
\alias{alt,VCF-method}
\alias{alt<-}
\alias{alt<-,VCF,CharacterList-method}
\alias{alt<-,VCF,DNAStringSetList-method}
\alias{qual}
\alias{qual,VCF-method}
\alias{qual<-}
\alias{qual<-,VCF,integer-method}
\alias{filt}
\alias{filt,VCF-method}
\alias{filt<-}
\alias{filt<-,VCF,character-method}
\alias{fixed}
\alias{fixed,VCF-method}
\alias{fixed<-}
\alias{fixed<-,VCF,DataFrame-method}
\alias{info}
\alias{info,VCF-method}
\alias{info<-}
\alias{info<-,VCF,DataFrame-method}
\alias{geno}
\alias{geno,VCF-method}
\alias{geno<-}
\alias{geno<-,VCF,character,matrix-method}
\alias{geno<-,VCF,numeric,matrix-method}
\alias{geno<-,VCF,missing,SimpleList-method}

% Subset methods
\alias{[,VCF-method}
\alias{[,VCF,ANY,ANY,ANY-method}
\alias{[<-,VCF,ANY,ANY,VCF-method}

% Other methods
\alias{genome,VCF-method}
\alias{seqlevels,VCF-method}
\alias{renameSeqlevels,VCF,character-method}
\alias{keepSeqlevels,VCF,character-method}

% show method:
\alias{show,VCF-method}

\title{VCF instances}

\description{
  The VCF class is an extension of the
  \code{\linkS4class{SummarizedExperiment}-class} with two additional slots,
  \code{fixed} and \code{info}.
}

\section{Extends}{
   Directly extends class \code{\link{SummarizedExperiment}}.
}

\section{Constructor}{
  \describe{
    \item{}{
      \code{VCF(rowData = GRanges(), colData = DataFrame(), 
                exptData = SimpleList(), fixed = DataFrame(), 
                info = DataFrame(), geno = SimpleList(), 
                ..., verbose = FALSE)
      }
    }
  }
}

\section{Accessors}{
  In the following code snippets \code{x} is a VCF object. All accessors except
  \code{geno} return the data as elementMetadata column(s) of the \code{rowData} 
  \code{GRanges} object. The \code{geno} accessor returns a \code{SimpleList}.

  \describe{
    \item{}{
      \code{ref(x)}, \code{ref(x) <- value}:
      Returns or sets the reference allele data from the REF column of the VCF
      file. \code{value} must be a \code{DNAStringSet}. 
    }
    \item{}{
      \code{alt(x)}, \code{alt(x) <- value}:
      Returns or sets the alternate allele data from the ALT column of the VCF
      file. \code{value} can be a \code{DNAStringSet} or a \code{CharacterList}
      (for a structural VCF file). 
    }
    \item{}{
      \code{qual(x)}, \code{qual(x) <- value}:
      Returns or sets the quality scores from the QUAL column of the VCF
      file. \code{value} must be an \code{integer(1L)}.
    }
    \item{}{
      \code{filt(x)}, \code{filt(x) <- value}:
      Returns or sets the filter data from the FILTER column of the VCF
      file. \code{value} must be a \code{character(1L)}.
    }
    \item{}{
      \code{fixed(x)}, \code{fixed(x) <- value}:
      Returns or sets a \code{DataFrame} of the REF, ALT, QUAL, and FILTER fields 
      from the VCF file. \code{value} must be a \code{DataFrame}.
    }
    \item{}{
      \code{info(x)}, \code{info(x) <- value}:
      Returns or sets \code{info} data. Contains the information
      stored in the INFO field of a VCF file. \code{value} must be a 
      \code{DataFrame}. If no \code{info} fields are present the
      \code{rowData} \code{GRanges} will be returned with no elementMetadata
      columns. 
    }
    \item{}{
      \code{geno(x)}, \code{geno(x) <- value}:
      Returns or sets \code{geno} data. Contains the genotype information 
      from the samples in a VCF file. \code{value} must be a \code{SimpleList} 
      of \code{matrices} or \code{arrays}. An optional \code{withDimnames}
      argument controls the return of dimnames (default = TRUE).

      \code{geno(x)[[i]]}, \code{geno(x)[[i]] <- value}:
      Returns or sets elements of geno. \code{value} can be a \code{matrix},
      or \code{array}. 
    }
    \item{}{
      \code{exptData(x)}, \code{exptData(x) <- value}:
      Returns or sets \code{exptData} data. Contains the header information 
      from a VCF file as well as any other experiment-specific information.
      \code{value} must be a \code{SimpleList}.
    }
    \item{}{
      \code{rowData(x)}, \code{rowData(x) <- value}:
      Returns or sets \code{rowData} data. Contains a \code{GRanges} constructed
      from the CHROM, POS and ID fields of the VCF file. The ID's serve as the
      \code{rownames}; if they are NULL, \code{rownames} are constructed from 
      CHROM:POS. \code{value} must be a \code{GRanges} with names representing 
      the ID's in the VCF file.
    }
    \item{}{
      \code{colData(x)}, \code{colData(x) <- value}:
      Returns or sets \code{colData} data. Contains a \code{DataFrame} of 
      sample-specific information. Each row represents a sample in the VCF
      file. \code{value} must be a \code{DataFrame} with rownames representing
      the samples in the VCF file.
    }
  }
}

\section{Subsetting}{
  In the following code snippets \code{x} is a VCF object.

  \describe{
    \item{}{
      \code{x[i, j]}, \code{x[i, j] <- value}:
      Gets or sets rows \code{i} and columns \code{j}. \code{i} and \code{j}
      can be integer or logical vectors. \code{value} is a replacement
      \code{VCF} object.
    }
  }
}

\section{Other methods}{
  \describe{
    \item{}{
      \code{genome(x)}:
      Extract the \code{genome} information from the \code{GRanges} in 
      the \code{rowData} slot.
    }
    \item{}{
      \code{seqlevels(x)}:
      Extract the \code{seqlevels} from the \code{GRanges} in the 
      \code{rowData} slot.
    }
    \item{}{
      \code{renameSeqlevels(x, value)}:
      Rename the seqlevels in the \code{GRanges} in the \code{rowData} slot
      of the VCF object. \code{value} is a named character vector where
      the names are the old seqlevels and the values are the new.
    }
    \item{}{
      \code{keepSeqlevels(x, value)}:
      Subset the \code{GRanges} in the \code{rowData} slot of the VCF object. 
      \code{value} is a character vector of seqlevels to keep.
    }
  }
}

\section{Arguments}{
  \describe{
    \item{geno}{A \code{list} or \code{SimpleList} of matrix elements,
      or a \code{matrix} containing the genotype information from a
      VCF file. If present, these data immediately follow the FORMAT
      field in the VCF. 

      Each element of the list must have the same dimensions, and dimension 
      names (if present) must be consistent across elements and with the row 
      names of \code{rowData}, \code{colData}. 
    }
    \item{info}{A \code{DataFrame} of data from the INFO field of a VCF file. 
      The number of rows must match that in the \code{rowData} object.
    }
    \item{fixed}{A \code{DataFrame} of REF, ALT, QUAL and FILTER 
      fields from a VCF file. The number of rows must match that
      of the \code{rowData} object.
    }
    \item{rowData}{A \code{GRanges} instance describing the ranges of interest. 
      Row names, if present, become the row names of the \code{VCF}. The length 
      of the \code{GRanges} must equal the number of rows of the matrices in 
      \code{geno}.
    }
    \item{colData}{A \code{DataFrame} describing the samples. Row names, if 
      present, become the column names of the \code{VCF}.
    }
    \item{exptData}{A \code{SimpleList} describing the header of the VCF file or 
      additional information for the overall experiment. 
    }
    \item{...}{Additional arguments passed to methods. 
    }
    \item{withDimnames}{A \code{logical(1)}, indicating whether dimnames
      should be applied to extracted assay elements. Applicable to the 
      \code{geno} accessor only.
    }
    \item{verbose}{A \code{logical(1)} indicating whether messages
      about data coercion during construction should be printed.
    }
  }
}


\details{
  The \code{VCF} class is designed to hold data from a Variant Call Format (VCF) 
  file. The class extends \code{\linkS4class{SummarizedExperiment}} with the 
  addition of two slots, \code{fixed} and \code{info}. See 
  ?\code{SummarizedExperiment} for a detailed description of the inherited
  slots. 

  Slots unique to the \code{VCF} class,
  \describe{
    \item{\code{fixed}}{A \link{DataFrame} containing information from the 
      REF, ALT, QUAL and FILTER fields from a VCF file.
    }
    \item{\code{info}}{A \link{DataFrame} containing information from the INFO 
      fields from a VCF file.
    }
  }

  Slots inherited from the \code{SummarizedExperiment} class,
  \describe{
    \item{\code{exptData}}{A \link{SimpleList}-class instance containing the 
      file header or other information about the overall experiment.
    }
    \item{\code{rowData}}{A \link{GRanges}-class instance defining the
      variant ranges and associated metadata columns of REF, ALT, QUAL and
      FILTER.
    }
    \item{\code{colData}}{A \link{DataFrame}-class instance describing the 
      samples and associated metadata.
    }
    \item{\code{geno}}{The \code{assays} slot from \code{SummarizedExperiment}
     has been renamed as \code{geno} for the VCF class. This slot contains
     the genotype information immediately following the FORMAT field in a VCF
     file. Each element of the \code{list} or \code{SimpleList} is a matrix or 
     array. 
    }
  }

}

\author{Valerie Obenchain <vobencha@fhcrc.org>}

\seealso{
  \link{GRanges},
  \link[IRanges]{DataFrame},
  \link[IRanges]{SimpleList},
  \link[GenomicRanges]{SummarizedExperiment},
  \code{\link{readVcf}},
  \code{\link{writeVcf}}
}

\examples{
  fl <- system.file("extdata", "structural.vcf", package="VariantAnnotation")
  vcf <- readVcf(fl, genome="hg19")

  ## Access all 'fixed' fields with fixed and individual fixed fields with
  ## ref, alt, qual and filt accessors
  fixed(vcf)
  ref(vcf)

  ##  Extract 'fixed', 'info' and 'geno' fields 
  head(fixed(vcf))
  head(info(vcf))
  geno(vcf)
  identical(geno(vcf)$DP, geno(vcf)[[3]])

  ## Rename and subset the seqlevels in the rowData GRanges object :
  ## The seqlevels (chromosome names) of the vcf are numbers instead of 
  ## being preceded by "chr*" or "ch*"
  seqlevels(vcf)
 
  ## Add "chr" in front of the seqlevels using renameSeqlevels()
  newnames <- paste("chr", seqlevels(vcf), sep="")
  names(newnames) <- seqlevels(vcf)
  vcf <- renameSeqlevels(vcf, newnames)
  seqlevels(vcf)
 
  ## To subset on seqlevels use the keepSeqlevels() helper function
  vcf_subset <- keepSeqlevels(vcf, "chr2")
  seqlevels(vcf_subset)
}
