\name{VCF-class}
\docType{class}

% Class
\alias{class:VCF}
\alias{VCF-class}
\alias{CollapsedVCF}
\alias{class:CollapsedVCF}
\alias{CollapsedVCF-class}
\alias{ExpandedVCF}
\alias{class:ExpandedVCF}
\alias{ExpandedVCF-class}

% Constructors:
\alias{VCF}
\alias{updateObject,VCF-method}

% Coercion methods:

% Accessor methods:
\alias{ref}
\alias{ref,VCF-method}
\alias{ref<-}
\alias{ref<-,VCF,DNAStringSet-method}
\alias{alt}
\alias{alt,VCF-method}
\alias{alt<-}
\alias{alt<-,CollapsedVCF,CharacterList-method}
\alias{alt<-,ExpandedVCF,character-method}
\alias{alt<-,CollapsedVCF,DNAStringSetList-method}
\alias{alt<-,ExpandedVCF,DNAStringSet-method}
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
\alias{strand,VCF-method}
\alias{strand<-,VCF-method}
\alias{header,VCF-method}

% Subset methods
\alias{[,VCF-method}
\alias{[,VCF,ANY,ANY-method}
\alias{[<-,VCF,ANY,ANY,VCF-method}

% Other methods
\alias{genome,VCF-method}
\alias{seqlevels,VCF-method}
\alias{renameSeqlevels,VCF,character-method}
\alias{keepSeqlevels,VCF,character-method}
\alias{expand,CollapsedVCF-method}

% show method:
\alias{show,VCF-method}
\alias{show,CollapsedVCF-method}
\alias{show,ExpandedVCF-method}

\title{VCF class objects}

\description{
  The VCF class is a virtual class extended from \code{SummarizedExperiment}. 
  The subclasses, \code{CompressedVCF} and \code{ExtendedVCF}, are containers 
  for holding data from Variant Call Format files.
}

\section{Constructors}{
  \describe{
    \item{}{
      \code{readVcf(rowData = GRanges(), colData = DataFrame(), 
                exptData = SimpleList(), fixed = DataFrame(), 
                info = DataFrame(), geno = SimpleList(), 
                ..., collapsed=TRUE, verbose = FALSE)}
      Creates a CollapsedVCF class from data in a Variant Call
      Format file on disk. 
    }
    \item{}{
      \code{expand(rowData = GRanges(), colData = DataFrame(), 
                exptData = SimpleList(), fixed = DataFrame(), 
                info = DataFrame(), geno = SimpleList(), 
                ..., collapsed=TRUE, verbose = FALSE)}
      Creates an ExpandedVCF from a CollapsedVCF.
    }
    \item{}{
      \code{VCF(rowData = GRanges(), colData = DataFrame(), 
                exptData = SimpleList(), fixed = DataFrame(), 
                info = DataFrame(), geno = SimpleList(), 
                ..., collapsed=TRUE, verbose = FALSE)}
      Low-level constructor. Generally not used directly. The
      \code{collapsed} argument dictates the creation of a 
      CollapsedVCF or ExpandedVCF.
    }
  }
}

\section{Accessors}{
  In the following code snippets \code{x} is a CollapsedVCF or ExpandedVCF 
  object.
  \describe{
    \item{}{
      \code{rowData(x)}, \code{rowData(x) <- value}:
      Gets or sets a \code{GRanges} constructed from the CHROM, POS and ID 
      fields of the VCF file with REF, ALT, QUAL and FILT as metadata columns. 
      The ID's serve as the \code{rownames}; if they are NULL, \code{rownames} 
      are constructed from  CHROM:POS. \code{value} must be a \code{GRanges} 
      with metadata columns.

      The \code{paramRangeID} metadata column groups ranges as
      specified with \code{ScanVcfParam}. If \code{ScanVcfParam} was
      not used with \code{readVcf} these values will be \code{NA}. 

      The metadata columns can be accessed with \code{mcols(x)$<variable>} or 
      the following accessors:
      \itemize{
        \item{\code{ref(x)}, \code{ref(x) <- value}:
          Gets or sets the reference allele (REF). \code{value} must 
          be a \code{DNAStringSet}. 
        }
        \item{\code{alt(x)}, \code{alt(x) <- value}:
          Gets or sets the alternate allele data (ALT). When \code{x} is 
          a CollapsedVCF, \code{value} must be a \code{DNAStringSetList}
          or \code{CompressedCharacterList}. For ExtendedVCF, \code{value}
          must be a \code{DNAStringSet} or \code{character}.
        }
        \item{\code{qual(x)}, \code{qual(x) <- value}:
          Returns or sets the quality scores (QUAL). \code{value} must 
          be an \code{integer(1L)}.
        }
        \item{\code{filt(x)}, \code{filt(x) <- value}:
          Returns or sets the filter data (FILT). \code{value} must 
          be a \code{character(1L)}.
        }
      }
    }
    \item{}{
      \code{info(x)}, \code{info(x) <- value}:
       Gets or sets a DataFrame of INFO variables. 
    }
    \item{}{
      \code{geno(x)}, \code{geno(x) <- value}:
       Gets a SimpleList of genotype data.
       \code{value} is a SimpleList. To replace a single variable in
       the SimpleList use \code{geno(x)$variable <- value}; in this 
       case \code{value} must be a matrix or array.
    }
    \item{}{
      \code{exptData(x)}:
       Gets a SimpleList of header information from the VCF. For
       detailed descriptions of info, geno, etc.,
       \itemize{
         \item{\code{info(exptData(x)$header)}}
         \item{\code{geno(exptData(x)$header)}}
         \item{\code{meta(exptData(x)$header)}}
         \item{\code{samples(exptData(x)$header)}}
      }
    }
    \item{}{
      \code{colData(x)}:
      Gets a \code{DataFrame} of sample-specific information. Each row 
      represents a sample in the VCf file. \code{value} must be a 
      \code{DataFrame} with rownames representing the samples in the VCF 
      file.
    }
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
      \code{strand(x)}:
      Extract the \code{strand} from the \code{GRanges} in the 
      \code{rowData} slot.
    }
    \item{}{
      \code{header(x)}:
      Extract the VCF header information. The slot names of the resulting
      \code{VCFHeader} object are accessors that will further extract
      the data, (i.e., info(), geno(), fixed(), samples(), meta()).
    }
  }
}

\section{Subsetting and related operations}{
  In the following code snippets \code{x} is a VCF object.
  \describe{
    \item{}{
      \code{x[i, j]}, \code{x[i, j] <- value}:
      Gets or sets rows and columns. \code{i} and \code{j} can be integer
      or logical vectors. \code{value} is a replacement \code{VCF} object.
    }
  }
}

\section{Other methods}{
  In the following code snippets \code{x} is a VCF object.
  \describe{
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

\section{expand}{
  In the following code snippets \code{x} is a CollapsedVCF object.
  \describe{
    \item{}{
      \code{expand(x)}:
      Expand (unlist) the ALT column of a CollapsedVCF object to one row 
      per ALT value. In addition to ALT, all variables with Number='A' in the 
      header are also expanded. The 'A' indicates one value per alternate
      allele. Data without Number='A' in the header are replicated
      the appropriate number of times. 

      The output is an ExpandedVCF object with ALT as a \code{DNAStringSet} 
      or \code{character} (structural variants). 
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
    \item{collapsed}{A \code{logical(1)} indicating whether a 
      CollapsedVCF or ExpandedVCF should be created. The ALT in a
      CollapsedVCF is a \code{DNAStringSetList} while in a
      ExpandedVCF it is a \code{DNAStringSet}. 
    }
    \item{verbose}{A \code{logical(1)} indicating whether messages
      about data coercion during construction should be printed.
    }
  }
}

\details{
  The \code{VCF} class is a virtual class with two concrete subclasses,
  \code{CollapsedVCF} and \code{ExtendedVCF}. 

  Slots unique to \code{VCF} and subclasses,
  \itemize{
    \item \code{fixed}: A \link{DataFrame} containing information from the 
           REF, ALT, QUAL and FILTER fields from a VCF file.
    \item \code{info}: A \link{DataFrame} containing information from the INFO 
           fields from a VCF file.
  }

  Slots inherited from \code{SummarizedExperiment},
  \itemize{
    \item \code{exptData}: A \link{SimpleList}-class instance containing the 
          file header or other information about the overall experiment.
    \item \code{rowData}: A \link{GRanges}-class instance defining the
           variant ranges and associated metadata columns of REF, ALT, QUAL 
           and FILTER.
    \item \code{colData}: A \link{DataFrame}-class instance describing the 
          samples and associated metadata.
    \item \code{geno}: The \code{assays} slot from \code{SummarizedExperiment}
          has been renamed as \code{geno} for the VCF class. This slot contains
          the genotype information immediately following the FORMAT field in a 
          VCF file. Each element of the \code{list} or \code{SimpleList} is a 
          matrix or array. 
  }

  It is expected that users will not create instances of the VCF class
  but instead one of the concrete subclasses, CollapsedVCF or ExpandVCF.
  CollapsedVCF contains the ALT data as a \code{DNAStringSetList} allowing 
  for multiple alleles per variant. The ExpandedVCF ALT data is a 
  \code{DNAStringSet} where the ALT column has been expanded to create a 
  flat form of the data with one row per variant-allele combination. In 
  the case of strucutral variants, ALT will be a \code{CompressedCharacterList}
  or \code{character} in the collapsed or expanded forms.
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

  ## ----------------------------------------------------------------
  ## Accessors 
  ## ----------------------------------------------------------------
  ## Variant locations are stored in a GRanges in the rowData slot.
  rowData(vcf)

  ## Individual fields can be extracted with ref(), alt(), qual() 
  ## and filt().
  qual(vcf)
  ref(vcf)

  ## All 'info' fields can be extracted together along
  ## with the GRanges of locations.
  head(info(vcf))

  ## All genotype fields can be seen with geno(). Individual
  ## fields are accessed with '$' or '[['.
  geno(vcf)
  identical(geno(vcf)$DP, geno(vcf)[[3]])

  ## ----------------------------------------------------------------
  ## Renaming seqlevels and subsetting 
  ## ----------------------------------------------------------------
  ## It is often the case that chromosome names need to be
  ## modified. In this VCF, the chromosome names (seqlevels) 
  ## are numbers instead of preceded by "chr" or "ch".
  seqlevels(vcf)
 
  ## seqlevels can be renamed with renameSeqlevels().
  newnames <- paste("chr", seqlevels(vcf), sep="")
  names(newnames) <- seqlevels(vcf)
  vcf <- renameSeqlevels(vcf, newnames)
  seqlevels(vcf)
 
  ## Subsetting on seqlevels can be done with keepSeqlevels().
  vcf_subset <- keepSeqlevels(vcf, "chr2")
  seqlevels(vcf_subset)

  ## ----------------------------------------------------------------
  ## CollapsedVCF 
  ## ----------------------------------------------------------------
  ## readVCF() produces a CollapsedVCF object.
  fl <- system.file("extdata", "ex2.vcf", 
                    package="VariantAnnotation")
  vcf <- readVcf(fl, genome="hg19")
  vcf

  ## The ALT column is a DNAStringSetList to allow for more
  ## than one alternate allele per variant.
  alt(vcf)

  ## When structural variants are present, ALT is 
  ## a CompressedCharacterList.
  fl <- system.file("extdata", "structural.vcf", 
                    package="VariantAnnotation")
  vcf <- readVcf(fl, genome="hg19")
  alt(vcf)

  ## ----------------------------------------------------------------
  ## ExpandedVCF
  ## ----------------------------------------------------------------
  ## ExpandedVCF is the 'flattened' counterpart of CollapsedVCF.
  ## The ALT and all variables with Number='A' in the header are
  ## expanded to one row per alternate allele.
  fl <- system.file("extdata", "ex2.vcf", package="VariantAnnotation")
  vcf <- readVcf(fl, genome="hg19")
  vcfLong <- expand(vcf)
  alt(vcfLong)
}
