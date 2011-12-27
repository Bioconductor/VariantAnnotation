\name{VCF-class}
\docType{class}

% Class
\alias{VCF-class}
\alias{VCF}

% Constructor-like function:
\alias{VCF}

% Coercion methods:

% Accessor methods:
\alias{info}
\alias{info,VCF,missing-method}
\alias{info,VCF,numeric-method}
\alias{info,VCF,character-method}
\alias{info<-}
\alias{info<-,VCF,missing,list-method}
\alias{info<-,VCF,missing,SimpleList-method}
\alias{info<-,VCF,numeric,Vector-method}
\alias{info<-,VCF,character,Vector-method}
\alias{geno}
\alias{geno,VCF,ANY-method}
\alias{geno,VCF,missing-method}
\alias{geno<-}
\alias{geno<-,VCF,ANY-method}
\alias{geno<-,VCF,missing-method}

% Subset methods
\alias{[,VCF-method}
\alias{[,VCF,ANY,ANY,ANY-method}
\alias{[<-,VCF,ANY,ANY,VCF-method}

% show method:
\alias{show,VCF-method}

\title{VCF instances - under construction -}

\description{
  The VCF class is an extension of the
  \code{linkS4class{SummarizedExperiment}-class} with an additional slot,
  \code{info}, which contains a \code{SimpleList}.
}

\section{Extends}{
   Directly extends class \code{\link{SummarizedExperiment}}.
}

\section{Constructor}{
  \describe{
    \item{}{
      \code{VCF(geno = SimpleList(), rowData = GRanges(), 
                colData = DataFrame(), exptData = SimpleList(), 
                info = SimpleList(), ..., verbose = FALSE)
      }
    }
  }
}

\section{Accessors}{
  In the following code snippets \code{x} is a VCF object.

  \describe{
    \item{}{
      \code{info(x, ..., withDimnames=TRUE)}, \code{info(x) <- value}:
      Returns or sets \code{info} data. Contains the information
      stored in the INFO field of a VCF file. \code{value} must be a 
      \code{SimpleList}.

      \code{info(x)[[i]]}, \code{info(x)[[i]] <- value}:
      Returns or sets elements of info. \code{value} can be a \code{Vector},
      \code{CompressedList} or \code{array}. 
    }
    \item{}{
      \code{geno(x, ..., withDimnames=TRUE)}, \code{geno(x) <- value}:
      Returns or sets \code{geno} data. Contains the genotype information 
      from the samples in a VCF file. \code{value} must be a \code{SimpleList} 
      of \code{matrices} or \code{arrays}.

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
      Returns or sets \code{rowData} data. Contains a \code{GRanges} of the 
      first 7 required fields of a VCF file, CHROM, POS, ID, REF, ALT, QUAL and 
      FILTER. REF, ALT, QUAL and FILTER are elementMetadata columns. 
      \code{value} must be a \code{GRanges} with names representing the ID's
      in the VCF file.
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

\section{Arguments}{
  \describe{
    \item{geno}{A \code{list} or \code{SimpleList} of matrix elements,
      or a \code{matrix}. Each element of the list must have the same
      dimensions, and dimension names (if present) must be consistent
      across elements and with the row names of \code{rowData},
      \code{colData}. This is equivalent to the \code{assays} argument
      in  \code{SummarizedExperiment}.
    }
    \item{info}{A \code{SimpleList}. This slot is intended to
      hold the fields from the INFO field of a VCF file. Each element of the
      list must have the same number of rows which must match the number of
      rows in \code{rowData}.
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
      should be applied to extracted assay elements.
    }
    \item{verbose}{A \code{logical(1)} indicating whether messages
      about data coercion during construction should be printed.
    }
    \item{x}{An instance of \code{VCF}-class.
    }
    \item{i, j}{Integer or logical

      For \code{VCF}, \code{VCF<-},
      \code{i}, \code{j} are instances that can act to subset the
      underlying \code{rowData}, \code{colData}, and \code{matrix}
      elements of \code{geno}.
    }
  }
}


\details{
  The \code{VCF} class is designed to hold data from a Variant Call 
  Format (VCF) file. The class extends \code{SummarizedExperiment}
  with the addition of a single slot, \code{info}. For a detailed 
  description of the slots other than \code{info} 
  see ?\code{SummarizedExperiment}.

  Slots unique to the \code{VCF} class,
  \describe{
    \item{\code{info}}{A \link{SimpleList}-class instance containing
      information from the INFO fields in a VCF file.
    }
  }

  Slots inherited from the \code{SummarizedExperiment} class,
  \describe{
    \item{\code{exptData}}{A \link{SimpleList}-class instance containing
      the header or other information about the overall experiment.
    }
    \item{\code{rowData}}{A \link{GRanges}-class instance defining the
      variant ranges and associated metadata columns of REF, ALT, QUAL and
      FILTER.
    }
    \item{\code{colData}}{A \link{DataFrame}-class instance describing
      the samples and associated metadata.
    }
    \item{\code{geno}}{A \link{SimpleList}-class instance, each
      element of which is a matrix or array summarizing data associated with 
      the corresponding range and sample. The \code{geno} accessor is a wrapper
      around the \code{assays} slot in \code{SummarizedExperiment}.
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
  fl <- system.file("extdata", "ex2.vcf", package="VariantAnnotation")
  vcf <- readVcf(fl, "hg19")
 
  names(info(vcf))
  identical(info(vcf)$AF, info(vcf)[[3]])
}
