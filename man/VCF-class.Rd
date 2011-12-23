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
      \code{VCF(assays = SimpleList(), rowData = GRanges(), 
                colData = DataFrame(), exptData = SimpleList(), 
                info = SimpleList(), ..., verbose = FALSE)
      }
    }
  }
}

\section{Accessors}{
  In the following code snippets, \code{x} is a VCF object.

  \describe{
    \item{}{
      \code{info(x, i, ..., withDimnames=TRUE)}
    }
    \item{}{
      \code{info(x, i, ...) <- value}
    }
    \item{}{
      \code{geno(x, i, ..., withDimnames=TRUE)}
    }
    \item{}{
      \code{geno(x, i, ...) <- value}
    }
    \item{}{
      \code{assays(x, ..., withDimnames=TRUE)}
    }
    \item{}{
      \code{assays(x, ...) <- value}
    }
    \item{}{
      \code{assay(x, i, ...)}
    }
    \item{}{
      \code{assay(x, i, ...) <- value}
    }
    \item{}{
      \code{exptData(x, ...)}
    }
    \item{}{
      \code{exptData(x, ...) <- value}
    }
    \item{}{
      \code{rowData(x, ...)}
    }
    \item{}{
      \code{rowData(x, ...) <- value}
    }
    \item{}{
      \code{colData(x, ...)}
    }
    \item{}{
      \code{colData(x, ...) <- value}
    }
  }
}

\section{Subsetting}{
  \describe{
    \item{}{
      \code{x[i, j]}, \code{x[i, j] <- value}:
      Gets or sets rows \code{i} and columns \code{j}.
    }
  }
}

\section{Arguments}{
  \describe{
    \item{assays}{A \code{list} or \code{SimpleList} of matrix elements,
      or a \code{matrix}. Each element of the list must have the same
      dimensions, and dimension names (if present) must be consistent
      across elements and with the row names of \code{rowData},
      \code{colData}.
    }
    \item{info}{An optional \code{SimpleList}. This slot is intended to
      hold the fields from the INFO field of a VCF file. If present, the
      number of fields must match those found in the header stored in
      \code{exptData(x)[["HEADER"]]$INFO}.
    }
    \item{rowData}{A \code{GRanges} or \code{GRangesList} instance describing 
      the ranges of interest. Row names, if present, become the row names of the
      \code{VCF}. The length of the \code{GRanges} or the
      \code{GRangesList} must equal the number of rows of the 
      matrices in \code{assays}.
    }
    \item{colData}{An optional \code{DataFrame} describing the
      samples. Row names, if present, become the column names of the
      \code{VCF}.
    }
    \item{exptData}{An optional \code{SimpleList} of arbitrary content
      describing the overall experiment.
    }
    \item{...}{For \code{VCF}, S4 methods \code{list}
      and \code{matrix}, arguments identical to those of the
      \code{SimpleList} method.
    }
    \item{withDimnames}{A \code{logical(1)}, indicating whether dimnames
      should be applied to extracted assay elements.
    }
    \item{verbose}{A \code{logical(1)} indicating whether messages
      about data coercion during construction should be printed.
    }
    \item{x}{An instance of \code{VCF}-class.
    }
    \item{i, j}{For \code{assay}, \code{assay<-}, \code{i} is a integer or
      numeric scalar; see \sQuote{Details} for additional constraints.

      For \code{VCF}, \code{VCF<-},
      \code{i}, \code{j} are instances that can act to subset the
      underlying \code{rowData}, \code{colData}, and \code{matrix}
      elements of \code{assays}.
    }
    \item{drop}{A \code{logical(1)}, ignored by these methods.
    }
    \item{value}{An instance of a class specified in the S4 method
      signature or as outlined in \sQuote{Details}.
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
      information about the overall experiment.
    }
    \item{\code{rowData}}{A \link{GRanges}-class instance defining the
      ranges of interest and associated metadata.
    }
    \item{\code{colData}}{A \link{DataFrame}-class instance describing
      the samples and associated metadata.
    }
    \item{\code{assays}}{A \link{SimpleList}-class instance, each
      element of which is a matrix summarizing data associated with the
      corresponding range and sample.
    }
  }

  The \code{readVcf} function parses the genotype information from a VCF file 
  into the \code{assays} slot. For this reason a \code{geno} accessor was 
  created for the \code{VCF} class for convenience. It is equivalent in 
  function to the \code{assays} accessor from \code{SummarizedExperiment}.
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
  
  showClass(class(vcf))
}
