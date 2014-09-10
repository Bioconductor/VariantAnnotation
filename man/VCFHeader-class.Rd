\name{VCFHeader-class}
\docType{class}

% Class
\alias{VCFHeader-class}

% Constructor-like function:
\alias{VCFHeader}

% Accessor methods:
\alias{reference}
\alias{reference,VCFHeader-method}
\alias{samples}
\alias{samples,VCFHeader-method}
\alias{header}
\alias{header,VCFHeader-method}

\alias{contig}
\alias{contig,VCFHeader-method}
\alias{meta}
\alias{meta,VCFHeader-method}
\alias{meta<-}
\alias{meta<-,VCFHeader,DataFrame-method}
\alias{fixed,VCFHeader-method}
\alias{fixed<-,VCFHeader,DataFrameList-method}
\alias{info,VCFHeader-method}
\alias{info<-,VCFHeader,DataFrame-method}
\alias{geno,VCFHeader-method}
\alias{geno<-,VCFHeader,missing,DataFrame-method}
\alias{seqinfo,VCFHeader-method}

% show method:
\alias{show,VCFHeader-method}

\title{VCFHeader instances}

\description{
  The \code{VCFHeader} class holds Variant Call Format (VCF) file 
  header information and is produced from a call to \code{scanVcfHeader}.
}

\section{Constructor}{
  \describe{
    \item{}{
      \code{VCFHeader(reference = character(), samples = character(), 
                      header = DataFrameList(), ...)
      }
    }
  }
}

\section{Accessors}{
  In the following code snippets \code{x} is a VCFHeader object.
  \describe{
    \item{}{
      \code{samples(x)}:
      Returns a character() vector of names of samples. 
    }
    \item{}{
      \code{header(x)}:
      Returns all information in the header slot which includes
      \code{meta}, \code{info} and \code{geno} if present.
    }
    \item{}{
      \code{meta(x)}, \code{meta(x)<- value}:
      Gets or sets a \code{DataFrame} of meta information. This includes 
      any information represented as simple key-value pairs in the
      VCF file header. Replacement value must be a \code{DataFrame}
      with a single column named \sQuote{Value}.
    }
    \item{}{
      \code{fixed(x), fixed(x)<- value}:
      Returns a \code{DataFrameList} of information pertaining to any of
      \sQuote{REF}, \sQuote{ALT}, \sQuote{FILTER} and \sQuote{QUAL}. 
      Replacement value must be a \code{DataFrameList} with names 
      \sQuote{QUAL} or \sQuote{FILTER}. 
    }
    \item{}{
      \code{info(x)}, \code{info(x)<- value}:
      Gets or sets a \code{DataFrame} of \sQuote{INFO} information.
      Replacement value must be a \code{DataFrame} with 3 columns
      named \sQuote{Number}, \sQuote{Type} and \sQuote{Description}. 
    }
    \item{}{
      \code{geno(x)}, \code{geno(x)<- value}:
      Returns a \code{DataFrame} of \sQuote{FORMAT} information.
      Replacement value must be a \code{DataFrame} with 3 columns
      named \sQuote{Number}, \sQuote{Type} and \sQuote{Description}. 
    }
    \item{}{
      \code{reference(x)}:
      Returns a character() vector with names of reference sequences.
      Not relevant for \code{scanVcfHeader}.
    }
  }
}

\section{Arguments}{
  \describe{
    \item{reference}{A character() vector of sequences. 
    }
    \item{sample}{A character() vector of sample names. 
    }
    \item{header}{A \code{DataFrameList} of parsed header lines
    (preceeded by \dQuote{##}) present in the VCF file.
    }
    \item{...}{Additional arguments passed to methods. 
    }
  }
}

\details{
  The \code{VCFHeader} class is holds header information from a
  from VCF file. 

  Slots :
  \describe{
    \item{\code{reference}}{character() vector
    }
    \item{\code{sample}}{character() vector 
    }
    \item{\code{header}}{\link{DataFrameList}-class
    }
  }
}

\author{Valerie Obenchain <vobencha@fhcrc.org>}

\seealso{
  \code{\link{scanVcfHeader}},
  \code{\link{DataFrameList}}
}

\examples{
  fl <- system.file("extdata", "structural.vcf", package="VariantAnnotation")
  hdr <- scanVcfHeader(fl)

  fixed(hdr)
  info(hdr)
  geno(hdr)
}
