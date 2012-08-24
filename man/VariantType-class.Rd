\name{VariantType-class}
\docType{class}

% Class
\alias{VariantType-class}

% Subclasses
\alias{CodingVariants-class}
\alias{IntronVariants-class}
\alias{FiveUTRVariants-class}
\alias{ThreeUTRVariants-class}
\alias{SpliceSiteVariants-class}
\alias{IntergenicVariants-class}
\alias{PromoterVariants-class}
\alias{AllVariants-class}

% Constructors for subclasses :
\alias{CodingVariants}
\alias{IntronVariants}
\alias{FiveUTRVariants}
\alias{ThreeUTRVariants}
\alias{SpliceSiteVariants}
\alias{IntergenicVariants}
\alias{PromoterVariants}
\alias{AllVariants}

% accessors :
\alias{upstream}
\alias{upstream,AllVariants-method}
\alias{upstream,PromoterVariants-method}
\alias{upstream<-}
\alias{upstream<-,AllVariants-method}
\alias{upstream<-,PromoterVariants-method}
\alias{downstream}
\alias{downstream,AllVariants-method}
\alias{downstream,PromoterVariants-method}
\alias{downstream<-}
\alias{downstream<-,AllVariants-method}
\alias{downstream<-,PromoterVariants-method}

% show method:
\alias{show,VariantType-method}
\alias{show,AllVariants-method}
\alias{show,PromoterVariants-method}

\title{VariantType subclasses}

\description{
  \code{VariantType} subclasses specify the type of variant to be located with
  \code{locateVariants}.
}

\usage{
    CodingVariants()
    IntronVariants()
    FiveUTRVariants()
    ThreeUTRVariants()
    SpliceSiteVariants()
    IntergenicVariants()
    PromoterVariants(upstream = 2000, downstream = 200)
    AllVariants(upstream = 2000, downstream = 200)
}

\details{
  \code{VariantType} is a virtual class inherited by the \code{CodingVariants},
  \code{IntronVariants}, \code{FiveUTRVariants}, \code{ThreeUTRVariants},
  \code{SpliceSiteVariants}, \code{IntergenicVariants} and \code{AllVariants}
  subclasses.

  The subclasses are used as the \code{region} argument to
  \code{locateVariants}. They designate the type of variant (i.e., region of
  the annotation to match) when calling \code{locateVariants}. 

  The majority of subclasses have no slots and require no arguments for an 
  instance to be created. The two exceptions are \code{PromoterVariants} and 
  \code{AllVariants}. These class have slots for the number of base pairs 
  \code{upstream} of the 5'-end and \code{downstream} of the 3'-end of the
  transcript region.
}

\arguments{
  \item{upstream}{A \code{numeric} indicating the number of base pairs 
    upstream of the 5'-end.
  }
  \item{downstream}{A \code{numeric} indicating the number of base pairs 
    downstream of the 3'-end.
  }
}

\section{Accessors}{
  In the following code, \code{x} is a \code{PromoterVariants} or a
  \code{AllVariants} object.
  \describe{
    \item{}{
      \code{upstream(x)}, \code{upstream(x) <- value}:
      Gets or sets the number of base pairs defining a range
      upstream of the 5' end (excludes 5' start value).
    }
    \item{}{
      \code{downstream(x)}, \code{downstream(x) <- value}:
      Gets or sets the number of base pairs defining a range
      downstream of the 3' end (excludes 3' end value).
    }
  }
}

\author{Valerie Obenchain <vobencha@fhcrc.org>}

\examples{
  CodingVariants()
  SpliceSiteVariants()
  PromoterVariants(upstream=1000, downstream=10000)
}
