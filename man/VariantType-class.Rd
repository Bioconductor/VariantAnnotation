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
\alias{upstream,IntergenicVariants-method}
\alias{upstream<-}
\alias{upstream<-,AllVariants-method}
\alias{upstream<-,PromoterVariants-method}
\alias{upstream<-,IntergenicVariants-method}
\alias{downstream}
\alias{downstream,AllVariants-method}
\alias{downstream,PromoterVariants-method}
\alias{downstream,IntergenicVariants-method}
\alias{downstream<-}
\alias{downstream<-,AllVariants-method}
\alias{downstream<-,PromoterVariants-method}
\alias{downstream<-,IntergenicVariants-method}
\alias{promoter}
\alias{promoter,AllVariants-method}
\alias{promoter<-}
\alias{promoter<-,AllVariants-method}
\alias{intergenic}
\alias{intergenic,AllVariants-method}
\alias{intergenic<-}
\alias{intergenic<-,AllVariants-method}

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
    IntergenicVariants(upstream = 100000L, downstream = 100000L)
    PromoterVariants(upstream = 2000L, downstream = 200L)
    AllVariants(promoter = PromoterVariants(),
                intergenic = IntergenicVariants())
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
  instance to be created. \code{PromoterVariants} and \code{IntergenicVariants} 
  and accept \code{upstream} and \code{downstream} arguments that define
  the number of base pairs upstream from the 5'-end and downstream from
  the 3'-end of the transcript region. See the ?\code{locateVariants} man
  page for details. \code{AllVariants} accepts \code{promoter} and
  \code{intergenic} arguments which are \code{PromoterVariants()} and
  \code{IntergenicVariants()} objects with the appropriate 
  \code{upstream} and \code{downstream} values. 
}

\arguments{
  \item{upstream, downstream}{
    Single \code{integer} values representing the number of base pairs 
    upstream of the 5'-end and downstream of the 3'-end. Used in contructing 
    \code{PromoterVariants()} and \code{IntergenicVariants()} objects only.
  }
  \item{promoter}{
    \code{PromoterVariants} object with appropriate \code{upstream} and 
    \code{downstream} values. Used when constructing \code{AllVariants} 
    objects only.
  }
  \item{intergenic}{
    \code{IntergenicVariants} object with appropriate \code{upstream} and 
    \code{downstream} values. Used when constructing \code{AllVariants} 
    objects only.
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
    \item{}{
      \code{promoters(x)}, \code{promoters(x) <- value}:
      Gets or sets the \code{PromoterVariants} in the 
      \code{AllVariants} object.
    }
    \item{}{
      \code{intergenic(x)}, \code{intergenic(x) <- value}:
      Gets or sets the \code{IntergenicVariants} in the 
      \code{AllVariants} object.
    }
  }
}

\author{Valerie Obenchain <vobencha@fhcrc.org>}

\examples{
  CodingVariants()
  SpliceSiteVariants()
  PromoterVariants(upstream=1000, downstream=10000)

  ## Default values for PromoterVariants and IntergenicVariants
  AllVariants()
  ## Modify 'upstream' and 'downstream' for IntergenicVariants
  AllVariants(intergenic=IntergenicVariants(500, 100))
  ## Reset PromoterVariants on existing AllVariants object
  av <- AllVariants()
  av
  promoter(av) <- PromoterVariants(100, 50)
  av
}
