\name{VariantType-class}
\docType{class}

% Class
\alias{VariantType-class}

% Subclasses
\alias{CodingVariants-class}
\alias{IntronVariants-class}
\alias{FiveUTRVariants-class}
\alias{ThreeUTRVariants-class}
\alias{IntergenicVariants-class}
\alias{SpliceSiteVariants-class}
\alias{AllVariants-class}

% Constructors for subclasses :
\alias{CodingVariants}
\alias{IntronVariants}
\alias{FiveUTRVariants}
\alias{ThreeUTRVariants}
\alias{IntergenicVariants}
\alias{SpliceSiteVariants}
\alias{AllVariants}

% show method:
\alias{show,VariantType-method}

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
    AllVariants()
}

\details{
  \code{VariantType} is a virtual class inherited by the \code{CodingVariants},
  \code{IntronVariants}, \code{FiveUTRVariants}, \code{ThreeUTRVariants},
  \code{SpliceSiteVariants}, \code{IntergenicVariants} and \code{AllVariants}
  subclasses.

  The subclasses are used as the \code{region} argument to
  \code{locateVariants}. They designate the type of variant (i.e., region of
  the annotation to match) when calling \code{locateVariants}. The subclasses
  themselves have no slots and require no arguments for an instance to be created.
}

\author{Valerie Obenchain <vobencha@fhcrc.org>}

\examples{
  CodingVariants()
  SpliceSiteVariants()
}
