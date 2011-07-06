\name{VAFilter-class}
\docType{class}
\alias{VAFilter-class}
\alias{vaFilter,VAFilter-method}
\alias{name}
\alias{name,VAFilter-method}
\alias{show,VAFilter-method}

\alias{c,VAFilter-method}
\alias{coerce,VAFilter,FilterRules-method}

\title{"VAFilter" for representing functions operating on genetic variants}

\description{

  Objects of this class are functions that operate on
  \code{\linkS4class{GRanges}} objects containing genetic variants.
  They return logical vectors indicating which parts of the object satisfy 
  the filter criterion and statistics describing the number of records that 
  passed the filter.

}

\section{Objects from the Class}{

  Objects are created through calls to constructors for predefined
  filters, as described on the \code{\link{vaFilter}} page.

}

\section{Slots}{
  \describe{

    \item{\code{.Data}:}{Object of class \code{"function"} taking a
      single named argument \code{x} corresponding to the 
      \code{\link{GRanges}} object that the filter will be applied to. 
      The return value of the filter function is expected to be a logical 
      vector that can be used to subset \code{x} to include those elements 
      of \code{x} satisfying the filter.
    }

    \item{\code{name}:}{Object of class \code{"ScalarCharacter"}
      representing the name of the filter. The name is useful for
      suggesting the purpose of the filter, and for debugging failed
      filters.
    }
  }
}

\section{Extends}{
Class \code{"\linkS4class{function}"}, from data part.
Class \code{"\linkS4class{.VAUtil}"}, directly.
Class \code{"\linkS4class{OptionalFunction}"}, by class "function", distance 2.
Class \code{"\linkS4class{PossibleMethod}"}, by class "function", distance 2.
}
\section{Methods}{
  \describe{

    \item{vaFilter}{\code{signature(fun = "VAFilter")}: Return the function
      representing the underlying filter; this is primarily for
      interactive use to understanding filter function; usually the
      filter is invoked as a normal function call
    }

    \item{name}{\code{signature(x = "VAFilter")}: Return, as a
      \code{ScalarCharacter}, the name of the function.
    }

    \item{show}{\code{signature(object = "VAFilter")}: display a brief
      summary of the filter
    }

  }
}

\author{Valerie Obenchain <vobencha@fhcrc.org>}

\seealso{
  \code{\link{vaFilter}} for predefined filters.
}

\examples{
## see ?vaFilter
}

\keyword{classes}
