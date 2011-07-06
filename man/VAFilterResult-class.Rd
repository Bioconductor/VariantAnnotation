\name{VAFilterResult-class}
\docType{class}
\alias{VAFilterResult-class}
\alias{VAFilterResult}
\alias{Logic,VAFilterResult,VAFilterResult-method}
\alias{!,VAFilterResult-method}
\alias{name,VAFilterResult-method}
\alias{show,VAFilterResult-method}
\alias{stats}
\alias{stats,VAFilterResult-method}

\title{"VAFilterResult" for VAFilter output and statistics}

\description{

  Objects of this class are logical vectors indicating records passing
  the applied filter, with an associated data frame summarizing the
  name, input number of records, records passing filter, and logical
  operation used for all filters in which the result participated.
  
}

\usage{
VAFilterResult(data = GRanges(), x = logical(), name = NA_character_, 
   subset = TRUE, input = length(x), passing = sum(x), op = NA_character_)
\S4method{Logic}{VAFilterResult,VAFilterResult}(e1, e2)
\S4method{name}{VAFilterResult}(x, ...)
stats(x, ...)
\S4method{show}{VAFilterResult}(object)
}

\arguments{

  \item{data}{A \code{linkS4class{GRanges}} of the data to be filtered.
  }
  \item{x, object, e1, e2}{For \code{VAFilterResult}, \code{logical()}
    indicating records that passed filter or, for others, an instance of
    \code{VAFilterResult} class.
  }
  \item{name}{\code{character()} indicating the name by which the filter
    is to be referred. Internally, \code{name}, \code{input},
    \code{passing}, and \code{op} may all be vectors representing
    columns of a \code{data.frame} summarizing the application of
    successive filters.
  }
  \item{subset}{\code{logical()} when TRUE, the function returns a 
    subset of the data that passed the filter as a
    \code{linkS4class{GRanges}} object. The \code{name} and \code{stats}
    information are in the metadata slot. When FALSE, an object of 
    \code{VAFilterResult} is returned.
  }
  \item{input}{\code{integer()} indicating the length of the original
    input.
  }
  \item{passing}{\code{integer()} indicating the number of records
    passing the filter.
  }
  \item{op}{\code{character()} indicating the logical operation, if any,
    associated with this filter.
  }
  \item{...}{Additional arguments, unused in methods documented on this
    page.
  }

}

\section{Objects from the Class}{
 
  Objects can be created through \code{\link{VAFilterResult}}, but these
  are automatically created by the application of \code{\link{vaFilter}}
  instances.

}

\section{Slots}{
  \describe{

    \item{\code{.Data}:}{Object of class \code{"logical"} indicating
      records that passed the filter. }

    \item{\code{name}:}{Object of class \code{"ScalarCharacter"}
      representing the name of the filter whose results are
      summarized. The name is either the actual name of the filter, or a
      combination of filter names and logical operations when the
      outcome results from application of several filters in a single
      logical expression. }

    \item{\code{stats}:}{Object of class \code{"data.frame"} summarizing
      the name, input number of records, records passing filter, and
      logical operation used for all filters in which the result
      participated. The \code{data.frame} rows correspond either to
      single filters, or to logical combinations of filters.}

  }
}

\section{Extends}{
Class \code{"\linkS4class{logical}"}, from data part.
Class \code{"\linkS4class{.VAUtil}"}, directly.
Class \code{"\linkS4class{vector}"}, by class "logical", distance 2.
Class \code{"\link{atomic}"}, by class "logical", distance 2.
Class \code{"vectorORfactor"}, by class "logical", distance 3.
}
\section{Methods}{
  \describe{

    \item{Logic}{\code{signature(e1 = "VAFilterResult", e2 =
      "VAFilterResult")}: logic operations on filters.}

    \item{!}{\code{signature(x = "VAFilterResult")}: Negate the outcome
      of the current filter results }

    \item{name}{\code{signature(x = "VAFilterResult")}: The name of the
      filter that the results are based on.}

    \item{stats}{\code{signature(x = "VAFilterResult")}: a
      \code{data.frame} as described in the \sQuote{Slots} section of
      this page.}

    \item{show}{\code{signature(object = "VAFilterResult")}: summary of
      filter results.}

  }

}

\author{Valerie Obenchain}

\seealso{\code{\link{vaFilter}}}

\examples{
## see ?vaFilter
}

\keyword{classes}
