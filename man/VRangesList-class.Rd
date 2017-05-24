\name{VRangesList-class}
\docType{class}

% Class:
\alias{class:VRangesList}
\alias{VRangesList-class}
\alias{CompressedVRangesList-class}
\alias{class:CompressedVRangesList}
\alias{SimpleVRangesList-class}
\alias{class:SimpleVRangesList}

% Constructors:
\alias{VRangesList}

% Accessors:
\alias{alt,VRangesList-method}
\alias{ref,VRangesList-method}

% Aggregation:
\alias{stackSamples}
\alias{stackSamples,VRangesList-method}

\title{VRangesList objects}

\description{
  VRangesList is a virtual class representing a list of
  \code{\linkS4class{VRanges}} objects and should behave much like any
  other derivative of \code{List}. It has both a simple and
  compressed implementation. VRangesList provides conveniences for
  manipulating sets of \code{VRanges} objects.
}

\section{Constructor}{
  \describe{
    \item{}{
      \code{VRangesList(...)}:
      Creates a VRangesList object from \code{VRanges} objects in \dots.
    }
  }
}

\section{Accessors}{
  \describe{
    \item{}{
      \code{alt(x)}: Returns a CharacterList or RleList, effectively by
      calling \code{alt(x[[i]])} on each element of \code{x}.
    }
    \item{}{
      \code{ref(x)}:
      Returns a CharacterList, effectively by calling \code{ref(x[[i]])}
      on each element of \code{x}.
    }
  }
}

\section{Utilities}{
  \describe{
    \item{}{
      \code{stackSamples(x)}: Concentrates the elements in \code{x},
      using \code{names(x)} to appropriately fill \code{sampleNames} in
      the result.
    }
  }
}
 
\author{Michael Lawrence}

\examples{
## construction
example(VRanges)
vrl <- VRangesList(sampleA = vr, sampleB = vr)
stackSamples(vrl)
}
