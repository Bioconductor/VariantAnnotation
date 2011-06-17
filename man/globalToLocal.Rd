\name{globalToLocal}
\alias{globalToLocal}

\title{globalToLocal}

\description{
  Converts reference-based (aka chromosome-based or genomic)
  locations into transcript-based locations
}

\usage{globalToLocal(global, ranges)}

\arguments{
  \item{global}{A \code{\linkS4class{GRanges}} object}
  \item{ranges}{A \code{\linkS4class{GRangesList}} object}
}

\details{
  This function translates coordinates from a reference location specified
  in \code{global} to the transcript location as defined in \code{ranges}.
}

\value{
  A \link[IRanges]{DataFrame} with three columns named 'global.ind', ranges.ind'
  and 'local'. The 'global.ind' and 'ranges.ind' columns are the indices
  of the \code{global} and \code{ranges} arguments, respectively. The
  'local' column contains the range coordinates of the elements inf
  \code{global} translated with respect to the ranges in 'ranges.ind'. 
}


\author{Michael Lawrence}


\seealso{
}

\examples{
}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
