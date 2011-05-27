\name{getTranscriptSeqs}
\alias{getTranscriptSeqs}

\title{Get transcript sequences}

\description{ }

\usage{getTranscriptSeqs(cds_by_tx, BSgenomeOrganism)}

\arguments{
  \item{cds_by_tx}{
  }
  \item{BSgenomeOrganism}{
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(cds_by_tx, organism) {
  cds_seqs <- getSeq(organism, unlist(cds_by_tx, use.names=FALSE))
  tx_id <- rep(seq(length(cds_by_tx)), elementLengths(cds_by_tx))
  unlist(lapply(split(cds_seqs, tx_id), paste, collapse=""))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
