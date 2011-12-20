\name{SIFTDb-class}
\docType{class}

\alias{SIFT}
\alias{SIFTDb}
\alias{class:SIFTDb}
\alias{SIFTDb-class}

\alias{metadata,SIFTDb-method}
\alias{cols,SIFTDb-method}
\alias{keys,SIFTDb-method}
\alias{select,SIFTDb-method}

\title{SIFTDb objects}

\description{
  The SIFTDb class is a container for storing a connection to a SIFT 
  sqlite database.
}

\section{Methods}{
  In the code below, \code{x} is a \code{SIFTDb} object.
  \describe{
    \item{}{
      \code{metadata(x)}:
      Returns \code{x}'s metadata in a data frame.
    }
    \item{}{
      \code{cols(x)}:
      Returns the names of the \code{cols} that can be used to subset the
      data columns.
    }
    \item{}{
      \code{keys(x)}:
      Returns the names of the \code{keys} that can be used to subset the
      data rows. The \code{keys} values are the rsid's.
    }
    \item{}{
      \code{select(x, keys = NULL, cols = NULL, ...)}:
      Returns a subset of data defined by the character vectors \code{keys} 
      and \code{cols}. If no \code{keys} are supplied, all rows are
      returned. If no \code{cols} are supplied, all columns
      are returned. For column descriptions see \code{?SIFTDbColumns}.
    }
  }
}


\details{
  SIFT is a sequence homology-based tool that sorts intolerant from tolerant 
  amino acid substitutions and predicts whether an amino acid substitution 
  in a protein will have a phenotypic effect. SIFT is based on the premise 
  that protein evolution is correlated with protein function. Positions 
  important for function should be conserved in an alignment of the protein 
  family, whereas unimportant positions should appear diverse in an alignment.

  SIFT uses multiple alignment information to predict tolerated 
  and deleterious substitutions for every position of the query sequence. 
  The procedure can be outlined in the following steps, 
  \itemize{
    \item search for similar sequences
    \item choose closely related sequences that may share similar
          function to the query sequence
    \item obtain the alignment of the chosen sequences
    \item calculate normalized probabilities for all possible
          substitutions from the alignment.
  } Positions with normalized probabilities less than 0.05 are predicted
  to be deleterious, those greater than or equal to 0.05 are predicted to be
  tolerated.
}

\references{
  SIFT Home:
  \url{http://sift.jcvi.org/}

  Kumar P, Henikoff S, Ng PC. Predicting the effects of coding non-synonymous
  variants on protein function using the SIFT algorithm. Nat Protoc.
  2009;4(7):1073-81

  Ng PC, Henikoff S. Predicting the Effects of Amino Acid Substitutions on
  Protein Function Annu Rev Genomics Hum Genet. 2006;7:61-80.

  Ng PC, Henikoff S. SIFT: predicting amino acid changes that affect protein
  function. Nucleic Acids Res. 2003 Jul 1;31(13):3812-4.
}

\author{Valerie Obenchain <vobencha@fhcrc.org>}

\examples{
library(SIFT.Hsapiens.dbSNP132)

## metadata
metadata(SIFT.Hsapiens.dbSNP132)

## available rsid's 
head(keys(SIFT.Hsapiens.dbSNP132))

## for column descriptions see ?SIFTDbColumns
cols(SIFT.Hsapiens.dbSNP132)

## subset on keys and cols 
rsids <- c("rs2142947", "rs17970171", "rs8692231", "rs3026284") 
subst <- c("RSID", "PREDICTION", "SCORE")
select(SIFT.Hsapiens.dbSNP132, keys=rsids, cols=subst)
select(SIFT.Hsapiens.dbSNP132, keys=rsids[1:2])
}

\keyword{classes}
\keyword{methods}

