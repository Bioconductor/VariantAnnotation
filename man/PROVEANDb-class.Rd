\name{PROVEANDb-class}
\docType{class}

\alias{PROVEAN}
\alias{PROVEANDb}
\alias{class:PROVEANDb}
\alias{PROVEANDb-class}

\alias{columns,PROVEANDb-method}
\alias{keys,PROVEANDb-method}
\alias{keytypes,PROVEANDb-method}
\alias{select,PROVEANDb-method}

\title{PROVEANDb objects}

\description{
  The PROVEANDb class is a container for storing a connection to a PROVEAN 
  sqlite database.
}

\section{Methods}{
  In the code below, \code{x} is a \code{PROVEANDb} object.
  \describe{
    \item{}{
      \code{metadata(x)}:
      Returns \code{x}'s metadata in a data frame.
    }
    \item{}{
      \code{columns(x)}:
      Returns the names of the \code{columns} that can be used to subset the
      data columns.
    }
    \item{}{
      \code{keys(x, keytype="DBSNPID", ...)}:
      Returns the names of the \code{keys} that can be used to subset the
      data rows. For SIFT.Hsapiens.dbSNP137 the \code{keys} are NCBI dbSNP ids.
    }
    \item{}{
      \code{keytypes(x)}:
      Returns the names of the \code{columns} that can be used as \code{keys}.
      For SIFT.Hsapiens.dbSNP137 the NCBI dbSNP ids are the only keytype.
    }
    \item{}{
      \code{select(x, keys = NULL, columns = NULL, keytype = "DBSNPID", ...)}:
      Returns a subset of data defined by the character vectors \code{keys} 
      and \code{columns}. If no \code{keys} are supplied, all rows are
      returned. If no \code{columns} are supplied, all columns
      are returned.
    }
  }
}

\details{
  The SIFT tool is no longer actively maintained. A few of the
  orginal authors have started the PROVEAN (Protein Variation
  Effect Analyzer) project. PROVEAN is a software tool which predicts
  whether an amino acid substitution or indel has an impact on the
  biological function of a protein. PROVEAN is useful for filtering
  sequence variants to identify nonsynonymous or indel variants that
  are predicted to be functionally important.

  See the web pages for a complete description of the methods.
  \itemize{
    \item PROVEAN Home: \url{http://provean.jcvi.org/index.php/}
    \item SIFT Home: \url{http://sift.jcvi.org/}
  }

  Though SIFT is not under active development, the PROVEAN team still
  provids the SIFT scores in the pre-computed downloads. This package,
  \code{SIFT.Hsapiens.dbSNP137}, contains both SIFT and PROVEAN scores.
  One notable difference between this and the previous SIFT database
  package is that \code{keys} in \code{SIFT.Hsapiens.dbSNP132} are
  rs IDs whereas in \code{SIFT.Hsapiens.dbSNP137} they are NCBI dbSNP IDs.
}

\references{
 The PROVEAN tool has replaced SIFT:
  \url{http://provean.jcvi.org/about.php}

  Choi Y, Sims GE, Murphy S, Miller JR, Chan AP (2012) Predicting the
  Functional Effect of Amino Acid Substitutions and Indels.
  PLoS ONE 7(10): e46688.

  Choi Y (2012) A Fast Computation of Pairwise Sequence Alignment Scores
  Between a Protein and a Set of Single-Locus Variants of Another Protein.
  In Proceedings of the ACM Conference on Bioinformatics,
  Computational Biology and Biomedicine (BCB '12). ACM, New York, NY, USA,
  414-417.

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
#  library(SIFT.Hsapiens.dbSNP137)
#
#  ## metadata
#  metadata(SIFT.Hsapiens.dbSNP137)
#
#  ## keys are the DBSNPID (NCBI dbSNP ID)
#  dbsnp <- keys(SIFT.Hsapiens.dbSNP137)
#  head(dbsnp)
#  columns(SIFT.Hsapiens.dbSNP137)
#
#  ## Return all columns. Note that the key, DBSNPID,
#  ## is always returned. 
#  select(SIFT.Hsapiens.dbSNP137, dbsnp[10])
#  ## subset on keys and cols 
#  cols <- c("VARIANT", "PROVEANPRED", "SIFTPRED")
#  select(SIFT.Hsapiens.dbSNP137, dbsnp[20:23], cols)
}

\keyword{classes}
\keyword{methods}

