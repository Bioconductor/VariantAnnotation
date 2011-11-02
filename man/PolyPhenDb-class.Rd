\name{PolyPhenDb-class}
\docType{class}

\alias{PolyPhen}
\alias{PolyPhenDb}
\alias{class:PolyPhenDb}
\alias{PolyPhenDb-class}

\alias{duplicateRSID}
\alias{metadata,PolyPhenDb-method}
\alias{cols,PolyPhenDb-method}
\alias{keys,PolyPhenDb-method}
\alias{select,PolyPhenDb-method}

\title{PolyPhenDb objects}

\description{
  The PolyPhenDb class is a container for storing a connection to a PolyPhen 
  sqlite database.
}

\section{Methods}{
  In the code below, \code{x} is a \code{PolyPhenDb} object.
  \describe{
    \item{}{
      \code{metadata(x)}:
      Returns \code{x}'s metadata in a data frame.
    }
    \item{}{
      \code{cols(x)}:
      Returns the names of the \code{cols} that can be used to subset the
      data columns. For column descriptions see \code{?PolyPhenDbColumns}.
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
      are returned. See \code{?PolyPhenDbColumns} for column descriptions.
    } 
    \item{}{
      \code{duplicateRSID(x)}:
      Returns a named list of duplicate rsid groups. The names are the 
      \code{keys}, the list elements are the rsid's that have been 
      reported as having identical chromosome position and alleles and 
      therefore translating into the same amino acid residue substitution.
    }
  }
}

\details{
  PolyPhen (Polymorphism Phenotyping) is a tool which predicts the possible 
  impact of an amino acid substitution on the structure and function of a 
  human protein by applying empirical rules to the sequence, phylogenetic 
  and structural information characterizing the substitution. 

  PolyPhen makes its predictions using UniProt features, PSIC profiles scores 
  derived from multiple alignment and matches to PDP or PQS structural 
  databases. The procedure can be roughly outlined in the following 
  steps, see the references for complete details,
  \itemize{
    \item sequence-based characterization of substitution site 
    \item calculation of PSIC profile scores for two amino acid variants 
    \item calculation of structural parameters and contacts 
    \item prediction 
  }
  PolyPhen uses empirically derived rules to predict that a non-synonymous 
  SNP is
  \itemize{
    \item probably damaging : it is with high confidence supposed to affect
          protein function or structure
    \item possibly damaging : it is supposed to affect protein function or
          structure
    \item benign : most likely lacking any phenotypic effect
    \item unknown : when in some rare cases, the lack of data do not allow PolyPhen to
          make a prediction 
  }
}
\seealso{
  \code{?PolyPhenDbColumns}
}

\references{
  PolyPhen Home:
  \url{http://genetics.bwh.harvard.edu/pph2/dokuwiki/}

  Adzhubei IA, Schmidt S, Peshkin L, Ramensky VE, Gerasimova A, Bork P,
  Kondrashov AS, Sunyaev SR. Nat Methods 7(4):248-249 (2010).

  Ramensky V, Bork P, Sunyaev S. Human non-synonymous SNPs: server and
  survey. Nucleic Acids Res 30(17):3894-3900 (2002).

  Sunyaev SR, Eisenhaber F, Rodchenkov IV, Eisenhaber B, Tumanyan VG,
  Kuznetsov EN. PSIC: profile extraction from sequence alignments with
  position-specific counts of independent observations. Protein
  Eng 12(5):387-394 (1999).
}

\author{Valerie Obenchain <vobencha@fhcrc.org>}

\examples{
library(PolyPhen.Hsapiens.dbSNP131)

## metadata
metadata(PolyPhen.Hsapiens.dbSNP131)

## available rsid's 
head(keys(PolyPhen.Hsapiens.dbSNP131))

## column descriptions found at ?PolyPhenDbColumns
cols(PolyPhen.Hsapiens.dbSNP131)

## subset on keys and cols 
subst <- c("AA1", "AA2", "PREDICTION")
rsids <- c("rs2142947", "rs4995127", "rs3026284")
select(PolyPhen.Hsapiens.dbSNP131, keys=rsids, cols=subst)

## retrieve substitution scores 
subst <- c("IDPMAX", "IDPSNP", "IDQMIN")
select(PolyPhen.Hsapiens.dbSNP131, keys=rsids, cols=subst)

## retrieve the PolyPhen-2 classifiers 
subst <- c("PPH2CLASS", "PPH2PROB", "PPH2FPR", "PPH2TPR", "PPH2FDR")
select(PolyPhen.Hsapiens.dbSNP131, keys=rsids, cols=subst)

## duplicate groups of rsid's
duplicateRSID(PolyPhen.Hsapiens.dbSNP131, c("rs71225486", "rs1063796"))

}

\keyword{classes}
\keyword{methods}

