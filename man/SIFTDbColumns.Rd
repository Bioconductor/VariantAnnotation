\name{SIFTDbColumns}

\alias{SIFTDbColumns}

\title{SIFTDb Columns}

\description{
  Description of the SIFT Sqlite Database Columns
}

\section{Column descriptions}{
  These column names are displayed when \code{columns} is called on a
  \code{SIFTDb} object.

    \itemize{
      \item RSID : rsid
      \item PROTEINID : NCBI RefSeq protein ID
      \item AACHANGE : amino acid substitution; reference aa is preceeding,
                       followed by the position and finally the snp aa
      \item METHOD : method of obtaining related sequences using PSI-BLAST
      \item AA : either the reference or snp residue amino acid
      \item PREDICTION : SIFT prediction
      \item SCORE : SIFT score (range 0 to 1)
        \itemize{
          \item TOLERATED : score is greater than 0.05
          \item DAMAGING : score is less than or equal to 0.05
          \item NOT SCORED : no prediction is made if there are less than 
                             2 homologous sequences that have an amino acid 
                             at the position of the given SNP or if the SIFT 
                             prediction is not available
       } 
      \item MEDIAN : diversity measurement of the sequences used for 
                     prediction (range 0 to 4.32) 
      \item POSITIONSEQS : number of sequences with an amino acide at the 
                           position of prediction
      \item TOTALSEQS : total number of sequences in alignment
    }
}

\seealso{
  \code{?SIFTDb}
}

\author{Valerie Obenchain <vobencha@fhcrc.org>}

\keyword{classes}
\keyword{methods}

