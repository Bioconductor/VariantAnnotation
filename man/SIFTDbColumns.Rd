\name{SIFTDbColumns}

\alias{SIFTDbColumns}

\title{SIFTDb Columns}

\description{
  Description of the SIFT Sqlite Database Columns
}

\section{Column descriptions}{
  These column names are displayed when \code{cols} is called on a
  \code{SIFTDb} object.

  \describe{
    \itemize{
      \item rsid : rsid
      \item protein_id : NCBI RefSeq protein ID
      \item aa_change : amino acid substitution; reference aa is preceeding,
                       followed by the position and finally the snp aa
      \item method : method of obtaining related sequences using PSI-BLAST
      \item aa : either the reference or snp residue amino acid
      \item prediction : SIFT prediction
      \item score : SIFT score (range 0 to 1)
        \itemize{
          \item TOLERATED : score is greater than 0.05
          \item DAMAGING : score is less than or equal to 0.05
          \item NOT SCORED : no prediction is made if there are less than 
                             2 homologous sequences that have an amino acid 
                             at the position of the given SNP or if the SIFT 
                             prediction is not available
       } 
      \item median : diversity measurement of the sequences used for 
                     prediction (range 0 to 4.32) 
      \item position_seqs : number of sequences with an amino acide at the 
                           position of prediction
      \item total_seqs : total number of sequences in alignment
    }
  }
}

\seealso{
  \code{?SIFTDb}
}

\author{Valerie Obenchain <vobencha@fhcrc.org>}

\keyword{classes}
\keyword{methods}

