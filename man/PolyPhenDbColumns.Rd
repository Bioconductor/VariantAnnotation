\name{PolyPhenDbColumns}

\alias{PolyPhenDbColumns}

\title{PolyPhenDb Columns}

\description{
  Description of the PolyPhen Sqlite Database Columns
}

\section{Column descriptions}{
  These column names are displayed when \code{cols} is called on a
  \code{PolyPhenDb} object.

    \itemize{
      \item rsid : rsid
    }
    Original query :
    \itemize{
      \item OSNPID : original SNP identifier from user input
      \item OSNPACC : original protein identifier from user input
      \item OPOS : original substitution position in the protein sequence
                    from user input
      \item OAA1 : original wild type (reference) aa residue from user input 
      \item OAA2 : original mutant (reference) aa residue from user input 
    }
    Mapped query :
    \itemize{
      \item SNPID : SNP identifier mapped to dbSNP rsID if available,
                     otherwise same as o_snp_id. This value was used as the 
                     rsid column 
      \item ACC : protein UniProtKB accession if known protein, otherwise
                  same as o_acc 
      \item POS : substitution position mapped to UniProtKB protein sequence
                  if known, otherwise same as o_pos 
      \item AA1 : wild type aa residue 
      \item AA2 : mutant aa residue 
      \item NT1 : wild type allele nucleotide 
      \item NT2 : mutant allele nucleotide 
    }
    PolyPhen-2 prediction :
    \itemize{
      \item PREDICTION : qualitative ternary classification FPR thresholds 
    }
    PolyPhen-1 prediction :
    \itemize{
      \item BASEDON : prediction basis 
      \item EFFECT : predicted substitution effect on the protein structure or
                     function
    }
    PolyPhen-2 classifiers : 
    \itemize{
      \item PPH2CLASS : binary classifier outcome ("damaging" or "neutral") 
      \item PPH2PROB : probability of the variation being dammaging 
      \item PPH2FPR : false positive rate at the pph2_prob level 
      \item PPH2TPR : true positive rate at the pph2_prob level 
      \item PPH2FDR : false discovery rate at the pph2_prob level 
    }
    UniProtKB-SwissProt derived protein sequence annotations :
    \itemize{
      \item SITE : substitution SITE annotation 
      \item REGION : substitution REGION annotation 
      \item PHAT : PHAT matrix element for substitution in the 
                    TRANSMEM region 
    }
    Multiple sequence alignment scores :
    \itemize{
      \item DSCORE : difference of PSIC scores for two aa variants 
                     (Score1 - Score2) 
      \item SCORE1 : PSIC score for wild type aa residue (aa1) 
      \item SCORE2 : PSIC score for mutant aa residue (aa2) 
      \item NOBS : number of residues observed at the substitution position
                   in the multiple alignment (sans gaps)
    }
    Protein 3D structure features : 
    \itemize{
      \item NSTRUCT : initial number of BLAST hits to similar proteins with
                      3D structures in PDB 
      \item NFILT : number of 3D BLAST hits after identity threshold filtering 
      \item PDBID : protein structure identifier from PDB 
      \item PDBPOS : position of substitution in PDB protein sequence 
      \item PDBCH : PDB polypeptide chain identifier 
      \item IDENT : sequence identity between query and aligned PDB sequences 
      \item LENGTH : PDB sequence alignment length 
      \item NORMACC : normalized accessible surface 
      \item SECSTR : DSSP secondary structure assignment 
      \item MAPREG : region of the phi-psi (Ramachandran) map derived from the
                     residue dihedral angles 
      \item DVOL : change in residue side chain volume 
      \item DPROP : change in solvent accessible surface propensity resulting
                    from the substitution 
      \item BFACT : normalized B-factor (temperature factor) for the residue
      \item HBONDS : number of hydrogen sidechain-sidechain and
                      sidechain-mainchain bonds formed by the residue 
      \item AVENHET : average number of contacts with heteroatoms per residue 
      \item MINDHET : closest contact with heteroatom 
      \item AVENINT : average number of contacts with other chains per residue 
      \item MINDINT : closest contact with other chain 
      \item AVENSIT : average number of contacts with critical sites per
                      residue 
      \item MINDSIT : closest contact with a critical site
    }
    Nucleotide sequence features (CpG/codon/exon junction) :
    \itemize{
      \item TRANSV : whether substitution is a transversion 
      \item CODPOS : position of the substitution within the codon 
      \item CPG : whether or not the substitution changes CpG context
      \item MINDJNC : substitution distance from exon/intron junction 
    }
    Pfam protein family :
    \itemize{
      \item PFAMHIT : Pfam identifier of the query protein
    }
    Substitution scores : 
    \itemize{
      \item IDPMAX : maximum congruency of the mutant aa residue to all 
                     sequences in multiple alignment 
      \item IDPSNP : maximum congruency of the mutant aa residue to the
                     sequence in alignment with the mutant residue 
      \item IDQMIN : query sequence identity with the closest homologue
                     deviating from the wild type aa residue
    }
    Comments : 
    \itemize{
      \item COMMENTS : Optional user comments 
    }
  }

\seealso{
  \code{?PolyPhenDb}
}

\author{Valerie Obenchain <vobencha@fhcrc.org>}

\keyword{classes}
\keyword{methods}

