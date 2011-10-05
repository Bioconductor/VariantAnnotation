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
      \item o_snp_id : original SNP identifier from user input
      \item o_acc : original protein identifier from user input
      \item o_pos : original substitution position in the protein sequence
                    from user input
      \item o_aa1 : original wild type (reference) aa residue from user input 
      \item o_aa2 : original mutant (reference) aa residue from user input 
    }
    Mapped query :
    \itemize{
      \item snp_id : SNP identifier mapped to dbSNP rsID if available,
                     otherwise same as o_snp_id. This value was used as the 
                     rsid column 
      \item acc : protein UniProtKB accession if known protein, otherwise
                  same as o_acc 
      \item pos : substitution position mapped to UniProtKB protein sequence
                  if known, otherwise same as o_pos 
      \item aa1 : wild type aa residue 
      \item aa2 : mutant aa residue 
      \item nt1 : wild type allele nucleotide 
      \item nt2 : mutant allele nucleotide 
    }
    PolyPhen-2 prediction :
    \itemize{
      \item prediction : qualitative ternary classification FPR thresholds 
    }
    PolyPhen-1 prediction :
    \itemize{
      \item based_on : prediction basis 
      \item effect : predicted substitution effect on the protein structure or
                     function
    }
    PolyPhen-2 classifiers : 
    \itemize{
      \item pph2_class : binary classifier outcome ("damaging" or "neutral") 
      \item pph2_prob : probability of the variation being dammaging 
      \item pph2_FPR : false positive rate at the pph2_prob level 
      \item pph2_TPR : true positive rate at the pph2_prob level 
      \item pph2_FDR : false discovery rate at the pph2_prob level 
    }
    UniProtKB-SwissProt derived protein sequence annotations :
    \itemize{
      \item site : substitution SITE annotation 
      \item region : substitution REGION annotation 
      \item PHAT : PHAT matrix element for substitution in the 
                    TRANSMEM region 
    }
    Multiple sequence alignment scores :
    \itemize{
      \item dScore : difference of PSIC scores for two aa variants 
                     (Score1 - Score2) 
      \item Score1 : PSIC score for wild type aa residue (aa1) 
      \item Score2 : PSIC score for mutant aa residue (aa2) 
      \item Nobs : number of residues observed at the substitution position
                   in the multiple alignment (sans gaps)
    }
    Protein 3D structure features : 
    \itemize{
      \item Nstruct : initial number of BLAST hits to similar proteins with
                      3D structures in PDB 
      \item Nfilt : number of 3D BLAST hits after identity threshold filtering 
      \item PDB_id : protein structure identifier from PDB 
      \item PDB_pos : position of substitution in PDB protein sequence 
      \item PDB_ch : PDB polypeptide chain identifier 
      \item ident : sequence identity between query and aligned PDB sequences 
      \item length : PDB sequence alignment length 
      \item NormAcc : normalized accessible surface 
      \item SecStr : DSSP secondary structure assignment 
      \item MapReg : region of the phi-psi (Ramachandran) map derived from the
                     residue dihedral angles 
      \item dVol : change in residue side chain volume 
      \item dProp : change in solvent accessible surface propensity resulting
                    from the substitution 
      \item B_fact : normalized B-factor (temperature factor) for the residue
      \item H_bonds : number of hydrogen sidechain-sidechain and
                      sidechain-mainchain bonds formed by the residue 
      \item AveNHet : average number of contacts with heteroatoms per residue 
      \item MinDHet : closest contact with heteroatom 
      \item AveNInt : average number of contacts with other chains per residue 
      \item MinDInt : closest contact with other chain 
      \item AveNSit : average number of contacts with critical sites per
                      residue 
      \item MinDSit : closest contact with a critical site
    }
    Nucleotide sequence features (CpG/codon/exon junction) :
    \itemize{
      \item Transv : whether substitution is a transversion 
      \item CodPos : position of the substitution within the codon 
      \item CpG : whether or not the substitution changes CpG context
      \item MinDJnc : substitution distance from exon/intron junction 
    }
    Pfam protein family :
    \itemize{
      \item PfamHit : Pfam identifier of the query protein
    }
    Substitution scores : 
    \itemize{
      \item IdPmax : maximum congruency of the mutant aa residue to all 
                     sequences in multiple alignment 
      \item IdPSNP : maximum congruency of the mutant aa residue to the
                     sequence in alignment with the mutant residue 
      \item IdQmin : query sequence identity with the closest homologue
                     deviating from the wild type aa residue
    }
    Comments : 
    \itemize{
      \item Comments : Optional user comments 
    }
  }

\seealso{
  \code{?PolyPhenDb}
}

\author{Valerie Obenchain <vobencha@fhcrc.org>}

\keyword{classes}
\keyword{methods}

