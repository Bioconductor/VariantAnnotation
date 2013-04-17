\name{locateVariants}
\alias{locateVariants}

\alias{locateVariants,Ranges,TranscriptDb,VariantType-method}
\alias{locateVariants,Ranges,GRangesList,VariantType-method}
\alias{locateVariants,GRanges,TranscriptDb,VariantType-method}
\alias{locateVariants,GRanges,GRangesList,VariantType-method}
\alias{locateVariants,VCF,TranscriptDb,VariantType-method}
\alias{locateVariants,VCF,GRangesList,VariantType-method}

\alias{locateVariants,GRanges,TranscriptDb,CodingVariants-method}
\alias{locateVariants,GRanges,GRangesList,CodingVariants-method}
\alias{locateVariants,GRanges,TranscriptDb,IntronVariants-method}
\alias{locateVariants,GRanges,GRangesList,IntronVariants-method}
\alias{locateVariants,GRanges,TranscriptDb,FiveUTRVariants-method}
\alias{locateVariants,GRanges,GRangesList,FiveUTRVariants-method}
\alias{locateVariants,GRanges,TranscriptDb,ThreeUTRVariants-method}
\alias{locateVariants,GRanges,GRangesList,ThreeUTRVariants-method}
\alias{locateVariants,GRanges,TranscriptDb,IntergenicVariants-method}
\alias{locateVariants,GRanges,GRangesList,IntergenicVariants-method}
\alias{locateVariants,GRanges,TranscriptDb,SpliceSiteVariants-method}
\alias{locateVariants,GRanges,GRangesList,SpliceSiteVariants-method}
\alias{locateVariants,GRanges,TranscriptDb,PromoterVariants-method}
\alias{locateVariants,GRanges,GRangesList,PromoterVariants-method}
\alias{locateVariants,GRanges,TranscriptDb,AllVariants-method}
\alias{locateVariants,GRanges,GRangesList,AllVariants-method}

\title{Locate variants}

\description{Variant location with respect to gene function}

\usage{
locateVariants(query, subject, region, ...)
\S4method{locateVariants}{VCF,TranscriptDb,VariantType}(query, subject, region, ..., cache=new.env(parent=emptyenv()), ignore.strand=FALSE, asHits=FALSE)
\S4method{locateVariants}{GRanges,TranscriptDb,VariantType}(query, subject, region, ..., cache=new.env(parent=emptyenv()), ignore.strand=FALSE, asHits=FALSE)
}

\arguments{
  \item{query}{A \link[IRanges]{Ranges}, \link[GenomicRanges]{GRanges} or 
    \linkS4class{VCF} object containing the variants. Metadata columns are 
    allowed but are ignored.

    NOTE: Circular chromosomes are not currently supported. 
  }
  \item{subject}{A \link[GenomicFeatures]{TranscriptDb} or \code{GRangesList}
    object that serves as the annotation. GFF files can be converted to
    \link[GenomicFeatures]{TranscriptDb} objects with 
    \code{makeTranscriptDbFromGFF()} in the \code{GenomicFeatures} package.
  }
  \item{region}{An instance of one of the 8 VariantType classes :
    \code{CodingVariants}, \code{IntronVariants}, \code{FiveUTRVariants}, 
    \code{ThreeUTRVariants}, \code{IntergenicVariants}, 
    \code{SpliceSiteVariants}, \code{PromoterVariants}, \code{AllVariants}.
    All objects can be instantiated with no arguments, e.g., CodingVariants()
    will create an object of \code{CodingVariants}. 

    \code{AllVariants}, \code{PromoterVariants} and \code{IntergenicVariants}
    have \code{upstream} and \code{downstream} arguments. For
    \code{PromoterVariants} and \code{IntergenicVariants} these are single
    integer values >= 0. For \code{AllVariants} these are integer vectors
    of length 2 named \sQuote{promoter} and \sQuote{intergenic}. See
    ?\code{upstream} for more details.

    When using \code{AllVariants}, a range in \code{query} may fall in
    multiple regions (e.g., 'intergenic' and 'promoter'). In this case
    the result will have a row for each match. All data in the
    row will be equivalent except the LOCATION column. 
  }
  \item{\dots}{Additional arguments passed to methods
  }
  \item{cache}{An \code{environment} into which required components
    of \code{subject} are loaded. Provide, and re-use, a cache to
    speed repeated queries to the same \code{subject} across
    different \code{query} instances.
  }
  \item{ignore.strand}{A \code{logical} indicating if strand should be
    igored when performing overlaps.
  }
  \item{asHits}{A \code{logical} indicating if the results should be 
    returned as a \link[IRanges]{Hits} object. Not applicable when
    \code{region} is AllVariants or IntergenicVariants.
  }
}

\details{
  \describe{
    \item{Range representation :}{
      The ranges in \code{query} should reflect the position(s) of the reference 
      allele. For snps the range will be of width 1. For range insertions or 
      deletions the reference allele could be a sequence such as GGTG in which 
      case the width of the range should be 4. 
    }
    \item{Location :}{
      Possible locations are \sQuote{coding}, \sQuote{intron}, 
      \sQuote{threeUTR}, \sQuote{fiveUTR}, \sQuote{intergenic},
      \sQuote{spliceSite}, or \sQuote{promoter}. 

      Overlap operations for \sQuote{coding}, \sQuote{intron}, 
      \sQuote{threeUTR}, and \sQuote{fiveUTR} require variants to fall 
      completely within the defined region to be classified as such. 

      To be classified as a \sQuote{spliceSite} the variant must overlap 
      with any portion of the first 2 or last 2 nucleotides in an intron. 

      \sQuote{intergenic} variants do not fall within a gene region as 
      defined by the annotaion. \sQuote{transcripts by gene} are extracted
      from the annotation and overlapped with the variant ranges. The 
      variants that do not have a \sQuote{hit} are considered 
      \code{intergenic}. When available, gene IDs for the flanking genes 
      are provided as \code{PRECEDEID} and \code{FOLLOWID}. The 
      \code{upstream} and \code{downstream} arguments to 
      \code{IntergenicVariants()} and \code{AllVariants()} control the 
      distance from the query for the flanking genes.

      \sQuote{promoter} variants fall within a specified range upstream and 
      downstream of the transcription start site. Ranges values can be set 
      with the \code{upstream} and \code{downstream} arguments when creating
      the \code{PromoterVariants()} or \code{AllVariants()} classes.
    }
    \item{Subject as GRangesList :}{
      The \code{subject} can be a \code{TranscriptDb} or \code{GRangesList}
      object. When using a \code{GRangesList} the type of data required
      is driven by the \code{VariantType} class. Below is a description of 
      the appropriate \code{GRangesList} for each \code{VariantType}. 
      \describe{
        \item{CodingVariants :}{coding (CDS) by transcript}
        \item{IntronVariants :}{introns by transcript}
        \item{FiveUTRVariants :}{five prime UTR by transcript}
        \item{ThreeUTRVariants :}{three prime UTR by transcript}
        \item{IntergenicVariants :}{transcripts by gene}
        \item{SpliceSiteVariants :}{introns by transcript}
        \item{PromoterVariants :}{list of transcripts}
        \item{AllVariants :}{no GRangeList method available}
      }
    }
    \item{Using the cache :}{
      When processing multiple VCF files performance is enhanced by specifying
      an enviornment as the \code{cache} argument. This cache is used to store 
      and reuse extracted components of the subject (TxDb) required by the 
      function. The first call to the function (i.e., processing the first
      VCF file in a list of many) populates the cache; repeated calls
      to \code{locateVariants} will access these objects from the cache vs
      re-extracting the same information.
    }
  }
}
\value{
  A \code{GRanges} object with a row for each variant-transcript match.
  Columns are \code{LOCATION}, \code{QUERYID}, \code{TXID}, \code{GENEID},
  \code{PRECEDEID}, \code{FOLLOWID} and \code{CDSID}. Each column is
  described in detail below. 

 \describe{
    \item{\code{LOCATION}}{
      Possible locations are \sQuote{coding}, \sQuote{intron}, 
      \sQuote{threeUTR}, \sQuote{fiveUTR}, \sQuote{intergenic},
      \sQuote{spliceSite} and \sQuote{promoter}.

      To be classified as \sQuote{coding}, \sQuote{intron}, \sQuote{threeUTR} 
      or \sQuote{fiveUTR} the variant must fall completely within the region.

      \sQuote{intergenic} variants do not fall within a transcript. The 
      \sQuote{GENEID} for these position are \code{NA} and instead the 
      \sQuote{PRECEDEID} and \sQuote{FOLLOWID} for the flanking genes are 
      given.

      A \sQuote{spliceSite} variant overlaps any portion of the first 2 or last 
      2 nucleotides of an intron.
    }
    \item{\code{QUERYID}}{
      The \code{QUERYID} column provides a map back to the row in the
      original \code{query}. If the \code{query} was a \code{VCF} object this
      index corresponds to the row in the \code{GRanges} in the \code{rowData}
      slot. 
    }
    \item{\code{TXID}}{
      The transcript id taken from the \code{TranscriptDb} object.
    }
    \item{\code{CDSID}}{
      The coding sequence id taken from the \code{TranscriptDb} object.
    }
    \item{\code{GENEID}}{
      The gene id taken from the \code{TranscriptDb} object.
    }
    \item{\code{PRECEDEID}}{
      The id of the gene the query precedes. Only present for
      \sQuote{intergenic} variants. \code{upstream} and
      and \code{downstream} argments control the distance
      of \code{PRECEDEID} and \code{FOLLOWID} from the query. 
    }
    \item{\code{FOLLOWID}}{
      The id of the gene the query follows. Only present for 
      \sQuote{intergenic} variants. \code{upstream} and
      and \code{downstream} argments control the distance
      of \code{PRECEDEID} and \code{FOLLOWID} from the query. 
    }
    All ID values will be \sQuote{NA} for variants with a location of
    \code{transcript_region} or \code{NA}.
  } 
}

\author{Valerie Obenchain <vobencha@fhcrc.org>}

\seealso{
  \code{\link{readVcf}},
  \code{\link{summarizeVariants}},
  \code{\link{predictCoding}}
  \code{\link{upstream}}
}

\examples{
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene 

  ## ---------------------------------------------------------------------
  ## TranscriptDb as subject 
  ## ---------------------------------------------------------------------
  ## Read variants from a VCF file 
  fl <- system.file("extdata", "chr7-sub.vcf.gz", 
                    package="VariantAnnotation")
  vcf <- readVcf(fl, "hg19")

  ## Variant seqlevels do not match the TxDb seqlevels
  head(seqlevels(vcf))
  head(seqlevels(txdb))
  intersect(seqlevels(vcf), seqlevels(txdb))

  ## Use the renameSeqlevels helper to rename seqlevels
  newnames <- paste("chr", seqlevels(vcf), sep="")
  names(newnames) <- seqlevels(vcf)
  vcf <- renameSeqlevels(vcf, newnames) 

  ## Confirm
  intersect(seqlevels(vcf), seqlevels(txdb))

  ## Overlaps for all possible variant locations
  loc_all <- locateVariants(vcf, txdb, AllVariants())
  table(loc_all$LOCATION)

  ## Overlaps for intergenic variants
  region <- IntergenicVariants(upstream=700000, downstream=500000)
  loc_int <- locateVariants(vcf, txdb, region)

  ## Intergenic variants do not overlap a gene range in the 
  ## annotation therefore no GENEID is reported. Instead, the
  ## flanking gene ids are reported as PRECEDEID and FOLLOWID.
  mcols(loc_int)[c("LOCATION", "PRECEDEID", "FOLLOWID")] 

  ## Compute distances from the variant to the flanking genes.
  ## See the ?'distance,GenomicRanges,TranscriptDb-method'
  ## man page in the GenomicRanges package for details on 
  ## computing distance with a GRanges and TranscriptDb. 
  p_dist <- distance(loc_int, txdb, id=loc_int$PRECEDEID, type="gene")
  f_dist <- distance(loc_int, txdb, id=loc_int$FOLLOWID, type="gene")
  DataFrame(p_dist, f_dist)

  ## ---------------------------------------------------------------------
  ## GRangesList as subject 
  ## ---------------------------------------------------------------------
  ## When 'subject' is a GRangesList the GENEID is unavailable and
  ## will always be reported as NA. This is because the GRangesList 
  ## objects are extractions of region-by-transcript, not region-by-gene.
  \dontrun{
  cdsbytx <- cdsBy(txdb)
  locateVariants(vcf, cdsbytx, CodingVariants())

  intbytx <- intronsByTranscript(txdb)
  locateVariants(vcf, intbytx, IntronVariants())
  }

  ## ---------------------------------------------------------------------
  ## Using the cache 
  ## ---------------------------------------------------------------------
  ## When processing multiple VCF files, the 'cache' can be used
  ## to store and reuse extracted components of the TranscriptDb
  ## thereby enhancing performance.
  \dontrun{
  myenv <- new.env()
  files <- list(vcf1, vcf2, vcf3)
  lapply(files, 
      function(fl) {
          vcf <- readVcf(fl, "hg19")
          ## modify seqlevels to match TxDb
          oldnms <- seqlevels(vcf)
          newnms <- paste("chr", oldnms, sep="")
          names(newnms) <- oldnms
          vcf_mod <- renameSeqlevels(vcf, newnms)
          locateVariants(vcf_mod, txdb, AllVariants(), cache=myenv)
      }) 
  }
}

\keyword{methods}
