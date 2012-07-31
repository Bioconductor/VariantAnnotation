\name{locateVariants}
\alias{locateVariants}

\alias{locateVariants,Ranges,TranscriptDb,ANY-method}
\alias{locateVariants,Ranges,GRangesList,ANY-method}
\alias{locateVariants,GRanges,TranscriptDb,ANY-method}
\alias{locateVariants,GRanges,GRangesList,ANY-method}
\alias{locateVariants,VCF,TranscriptDb,ANY-method}
\alias{locateVariants,VCF,GRangesList,ANY-method}

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
\alias{locateVariants,GRanges,TranscriptDb,AllVariants-method}
\alias{locateVariants,GRanges,GRangesList,AllVariants-method}

\title{Locate variants}

\description{Variant location with respect to gene function}

\usage{
locateVariants(query, subject, region, ...)
\S4method{locateVariants}{VCF,TranscriptDb,ANY}(query, subject, region, ..., cache = new.env(parent=emptyenv()))
\S4method{locateVariants}{Ranges,TranscriptDb,ANY}(query, subject, region, ..., cache = new.env(parent=emptyenv()))
\S4method{locateVariants}{GRanges,TranscriptDb,ANY}(query, subject, region, ..., cache = new.env(parent=emptyenv()))
\S4method{locateVariants}{GRanges,TranscriptDb,CodingVariants}(query, subject, region, ..., cache = new.env(parent=emptyenv()), ignore.strand = FALSE)
\S4method{locateVariants}{GRanges,TranscriptDb,IntronVariants}(query, subject, region, ..., cache = new.env(parent=emptyenv()), ignore.strand = FALSE)
\S4method{locateVariants}{GRanges,TranscriptDb,FiveUTRVariants}(query, subject, region, ..., cache = new.env(parent=emptyenv()), ignore.strand = FALSE)
\S4method{locateVariants}{GRanges,TranscriptDb,ThreeUTRVariants}(query, subject, region, ..., cache = new.env(parent=emptyenv()), ignore.strand = FALSE)
\S4method{locateVariants}{GRanges,TranscriptDb,IntergenicVariants}(query, subject, region, ..., cache = new.env(parent=emptyenv()), ignore.strand = FALSE)
\S4method{locateVariants}{GRanges,TranscriptDb,SpliceSiteVariants}(query, subject, region, ..., cache = new.env(parent=emptyenv()), ignore.strand = FALSE)
\S4method{locateVariants}{GRanges,TranscriptDb,AllVariants}(query, subject, region, ..., cache = new.env(parent=emptyenv()), ignore.strand = FALSE)
}

\arguments{
  \item{query}{A \link[IRanges]{Ranges}, \link[GenomicRanges]{GRanges} or 
    \linkS4class{VCF} object containing the variants. Metadata columns are 
    allowed but are ignored. 
  }
  \item{subject}{A \link[GenomicFeatures]{TranscriptDb} or \code{GRangesList}
    object that serves as the annotation. GFF files can be converted to
    \link[GenomicFeatures]{TranscriptDb} objects with 
    \code{makeTranscriptDbFromGFF()} in the \code{GenomicFeatures} package.
  }
  \item{region}{An instance of one of \code{CodingVariants},
    \code{IntronVariants}, \code{FiveUTRVariants}, \code{ThreeUTRVariants},
    \code{IntergenicVariants}, \code{SpliceSiteVariants}, \code{AllVariants}.
    These objects are instantiated with no arguments, e.g., CodingVariants()
    will create an object of \code{CodingVariants}. 
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
      \sQuote{threeUTR}, \sQuote{fiveUTR}, \sQuote{intergenic}, or
      \sQuote{spliceSite}. 

      Overlap operations for \sQuote{coding}, \sQuote{intron}, 
      \sQuote{threeUTR}, and \sQuote{fiveUTR} require variants to fall 
      completely within the defined region to be classified as such. 

      To be classified as a \sQuote{spliceSite} the variant must overlap with 
      any portion of the first 2 or last 2 nucleotides in an intron. 

      \sQuote{intergenic} variants are those that do not fall within a 
      transcript associated with a gene. If available, gene IDs for the flanking 
      genes are give as the \code{PRECEDEID} and \code{FOLLOWID}.
    }
    \item{Subject as GRangesList :}{
      The \code{subject} can be a \code{TranscriptDb} or \code{GRangesList}
      object. When using a \code{GRangesList} the type of data required
      is driven by the \code{VariantType} class in the \code{region} argument.
      Below is a description of the appropriate \code{GRangesList}
      \code{subject} for each type of \code{region}. 
      \describe{
        \item{CodingVariants :}{coding regions by transcript}
        \item{IntronVariants :}{intron regions by transcript}
        \item{IntergenicVariants :}{transcript regions by gene}
        \item{SpliceSiteVariants :}{intron regions by transcript}
        \item{FiveUTRVariants :}{five prime UTR regions by transcript}
        \item{ThreeUTRVariants :}{three prime UTR regions by transcript}
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
  A \link[IRanges]{DataFrame} with a row for each variant-transcript match.
  Columns are \code{LOCATION}, \code{QUERYID}, \code{TXID}, \code{GENEID} and 
  \code{CDSID}. 

 \describe{
    \item{\code{LOCATION}}{
      Possible locations are \sQuote{coding}, \sQuote{intron}, 
      \sQuote{threeUTR}, \sQuote{fiveUTR}, \sQuote{intergenic}, and
      \sQuote{spliceSite}.

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
      \sQuote{intergenic} variants.
    }
    \item{\code{FOLLOWID}}{
      The id of the gene the query follows. Only present for 
      \sQuote{intergenic} variants.
    }
    All ID values will be \sQuote{NA} for variants with a location of
    \code{transcript_region} or \code{NA}.
  } 
}

\author{Valerie Obenchain <vobencha@fhcrc.org>}

\seealso{
  \code{\link{readVcf}},
  \code{\link{predictCoding}}
}

\examples{
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene 

  ## --------------------------------
  ## TranscriptDb object as subject 
  ## --------------------------------
  ## Read variants from a VCF file 
  fl <- system.file("extdata", "structural.vcf", package="VariantAnnotation")
  vcf <- readVcf(fl, "hg19")

  ## Variant seqlevels do not match the TxDb seqlevels
  rd <- rowData(vcf)
  head(seqlevels(rd))
  head(seqlevels(txdb))
  intersect(seqlevels(rd), seqlevels(txdb))

  ## Use the renameSeqlevels helper to rename seqlevels
  newnames <- paste("chr", seqlevels(vcf), sep="")
  names(newnames) <- seqlevels(vcf)
  vcf_adj <- renameSeqlevels(vcf, newnames) 
  rd_adj <- renameSeqlevels(rd, newnames)

  ## Confirm
  intersect(seqlevels(vcf_adj), seqlevels(txdb))
  intersect(seqlevels(rd_adj), seqlevels(txdb))

  ## VCF object as the query
  loc_coding <- locateVariants(vcf_adj, txdb, CodingVariants())

  ## GRanges as the query
  loc_intron <- locateVariants(rd_adj, txdb, IntronVariants())
  loc_splice <- locateVariants(rd_adj, txdb, SpliceSiteVariants())

  ## --------------------------------
  ## GRangesList object as subject 
  ## --------------------------------
  ## Note the results do not include GENEID. This information
  ## is not available when a GRangesList is used in place of
  ## a TranscriptDb.
  cdsbytx <- cdsBy(txdb)
  locateVariants(rd_adj, cdsbytx, CodingVariants())

  intbytx <- intronsByTranscript(txdb)
  locateVariants(rd_adj, intbytx, IntronVariants())

  ## --------------------------------
  ## Using the cache 
  ## --------------------------------
  ## myenv <- new.env()
  ## files <- list(vcf1, vcf2, vcf3)
  ## lapply(files, 
  ##        function(fl) 
  ##        {
  ##            vcf <- readVcf(fl, "hg19")
  ##            ## modify seqlevels to match TxDb
  ##            oldnms <- seqlevels(vcf)
  ##            newnms <- paste("chr", oldnms, sep="")
  ##            names(newnms) <- oldnms
  ##            vcf_mod <- renameSeqlevels(vcf, newnms)
  ##            locateVariants(vcf_mod, txdb, AllVariants(), cache=myenv)
  ##        }
  ## ) 
}

\keyword{methods}
