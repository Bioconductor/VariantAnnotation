\name{locateVariants}
\alias{locateVariants}

\alias{locateVariants,Ranges,TxDb,VariantType-method}
\alias{locateVariants,Ranges,GRangesList,VariantType-method}
\alias{locateVariants,GRanges,TxDb,VariantType-method}
\alias{locateVariants,GRanges,GRangesList,VariantType-method}
\alias{locateVariants,VCF,TxDb,VariantType-method}
\alias{locateVariants,VCF,GRangesList,VariantType-method}

\alias{locateVariants,GRanges,TxDb,CodingVariants-method}
\alias{locateVariants,GRanges,GRangesList,CodingVariants-method}
\alias{locateVariants,GRanges,TxDb,IntronVariants-method}
\alias{locateVariants,GRanges,GRangesList,IntronVariants-method}
\alias{locateVariants,GRanges,TxDb,FiveUTRVariants-method}
\alias{locateVariants,GRanges,GRangesList,FiveUTRVariants-method}
\alias{locateVariants,GRanges,TxDb,ThreeUTRVariants-method}
\alias{locateVariants,GRanges,GRangesList,ThreeUTRVariants-method}
\alias{locateVariants,GRanges,TxDb,IntergenicVariants-method}
\alias{locateVariants,GRanges,GRangesList,IntergenicVariants-method}
\alias{locateVariants,GRanges,TxDb,SpliceSiteVariants-method}
\alias{locateVariants,GRanges,GRangesList,SpliceSiteVariants-method}
\alias{locateVariants,GRanges,TxDb,PromoterVariants-method}
\alias{locateVariants,GRanges,GRangesList,PromoterVariants-method}
\alias{locateVariants,GRanges,TxDb,AllVariants-method}
\alias{locateVariants,GRanges,GRangesList,AllVariants-method}

\title{Locate variants}

\description{Variant location with respect to gene function}

\usage{
locateVariants(query, subject, region, ...)
\S4method{locateVariants}{VCF,TxDb,VariantType}(query, subject, region, ..., 
    cache=new.env(parent=emptyenv()), ignore.strand=FALSE, asHits=FALSE)
\S4method{locateVariants}{GRanges,TxDb,VariantType}(query, subject, region, ..., 
    cache=new.env(parent=emptyenv()), ignore.strand=FALSE, asHits=FALSE)
}

\arguments{
  \item{query}{A \link[IRanges]{Ranges}, \link[GenomicRanges]{GRanges} or 
    \linkS4class{VCF} object containing the variants. Metadata columns are 
    allowed but ignored.

    NOTE: Zero-width ranges are treated as width-1 ranges; start values
    are decremented to equal the end value.
  }
  \item{subject}{A \link[GenomicFeatures]{TxDb} or \code{GRangesList}
    object that serves as the annotation. GFF files can be converted to
    \link[GenomicFeatures]{TxDb} objects with 
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
    returned as a \link[S4Vectors]{Hits} object. Not applicable when
    \code{region} is AllVariants or IntergenicVariants.
  }
}

\details{
  \describe{
    \item{Range representation :}{
      The ranges in \code{query} should reflect the position(s) of the 
      reference allele. For snps the range will be of width 1. For range 
      insertions or deletions the reference allele could be a sequence 
      such as GGTG in which case the width of the range should be 4. 
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

      \sQuote{intergenic} variants are ranges that do not fall within a 
      defined gene region. \sQuote{transcripts by gene} are extracted from 
      the annotation and overlapped with the variant positions. Variants with 
      no overlaps are classified as \code{intergenic}. When available, gene 
      IDs for the flanking genes are provided as \code{PRECEDEID} and 
      \code{FOLLOWID}. \code{upstream} and \code{downstream} arguments define 
      the acceptable distance from the query for the flanking genes. 
      \code{PRECEDEID} and \code{FOLLOWID} results are lists and contain all 
      genes that fall within the defined distance. See examples for how
      to compute the distance from ranges to PRECEDEID and FOLLOWID.

      \sQuote{promoter} variants fall within a specified range upstream and 
      downstream of the transcription start site. Ranges values can be set 
      with the \code{upstream} and \code{downstream} arguments when creating
      the \code{PromoterVariants()} or \code{AllVariants()} classes.
    }
    \item{Subject as GRangesList :}{
      The \code{subject} can be a \code{TxDb} or \code{GRangesList}
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
  The strand of the output matches the strand of the \code{subject} hit 
  except in the case of IntergenicVariants. For intergenic, multiple precede
  and follow gene ids are returned for each variant. When 
  \code{ignore.strand=TRUE} the return strand is \code{*} because
  genes on both strands are considered and it is possible to have a mixture.
  When \code{ignore.strand=FALSE} the strand will match the \code{query}
  because only genes on the same strand are considered.

  Metadata columns are \code{LOCATION}, \code{QUERYID}, 
  \code{TXID}, \code{GENEID}, \code{PRECEDEID}, \code{FOLLOWID} and 
  \code{CDSID}. Results are ordered by \code{QUERYID}, \code{TXID} and
  \code{GENEID}. Columns are described in detail below. 

 \describe{
    \item{\code{LOCATION}}{
      Possible locations are \sQuote{coding}, \sQuote{intron}, 
      \sQuote{threeUTR}, \sQuote{fiveUTR}, \sQuote{intergenic},
      \sQuote{spliceSite} and \sQuote{promoter}.

      To be classified as \sQuote{coding}, \sQuote{intron}, \sQuote{threeUTR} 
      or \sQuote{fiveUTR} the variant must fall completely within the region.

      \sQuote{intergenic} variants do not fall within a transcript. The 
      \sQuote{GENEID} for these positions are \code{NA}. Lists of flanking
      genes that fall within the distance defined by \code{upstream} and 
      \code{downstream} are given as \sQuote{PRECEDEID} and \sQuote{FOLLOWID}. 

      A \sQuote{spliceSite} variant overlaps any portion of the first 2 or last 
      2 nucleotides of an intron.
    }
    \item{\code{LOCSTART, LOCEND}}{
      Genomic position in LOCATION-centric coordinates. 
      If LOCATION is `intron`, these are intron-centric coordinates,
      if LOCATION is `coding` then cds-centric. All coordinates are
      relative to the start of the transcript. SpliceSiteVariants,
      IntergenicVariants and PromoterVariants have no formal
      extraction `by transcript` so for these variants LOCSTART and
      LOCEND are NA. Coordinates are computed with \code{mapToTranscripts}; 
      see ?\code{mapToTranscripts} in the GenomicFeatures package for details.
    }
    \item{\code{QUERYID}}{
      The \code{QUERYID} column provides a map back to the row in the
      original \code{query}. If the \code{query} was a \code{VCF} object this
      index corresponds to the row in the \code{GRanges} in the \code{rowData}
      slot. 
    }
    \item{\code{TXID}}{
      The transcript id taken from the \code{TxDb} object.
    }
    \item{\code{CDSID}}{
      The coding sequence id(s) taken from the \code{TxDb} object.
    }
    \item{\code{GENEID}}{
      The gene id taken from the \code{TxDb} object.
    }
    \item{\code{PRECEDEID}}{
      IDs for all genes the query precedes within the defined
      \code{upstream} and \code{downstream} distance. Only applicable
      for \sQuote{intergenic} variants.
    }
    \item{\code{FOLLOWID}}{
      IDs for all genes the query follows within the defined
      \code{upstream} and \code{downstream} distance. Only applicable
      for \sQuote{intergenic} variants.
    }
    All ID values will be \sQuote{NA} for variants with a location of
    \code{transcript_region} or \code{NA}.
  } 
}

\author{Valerie Obenchain <vobencha@fhcrc.org>}

\seealso{
  \itemize{
    \item The \link{readVcf} function. 
    \item The \link{predictCoding} function. 
    \item The promoters function on the 
          \link[GenomicRanges]{intra-range-methods} man page in the 
          GenomicRanges package.
  }
}

\examples{
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene 

  ## ---------------------------------------------------------------------
  ## Variants in all gene regions 
  ## ---------------------------------------------------------------------
  ## Read variants from a VCF file.
  fl <- system.file("extdata", "gl_chr1.vcf", 
                    package="VariantAnnotation")
  vcf <- readVcf(fl, "hg19")

  ## Often the seqlevels in the VCF file do not match those in the TxDb.
  head(seqlevels(vcf))
  head(seqlevels(txdb))
  intersect(seqlevels(vcf), seqlevels(txdb))

  ## Rename seqlevels with renameSeqlevesl().
  vcf <- renameSeqlevels(vcf, paste0("chr", seqlevels(vcf))) 

  ## Confirm.
  intersect(seqlevels(vcf), seqlevels(txdb))

  ## Overlaps for all possible variant locations.
  loc_all <- locateVariants(vcf, txdb, AllVariants())
  table(loc_all$LOCATION)

  ## ---------------------------------------------------------------------
  ## Variants in intergenic regions 
  ## ---------------------------------------------------------------------
  ## Intergenic variants do not overlap a gene range in the 
  ## annotation and therefore 'GENEID' is always NA. Flanking genes 
  ## that fall within the 'upstream' and 'downstream' distances are 
  ## reported as PRECEDEID and FOLLOWID.
  region <- IntergenicVariants(upstream=70000, downstream=70000)
  loc_int <- locateVariants(vcf, txdb, region)
  mcols(loc_int)[c("LOCATION", "PRECEDEID", "FOLLOWID")] 

  ## Distance to the flanking genes can be computed for variants that 
  ## have PRECEDEID(s) or FOLLOWID(s). Each variant can have multiple 
  ## flanking id's so we first expand PRECEDEID and the corresponding
  ## variant ranges.
  p_ids <- unlist(loc_int$PRECEDEID, use.names=FALSE) 
  exp_ranges <- rep(loc_int,  elementLengths(loc_int$PRECEDEID))

  ## Provide the variant GRanges as 'x' and the TxDb annotation as 'y' 
  ## to distance(). The 'id' and 'type' arguments describe an id found 
  ## in the TxDb. See the ?'nearest-methods' man page in GenomicFeatures
  ## for details.
  p_dist <- distance(exp_ranges, txdb, id=p_ids, type="gene")
  head(p_dist)

  ## Expanded view of ranges, gene id and distance:
  exp_ranges$PRECEDE_DIST <- p_dist
  exp_ranges

  ## Collapsed view of ranges, gene id and distance:
  loc_int$PRECEDE_DIST <- relist(p_dist, loc_int$PRECEDEID)
  loc_int 

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
  ## to store the extracted components of the TxDb
  ## (i.e., cds by tx, introns by tx etc.). This avoids having to 
  ## re-extract these GRangesLists during each loop.
  \dontrun{
  myenv <- new.env()
  files <- list(vcf1, vcf2, vcf3)
  lapply(files, 
      function(fl) {
          vcf <- readVcf(fl, "hg19")
          ## modify seqlevels to match TxDb
          seqlevels(vcf_mod) <- paste0("chr", seqlevels(vcf)) 
          locateVariants(vcf_mod, txdb, AllVariants(), cache=myenv)
      }) 
  }

  ## ---------------------------------------------------------------------
  ## Parallel implmentation 
  ## ---------------------------------------------------------------------
  \dontrun{
  library(BiocParallel)

  ## A connection to a TxDb object is established when
  ## the package is loaded. Because each process reading from an 
  ## sqlite db must have a unique connection the TxDb
  ## object cannot be passed as an argument when running in 
  ## parallel. Instead the package must be loaded on each worker.

  ## The overhead of the multiple loading may defeat the
  ## purpose of running the job in parallel. An alternative is
  ## to instead pass the appropriate GRangesList as an argument.
  ## The details section on this man page under the heading
  ## 'Subject as GRangesList' explains what GRangesList is
  ## appropriate for each variant type.

  ## A. Passing a GRangesList:

  fun <- function(x, subject, ...)
      locateVariants(x, subject, IntronVariants())

  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  grl <- intronsByTranscript(TxDb.Hsapiens.UCSC.hg19.knownGene)
  mclapply(c(vcf, vcf), fun, subject=grl)


  ## B. Passing a TxDb:
 
  ## Forking:
  ## In the case of forking, the TxDb cannot be loaded
  ## in the current workspace.
  ## To detach the NAMESPACE:
  ##     unloadNamespace("TxDb.Hsapiens.UCSC.hg19.knownGene")

  fun <- function(x) {
      library(TxDb.Hsapiens.UCSC.hg19.knownGene)
      locateVariants(x, TxDb.Hsapiens.UCSC.hg19.knownGene, 
                     IntronVariants())
  }
  mclapply(c(vcf, vcf), fun)

  ## Clusters:
  cl <- makeCluster(2, type = "SOCK")
  fun <- function(query, subject, region) {
      library(VariantAnnotation)
      library(TxDb.Hsapiens.UCSC.hg19.knownGene)
      locateVariants(query, TxDb.Hsapiens.UCSC.hg19.knownGene, region)
  }
  parLapply(cl, c(vcf, vcf), fun, region=IntronVariants())
  stopCluster(cl)
  }
}
\keyword{methods}
