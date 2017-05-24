### =========================================================================

### VRanges: Variants in a GRanges
### -------------------------------------------------------------------------
###

### Thoughts on the gVCF format.
## The gVCF format is essentially a run-length encoding of wildtype
## calls. We could represent the run as a range, but we would need
## somewhat to indicate "WT". Probably best to go the VCF route: have
## the first ref base, and an NA alt. We can distinguish this from an
## SNV or indel, because the width of the range is > 1, but there is
## only a single ref base. This will break the assertion that the ref
## length always equal the range width, but that is not a big deal.

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessors
###

setMethod("ref", "VRanges", function(x) x@ref)
setReplaceMethod("ref", "VRanges", function(x, value) {
    x@ref <- S4Vectors:::recycleVector(value, length(x))
    x
})
setMethod("alt", "VRanges", function(x) x@alt)
setReplaceMethod("alt", "VRanges", function(x, value) {
  x@alt <- as(.rleRecycleVector(value, length(x)), "characterOrRle")
  x
})
setMethod("sampleNames", "VRanges", function(object) object@sampleNames)
setReplaceMethod("sampleNames", "VRanges", function(object, value) {
  object@sampleNames <- as(.rleRecycleVector(value, length(object)),
                           "factorOrRle")
  object
})
setMethod("totalDepth", "VRanges", function(x) x@totalDepth)
`totalDepth<-` <- function(x, value) {
  x@totalDepth <- as(.rleRecycleVector(value, length(x)), "integerOrRle")
  x
}
setMethod("altDepth", "VRanges", function(x) x@altDepth)
`altDepth<-` <- function(x, value) {
  x@altDepth <- as(.rleRecycleVector(value, length(x)), "integerOrRle")
  x
}
setMethod("refDepth", "VRanges", function(x) x@refDepth)
`refDepth<-` <- function(x, value) {
  x@refDepth <- as(.rleRecycleVector(value, length(x)), "integerOrRle")
  x
}
setMethod("softFilterMatrix", "VRanges", function(x) x@softFilterMatrix)
setReplaceMethod("softFilterMatrix", "VRanges", function(x, value) {
  x@softFilterMatrix <- value
  x
})
setMethod("hardFilters", "VRanges", function(x) x@hardFilters)
setReplaceMethod("hardFilters", "VRanges", function(x, value) {
  x@hardFilters <- value
  x
})
setMethod("called", "VRanges", function(x) {
  rowSums(softFilterMatrix(x)) == ncol(softFilterMatrix(x))
})

setMethod("altFraction", "VRanges", function(x) {
  altDepth(x) / totalDepth(x)
})

setMethod("refFraction", "VRanges", function(x) {
    refDepth(x) / totalDepth(x)
})

setMethod("isDeletion", "VRanges", function(x, ...) 
    .dispatchSNV_ExpandedVCF(.isDeletion, x)
)

setMethod("isInsertion", "VRanges", function(x, ...) 
    .dispatchSNV_ExpandedVCF(.isInsertion, x)
)

setMethod("isIndel", "VRanges", function(x, ...)
    .dispatchSNV_ExpandedVCF(.isIndel, x)
)

setMethod("isDelins", "VRanges", function(x, ...)
    .dispatchSNV_ExpandedVCF(.isDelins, x)
)

setMethod("isSNV", "VRanges", function(x, ...)
    .dispatchSNV_ExpandedVCF(.isSNV, x)
)

setMethod("isSubstitution", "VRanges", function(x, ...)
    .dispatchSNV_ExpandedVCF(.isSubstitution, x)
)

setMethod("isTransition", "VRanges", function(x, ...)
    .dispatchSNV_ExpandedVCF(.isTransition, x)
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

VRanges <-
  function(seqnames = Rle(), ranges = IRanges(),
           ref = character(), alt = NA_character_,
           totalDepth = NA_integer_, refDepth = NA_integer_,
           altDepth = NA_integer_, ..., sampleNames = NA_character_,
           softFilterMatrix = FilterMatrix(matrix(nrow = length(gr), ncol = 0L),
             FilterRules()),
           hardFilters = FilterRules())
{
  gr <- GRanges(seqnames, ranges,
                strand = .rleRecycleVector("*", length(ranges)), ...)
  if (length(gr) == 0L && length(ref) == 0L) {
    maxLen <- 0L
  } else {
    maxLen <- max(length(gr), length(ref), length(alt),
                  length(totalDepth), length(refDepth), length(altDepth),
                  length(sampleNames), nrow(softFilterMatrix))
  }
  if (length(gr) != maxLen)
    gr <- rep(gr, length.out = maxLen)
  ref <- as.character(ref)
  ref <- S4Vectors:::recycleVector(ref, maxLen)
  alt <- .rleRecycleVector(alt, maxLen)
  alt <- as(alt, "characterOrRle")
  refDepth <- .rleRecycleVector(refDepth, maxLen)
  altDepth <- .rleRecycleVector(altDepth, maxLen)
  totalDepth <- .rleRecycleVector(totalDepth, maxLen)
  softFilterMatrix <- as.matrix(softFilterMatrix)
  storage.mode(softFilterMatrix) <- "logical"
  softFilterMatrix.ind <-
    S4Vectors:::recycleVector(seq_len(nrow(softFilterMatrix)), maxLen)
  softFilterMatrix <- softFilterMatrix[softFilterMatrix.ind,,drop=FALSE]
  sampleNames <- .rleRecycleVector(sampleNames, maxLen)
  totalDepth <- as(totalDepth, "integerOrRle")
  refDepth <- as(refDepth, "integerOrRle")
  altDepth <- as(altDepth, "integerOrRle")
  if (any(naToZero(refDepth) + naToZero(altDepth) > totalDepth, na.rm = TRUE))
    warning("'refDepth' + 'altDepth' exceeds 'totalDepth'; using GATK?")
  new("VRanges", gr, ref = ref, alt = alt,
      totalDepth = totalDepth,
      refDepth = refDepth,
      altDepth = altDepth,
      softFilterMatrix = softFilterMatrix,
      sampleNames = as(sampleNames, "factorOrRle"),
      hardFilters = hardFilters)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion
###

setMethod("relistToClass", "VRanges", function(x) "CompressedVRangesList")

parseFilterStrings <- function(x) {
  x[x == "."] <- NA
  if (all(x == "PASS", na.rm=TRUE))
    return(FilterMatrix(matrix(nrow = length(x), ncol = 0L), FilterRules()))
  filterSplit <- strsplit(x, ";", fixed=TRUE)
  filters <- unlist(filterSplit)
  filterNames <- setdiff(unique(filters[!is.na(filters)]), "PASS")
  filterMat <- matrix(TRUE, length(x), length(filterNames),
                      dimnames = list(NULL, filterNames))
  filterMat[cbind(togroup(PartitioningByWidth(filterSplit)), match(filters, filterNames))] <- FALSE
  filterMat[is.na(x),] <- NA
  filterMat
}

genoToMCol <- function(x) {
  if (length(dim(x)) == 3) {
    I(matrix(x, nrow(x) * max(ncol(x), 1L), dim(x)[3]))
  } else {
    if (ncol(x) == 0L) {
      if (is.list(x)) {
        x <- as.logical(x)
      }
      x <- Rle(x[NA], nrow(x))
    } else {
      dim(x) <- NULL
      if (is.list(x))
        x <- as(x, "List")
      x
    }
  }
}

setAs("VCF", "VRanges", function(from) {
  from <- expand(from)
  rd <- rowRanges(from)
  seqnames <- seqnames(rd)
  ranges <- ranges(rd)
  ref <- rd$REF
  alt <- rd$ALT
  if (is(alt, "DNAStringSetList") || is(alt, "CharacterList"))
    alt <- unlist(alt)
  alt <- as.character(alt)
  alt[!nzchar(alt)] <- NA
  alt[alt == "."] <- NA
  ad <- geno(from)$AD
  refDepth <- NA_integer_
  altDepth <- NA_integer_
  if (is.array(ad) && ncol(from) > 0L) {
    if (length(dim(ad)) == 3L) {
        refDepth <- ad[,,1,drop=FALSE]
        altDepth <- ad[,,2,drop=FALSE]
    } else if (length(dim(ad)) == 2L) {
        RD <- geno(from)$RD
        if (is.matrix(RD)) {
            refDepth <- RD
        }
        altDepth <- ad
    }
  }
  totalDepth <- geno(from)$DP
  if (is.null(totalDepth) || ncol(from) == 0L)
    totalDepth <- NA_integer_
  nsamp <- ncol(from)
  if (ncol(from) > 0L) {
    if (!is.null(colnames(from))) {
      sampleNames <- Rle(colnames(from), rep(nrow(from), nsamp))
    } else {
      sampleNames <- as.character(seq_len(ncol(from)))
    }
  } else sampleNames <- NA_character_
  meta <- info(from)
  if (!is.null(mcols(rd)$QUAL)) {
    meta <- cbind(mcols(rd)["QUAL"], meta)
  }
  if (!is.null(meta$END)) {
    end(ranges)[!is.na(meta$END)] <- meta$END[!is.na(meta$END)]
    meta$END <- NULL
  }
  rownames(meta) <- NULL
  if (!is.null(rd$FILTER))
    filter <- parseFilterStrings(rd$FILTER)
  else filter <- FilterMatrix(matrix(nrow = nrow(from), ncol = 0L),
                              FilterRules())
  if (nsamp > 1L) {
    meta <- meta[rep(seq_len(nrow(meta)), nsamp),,drop=FALSE]
    filter <- filter[rep(seq_len(nrow(filter)), nsamp),,drop=FALSE]
    seqnames <- rep(seqnames, nsamp)
    ranges <- rep(ranges, nsamp)
    alt <- rep(alt, nsamp)
    ref <- rep(ref, nsamp)
  } else if (nsamp == 0L && !is.null(meta$DP)) {
    totalDepth <- meta$DP
    meta$DP <- NULL
  }
  if (!is.null(geno(from)[["FT"]]))
    filter <- cbind(filter, parseFilterStrings(as.vector(geno(from)[["FT"]])))
  otherGeno <- geno(from)[setdiff(names(geno(from)), c("AD", "DP", "FT"))]
  if (length(otherGeno) > 0L)
    meta <- DataFrame(meta, lapply(otherGeno, genoToMCol))
  vr <- VRanges(seqnames, ranges, ref, alt,
                totalDepth, refDepth, altDepth,
                hardFilters = FilterRules(), sampleNames = sampleNames,
                softFilterMatrix = filter, meta)
  seqinfo(vr) <- seqinfo(from)
  vr
})

setAs("VRanges", "VCF", function(from) {
  asVCF(from)
})

optionalDescriptions <- c(
  GT = "Genotype",
  GQ = "Genotype quality",
  PL = "Normalized, Phred-scaled likelihoods for genotypes",
  MIN_DP = "Minimum DP observed within the gVCF block",
  END = "End position of the gVCF WT run"
  )
  
makeFORMATheader <- function(x) {
  fixed <-
    DataFrame(row.names = c("AD", "DP", "FT"),
              Number = c(2L, 1L, 1L),
              Type = c("Integer", "Integer", "String"),
              Description =
              c("Allelic depths (number of reads in each observed allele)",
                "Total read depth", "Variant filters"))
  df <- rbind(fixed, makeINFOheader(x))
  df$Type[df$Type == "Flag"] <- "Integer" # Flags not allowed in FORMAT
  df
}

makeINFOheader <- function(x) {
  df <- mcols(x)
  if (is.null(df))
    df <- new("DataFrame", nrows = length(x))
  rownames(df) <- names(x)
  numberForColumn <- function(xi) {
    if (length(dim(xi)) == 2)
      ncol(xi)
    else if (is.list(xi) || is(xi, "List"))
      "."
    else 1L
  }
  typeForColumn <- function(xi) {
    if (is.integer(xi))
      "Integer"
    else if (is.numeric(xi))
      "Float"
    else if (is.logical(xi))
      "Flag"
    else if (is.character(xi) || is.factor(xi)) {
      if (all(nchar(as.character(xi)) == 1L, na.rm=TRUE))
        "Character"
      else "String"
    }
  }
  df$Number <- as.character(lapply(x, numberForColumn))
  df$Type <- as.character(lapply(x, typeForColumn))
  if (is.null(df$Description))
    df$Description <- optionalDescriptions[rownames(df)]
  df[c("Number", "Type", "Description")]
}

makeFILTERheader <- function(x) {
  mat <- softFilterMatrix(x)
  if (!is(mat, "FilterMatrix"))
    return(DataFrame())
  rules <- filterRules(mat)
  df <- mcols(rules)
  if (is.null(df)) {
    df <- new("DataFrame", nrows = length(rules))
  }
  if (is.null(df$Description)) {
    exprs <- as.logical(lapply(rules, is.expression))
    df <- df[exprs,]
    df$Description <- as.character(lapply(rules[exprs], function(expr) {
      gsub("\"", "\\\"", gsub("\\", "\\\\", deparse(expr[[1]]), fixed=TRUE),
           fixed=TRUE)
    }))
  }
  rownames(df) <- names(rules[exprs])
  df
}

## If a row contains a FALSE, all FALSE filters are listed;
## otherwise, if a row contains an NA or it is empty, the result is
## NA (.), otherwise PASS.
makeFILTERstrings <- function(x) {
  failed <- which(!x)
  ftNames <- as.character(colnames(x)[col(x)][failed])
  ftStrings <- rep(NA_character_, nrow(x))
  passedRows <- rowSums(x) > 0L
  failedRows <- rowSums(!x, na.rm=TRUE) > 0L
  ftStrings[passedRows & !failedRows] <- "PASS"
  ftStrings[failedRows] <-
    unstrsplit(splitAsList(ftNames, row(x)[failed]), ";")
  ftStrings
}

isGVCFRun <- function(x) {
  nchar(ref(x)) == 1L & width(x) > 1L
}

gVCFRunEnds <- function(x) {
  ifelse(isGVCFRun(x), end(x), NA_integer_)
}

vranges2Vcf <- function(x, info = character(), filter = character(),
                        meta = character())
{
  if (!is.character(info) || any(is.na(info)) ||
      any(!info %in% colnames(mcols(x))))
    stop("'info' must be a character vector naming the mcols to include",
         " in the INFO column of the VCF file; the rest become FORMAT fields.")
  if (!is.character(filter) || any(is.na(filter)) ||
      any(!filter %in% colnames(softFilterMatrix(x))))
    stop("'filter' must be a character vector naming the filters to include",
         " in the FILTER column of the VCF file; the rest are encoded as the FT",
         " FORMAT field.")
  if (!is.character(meta) || any(is.na(meta)) ||
      any(!meta %in% names(metadata(x))))
    stop("'filter' must be a character vector naming metadata(x) elements",
         " to include as metadata in the VCF header.")
  metaStrings <- as.character(sapply(metadata(x)[meta], as.character))
  if (any(elementNROWS(metaStrings) != 1L))
    stop("The elements named in 'meta' must be of length 1")

  if (any(is.na(sampleNames(x)))) {
    stop("sampleNames(x) must not contain missing values (for VCF export)")
  }
  sampleLevels <- levels(sampleNames(x))
  if (length(sampleLevels) > 1) {
    xUniq <- unique(x)
  } else {
    xUniq <- x
  }

  END <- gVCFRunEnds(xUniq)
  if (!all(is.na(END))) {
    if (!is.null(xUniq$END)) {
      warning("Replacing 'END' metadata/info column with computed gVCF END")
    }
    xUniq$END <- END
    info <- union(info, "END")
  }
  
  rowRanges <- xUniq
  mcols(rowRanges) <- NULL
  rowRanges <- as(rowRanges, "GRanges", strict=TRUE)

  colData <- DataFrame(Samples = seq_along(sampleLevels),
                       row.names = sampleLevels)

  meta_vector <- c(fileformat = "VCFv4.1",
                   source = paste("VariantAnnotation",
                     packageVersion("VariantAnnotation")),
                   phasing = "unphased", 
                   metaStrings)
  meta_header <- DataFrame(Value = meta_vector, row.names = names(meta_vector))
  genoMCols <- setdiff(names(mcols(x)), info)
  header <- VCFHeader(reference = seqlevels(x), samples = sampleLevels,
                      header = DataFrameList(META = meta_header,
                        FORMAT = makeFORMATheader(mcols(x)[genoMCols]),
                        INFO = makeINFOheader(mcols(xUniq)[info]),
                        FILTER = makeFILTERheader(x)))
  metadata <- list(header = header)

  alt <- as.character(alt(xUniq))
  qual <- xUniq$QUAL
  if (is.null(qual) || !is.numeric(qual))
    qual <- rep.int(NA_real_, length(xUniq))
  filtMat <- softFilterMatrix(xUniq)
  if (ncol(filtMat) > 0L)
    filtMat <- filtMat[,filter,drop=FALSE]
  filterStrings <- makeFILTERstrings(filtMat)
  fixed <- DataFrame(REF = DNAStringSet(ref(xUniq)),
                     ALT = alt,
                     QUAL = qual,
                     FILTER = filterStrings)

  if (length(sampleLevels) > 1) {
    samplesToUniq <- tapply(x, sampleNames(x), function(xi) {
      match(xi, xUniq)
    }, simplify=FALSE)
  }
  
  genoArray <- function(v) {
    v <- as.vector(v) # handles e.g. Rle
    if (is.logical(v)) {
      v <- as.integer(v)
    }
    width <- length(v)/length(x)
    if (width > 1L) {
      v <- matrix(v, nrow=length(x), ncol=width)
    }
    if (length(sampleLevels) > 1L) {
      sample.ord <- order(sampleNames(x))
      if (is.matrix(v)) {
        v <- v[sample.ord,,drop=FALSE]
      } else {
        v <- v[sample.ord]
      }
      ind <- unlist(samplesToUniq, use.names=FALSE) +
        (togroup(PartitioningByWidth(samplesToUniq))-1L)*length(xUniq)
      a <- matrix(NA_integer_, nrow=length(xUniq)*length(sampleLevels),
                  ncol=width)
      a[ind,] <- v
    } else {
      a <- v
    }
    if (width > 1L) {
      array(a, c(length(xUniq), length(sampleLevels), ncol(a)))
    } else {
      matrix(a, length(xUniq), length(sampleLevels))
    }
  }

  alleleDepth <- c(refDepth(x), altDepth(x))
  filtMat <- softFilterMatrix(x)
  if (ncol(filtMat) > 0L)
    filtMat <- filtMat[,setdiff(colnames(filtMat), filter), drop=FALSE]
  ftStrings <- makeFILTERstrings(filtMat)
  geno <-
    SimpleList(AD = genoArray(alleleDepth),
               DP = genoArray(totalDepth(x)),
               FT = genoArray(ftStrings))
  if (length(genoMCols) > 0L)
    geno <- c(geno, SimpleList(lapply(mcols(x)[genoMCols], genoArray)))

  info <- mcols(xUniq)[info]
 
  VCF(rowRanges = rowRanges, colData = colData, metadata = metadata,
      fixed = fixed, geno = geno, info = info, collapsed = FALSE)
}

setMethod("asVCF", "VRanges", vranges2Vcf)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Reading from VCF
###

VRangesScanVcfParam <- function(fixed="ALT", info=NA, geno="AD", ...) {
  .Defunct("ScanVcfParam")
}

readVcfAsVRanges <- function(x, genome, param=ScanVcfParam(),
                             use.names=FALSE, ...)
{
  if (!isTRUEorFALSE(use.names)) {
    stop("'use.names' must be TRUE or FALSE")
  }
  as(readVcf(x, genome, param=param, row.names=use.names, ...), "VRanges")
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Writing to VCF (see methods-writeVcf.R)
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### gVCF support
###

addTotalDepthRuns <- function(x, runs, genome) {
  mins <- viewMins(runs)
  gr <- as(ranges(runs), "GRanges")
  ref <- getSeq(genome, resize(gr, 1L))
  vr <- VRanges(seqnames(gr), ranges, ref, alt, totalDepth=mins)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Filtering
###

softFilter <- function(x, filters, ...) {
  softFilterMatrix(x) <-
    cbind(softFilterMatrix(x),
          FilterMatrix(matrix = evalSeparately(filters, x, ...),
                       filterRules = filters))
  x
}

resetFilter <- function(x) {
  softFilterMatrix(x) <- FilterMatrix(matrix(nrow = length(x), ncol = 0L),
                                      filterRules = FilterRules())
  x
}

setMethod("subsetByFilter", c("VRanges", "FilterRules"), function(x, filter) {
  ans <- callNextMethod(x, filter)
  hardFilters(ans) <- c(hardFilters(ans), filter)
  ans
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Utilities
###

setMethod("duplicated", "VRanges",
          function (x, incomparables = FALSE, fromLast = FALSE,
                    method = c("auto", "quick", "hash"))
          {
            if (!identical(incomparables, FALSE))
              stop("\"duplicated\" method for VRanges objects ",
                   "only accepts 'incomparables=FALSE'")
            duplicatedIntegerQuads(as.factor(seqnames(x)),
                                   as.factor(alt(x)),
                                   start(x), width(x),
                                   fromLast = fromLast,
                                   method = method)
          })

setMethod("match", c("VRanges", "VRanges"),
          function(x, table, nomatch = NA_integer_, incomparables = NULL,
                   method = c("auto", "quick", "hash"))
          {
            if (!isSingleNumberOrNA(nomatch))
              stop("'nomatch' must be a single number or NA")
            if (!is.integer(nomatch))
              nomatch <- as.integer(nomatch)
            if (!is.null(incomparables))
              stop("\"match\" method for VRanges objects ",
                   "only accepts 'incomparables=NULL'")
            merge(seqinfo(x), seqinfo(table))
            altLevels <- as.character(union(alt(x), alt(table)))
            x_seqnames <- GenomicRanges:::relevelSeqnamesForMatch(x, table)
            matchIntegerQuads(x_seqnames, factor(as.character(alt(x)),
                                                 altLevels),
                              start(x), width(x),
                              as.factor(seqnames(table)),
                              factor(as.character(alt(table)),
                                     altLevels),
                              start(table), width(table),
                              nomatch = nomatch, method = method)
          })

setMethod("%in%", c("VRanges", "TabixFile"),
          function(x, table)
          {
            table <- readVcfAsVRanges(table, genome=genome(x), param=x)
            x %in% table
          })

setMethod("tabulate", "VRanges", function(bin, nbins) {
  if (!missing(nbins))
    stop("'nbins' argument not relevant")
  m <- match(bin, bin)
  ## for dupes, 'm' always points to the subject with the lowest index
  tab <- tabulate(m, length(bin))
  ans <- bin[tab > 0]
  ans$sample.count <- tab[tab > 0]
  ans
})

setMethod("merge", c("VRanges", "VRanges"), function(x, y, ...) {
### FIXME: support softFilterMatrix
  xdf <- as.data.frame(x)
  ydf <- as.data.frame(y)
  bothAllNA <- function(nm) all(is.na(slot(x, nm))) && all(is.na(slot(y, nm)))
  ignore.cols <- Filter(bothAllNA, c("refDepth", "altDepth", "totalDepth"))
  ignore.cols <- c(ignore.cols, c("width", "strand"))
  xdf <- xdf[setdiff(colnames(xdf), ignore.cols)]
  ydf <- ydf[setdiff(colnames(ydf), ignore.cols)]
  by <- c("seqnames", "start", "end", "ref", "alt", "sampleNames")
  merged <- merge(xdf, ydf, by = by, ...)
  with(merged, VRanges(seqnames, IRanges(start, end), ref, alt, NA, NA, NA,
                       sampleNames = sampleNames,
                       merged[setdiff(colnames(merged), by)]))
})

setMethod("liftOver", c("VRanges", "Chain"), function(x, chain, ...) {
  grl <- liftOver(GRanges(seqnames(x), ranges(x)), chain, ...)
  gr <- unlist(grl, use.names=FALSE)
  ans <- x[togroup(PartitioningByWidth(grl))]
  ans <- update(ans, seqinfo = seqinfo(gr), seqnames = seqnames(gr),
                ranges = ranges(gr))
  relist(ans, grl)
})

pileupGRanges <- function(x) {
  x <- x[!isIndel(x)]
  gr <- GRanges(seqnames(x), ranges(x), strand(x))
  sm <- selfmatch(gr)
  uniq <- sm == seq_len(length(gr))
  map <- integer()
  map[sm[uniq]] <- seq_len(sum(uniq))
  pos <- map[sm]
  base <- match(c(alt(x), ref(x)), DNA_BASES)
  depth <- c(altDepth(x), refDepth(x))
  samp <- as.factor(sampleNames(x))
  pileup <- array(0L, c(sum(uniq), length(levels(samp)), length(DNA_BASES)),
                  dimnames=list(NULL, levels(samp), base=DNA_BASES))
  pileup[cbind(pos, samp, base)[!is.na(base),]] <- depth[!is.na(base)]
  ugr <- gr[uniq]
  ugr$ref[pos] <- ref(x)
  ugr$pileup <- pileup
  ugr
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Show
###

.showHardFilters <- function(object) {
  cat(S4Vectors:::labeledLine("  hardFilters", names(hardFilters(object))))
}

setMethod("show", "VRanges", function(object) {
  callNextMethod()
  .showHardFilters(object)
})

setMethod(GenomicRanges:::extraColumnSlotNames, "VRanges",
          function(x) {
            c("ref", "alt", "totalDepth", "refDepth", "altDepth",
              "sampleNames", "softFilterMatrix")
          })
