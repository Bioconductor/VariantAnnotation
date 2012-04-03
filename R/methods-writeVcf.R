### =========================================================================
### writeVcf methods 
### =========================================================================

setMethod(writeVcf, c("VCF", "character"),
    function(obj, filename, ...)
{
    hdr <- .makeVcfHeader(obj)
    mat <- .makeVcfMatrix(obj)
 
    con = file(filename, open="w")
    writeLines(hdr, con)
    writeLines(mat, con)
 
    close(con)
})

.makeVcfMatrix <- function(obj)
{
    rd <- rowData(obj)

    CHROM <- as.vector(seqnames(rd))
    POS <- start(rd)
    ID <- .makeVcfID(names(rd))
    REF <- as.character(values(ref(obj))[["REF"]])
    ALT <- values(alt(obj))[["ALT"]]
    if (is(ALT, "DNAStringSetList")) {
        lt <- rep(seq_len(length(ALT)), width(partitioning(ALT)))
        wd <- do.call(c, (lapply(ALT, width)))
        sp <- split(as.character(unlist(ALT, use.names=FALSE)), lt)
        sp[lt[wd == 0]] <- "." 
        ALT <- unlist(lapply(sp, function(i) paste(i, collapse=",")))
    } else if (is(ALT, "CharacterList")) {
        ALT <- unlist(ALT)
    }
    QUAL <- values(qual(obj))[["QUAL"]]
    QUAL[is.na(QUAL)] <- "."
    FILTER <- values(filt(obj))[["FILTER"]]
    INFO <- .makeVcfInfo(values(info(obj))[-1])
    GENO <- .makeVcfGeno(geno(obj))
    dat <- c(CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, GENO)
    lst <- split(dat, seq_len(nrow(obj)))
    .pasteCollapse(CharacterList(lst), collapse="\t")
}

.makeVcfID <- function(id, ...)
{
    idx <- grep(":", id, fixed=TRUE)
    id[idx] <- "."
    id
}

.makeVcfFormat <- function(geno, cls, idx, ...)
{
    if (length(idx) > 0) 
        geno[[idx]] <- matrix(unlist(geno[[idx]], use.names=FALSE), 
            nrow=nrow(geno[[idx]])) 
    map <- Map(function(elt, nms) 
               {
                   apply(elt, 1, 
                     function(i, nms) 
                     {
                       if (!all(is.na(i)))
                           nms
                       else
                           NA
                     }, nms)
               }, as.list(geno), names(geno))

    lst <- split(unlist(map, use.names=FALSE), 
                 rep(seq_len(nrow(geno[[1]])), length(cls)))
    fmt <- lapply(lst, na.omit)
    .pasteCollapse(CharacterList(fmt), collapse=":")
}

.makeVcfGeno <- function(geno, ...)
{
    cls <- lapply(geno, class) 
    idx <- which(cls == "array") 
    FORMAT <- .makeVcfFormat(geno, cls, idx)

    nsub <- ncol(geno[[1]])
    nrec <- nrow(geno[[1]])
    arylst <- lapply(geno[idx], 
                  function(elt, nsub) 
                  {
                      m <- matrix(c(elt), ncol=nsub)
                      matrix(split(m, seq_len(nrec*nsub)), ncol=nsub)
                  }, nsub)
    geno[idx] <- SimpleList(arylst)

    fmtidx <- lapply(strsplit(FORMAT, ":"), function(x) names(geno) %in% x)
    subj <- lapply(seq_len(nsub), 
                function(i, geno) 
                {
                    dat <- unlist(lapply(geno, function(fld) fld[,i]), 
                                  use.names=FALSE, recursive=FALSE)
                    dat <- lapply(dat, function(elt)
                                       {
                                           if (na <- any(is.na(elt)))
                                               elt[na] <- "."
                                           else
                                               elt
                                           paste(elt, collapse=",")
                                        })
                    mat <- matrix(dat, ncol=length(geno))
                    lst <- split(mat, rep(seq_len(nrec), length(geno)))
                    map <- Map("[", lst, fmtidx)
                    .pasteCollapse(CharacterList(map), collapse=":")
                }, geno)
    cbind(FORMAT, do.call(cbind, subj))
}

.makeVcfInfo <- function(info, ...)
{
    cls <- lapply(info, class)
    ## CompressedLists
    cmp <- grep("Compressed", cls, fixed=TRUE)
    for (i in cmp)
        info[[i]] <- lapply(info[[i]], 
            function(elt) 
            {
                if (all(is.na(elt))) 
                    NA
                else
                    paste(as.list(elt), collapse=",")
            }) 
    ## arrays
    arr <- which(cls == "array") 
    for (i in arr) {
        mat <- matrix(info[[i]], nrow=nrow(info[[i]]))
        info[[i]] <- apply(mat, 1, 
            function(elt) 
            {
                if (all(is.na(elt)))
                    NA
                else
                    paste(elt, collapse=",")
            })
    }
    map <- Map(function(elt, nms, cls) {
        if (is(elt, "logical")) {
            lapply(elt, function(x) if(x) nms)
        } else {
            lapply(elt, function(i) 
                        {
                            if (!is.na(i))
                                paste(nms, "=", i, sep="")
                        })
        }
    }, as.list(info), names(info), cls)

    lapply(seq_len(length(map[[1]])), 
        function(i, map) 
        {
            paste(unlist(lapply(map, "[", i), use.names=FALSE), collapse=";")
        }, map)
}

.pasteCollapse <- rtracklayer:::pasteCollapse
 
.makeVcfHeader <- function(obj, ...)
## FIXME : not writing out reference or sample
{
    hdr <- exptData(obj)[["header"]]
    header <- Map(.formatHeader, as.list(header(hdr)),
                  as.list(names(header(hdr))))
    samples <- samples(hdr) 
    colnms <- paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",
                    "INFO", "FORMAT", samples[!is.null(samples)]), collapse="\t")

    unlist(c(header, colnms), use.names=FALSE) 
}

.formatHeader <- function(df, nms)
{
    if (nms == "META") {
        fd <- format(Sys.time(), "%Y%m%d")
        if ("fileDate" %in% rownames(df))
            df[rownames(df) == "fileDate", ] <- fd
        else
            df <- rbind(df, DataFrame(Value=fd, row.names="fileDate")) 
        paste("##", rownames(df), "=", df[,1], sep="")
    } else {
        if ("Description" %in% colnames(df)) {
            df$Description <- paste("\"", df$Description, "\"", sep="")
            prs <- paste(rep(colnames(df), each=nrow(df)), "=",
                         unlist(df, use.names=FALSE), sep="")
            lst <- split(prs, seq_len(nrow(df)))
            lns <- .pasteCollapse(CharacterList(lst), collapse=",") 
            paste("##", nms, "=<ID=", rownames(df), ",", lns, ">", sep="")
        }
    }
}

### ------------------------------------------------------------
#codeAlleleObservations <- function(ref, alt, observed, ...) {
#    if (length(ref) != length(alt))
#        stop("ref and alt must be the same length")
#
#    len <- lapply(ref, length) 
#    if (any(len > 1))
#        stop("each element of ref must be a single nucleotide") 
#
#    options <- .codeAlleleOptions(ref, alt)
#    if (length(options) != length(observed))
#        stop("ref, alt and observed must all be the same length")
#    .Call(.code_allele_observations, options, observed)
#}
#
#.codeAlleleOptions <- function(ref, alt) {
#    mapply(function(x, y) {
#      unique(c(x, y))
#    },ref, alt, USE.NAMES=FALSE)
#}


############
## Functions
############

#codeAlleleOptions <- function(refbase, alleles) {
#  mapply(function(x, y) {
#      unique(c(x, y))
#    },refbase, alleles, USE.NAMES=FALSE)
#}
#
#codeAlleleObservations <- function( allele.options, observed) {
#  if (length(allele.options) != length(observed)) {
#    stop("allele.options and observed must be the same length\n")
#  }
#  coded = .Call(.code_allele_observations, allele.options, observed)
#  return(coded)
#}

##### Functions to create parts of a VCF file
##' @import BSgenome.Hsapiens.UCSC.hg19 IRanges
##' @importFrom Biostrings getSeq
##' @useDynLib Rvcftoolss

##' Read probeset info file
##' Read file from Illumina with probe locs and possible nucleotides
##' @param pinfo.file character, name of file from Illumina with probe location and alleles
##' @param mapping character, genome version (e.g. hg19)
##' @return data.frame of pinfo.file contents with a little tidying up.
##' @author Peter M. Haverty \email{phaverty@@gene.com}
##' @export
loadPinfo <- function(pinfo.file,mapping="hg19") {
  if (mapping == "hg19") {
    pinfo = read.delim(
      pinfo.file,
      colClasses=c("character","character","integer","character","character","character","character"),  # hg19 version
      as.is=TRUE, row.names=1)
    pinfo[ pinfo[,"IlmnDesStrand"] %in% "T","IlmnDesStrand"] = "+"
    pinfo[ pinfo[,"IlmnDesStrand"] %in% "B","IlmnDesStrand"] = "-"
    pinfo = pinfo[ ! pinfo$MapInfo == 0 | ! pinfo$MapInfo == 0, ]
  }  else if (mapping == "hg18") {
    pinfo = read.delim(
      pinfo.file,
      colClasses=c("character","numeric","character","integer","character","character","character","character","character"),
      as.is=TRUE, row.names=1)
    pinfo[ pinfo[,"IlmnDesignStrand"] %in% "T","IlmnDesignStrand"] = "+"
    pinfo[ pinfo[,"IlmnDesignStrand"] %in% "B","IlmnDesignStrand"] = "-"
    pinfo = pinfo[ ! pinfo$MapInfo == 0 | ! pinfo$MapInfo == 0, ]
  } else {
    stop("Mapping ", mapping, "not supported\n")
  }
  return(pinfo)
}

##' Parse probe info file to get nucleotide options
##' Parse probe info file to get nucleotide options on top strand.
##' @param pinfo data.frame, as read by load.pinfo()
##' @return list of characters, e.g. list( c("A","G"), c("C","A"))
##' @author Peter M. Haverty \email{phaverty@@gene.com}
##' @export
genotypesFromPinfo <- function(pinfo) {
  alleles = sub("^.*\\[(.*)\\].*$","\\1",pinfo[,"TopGenomicSeq"])
  alleles = strsplit(alleles, "/")
  # Determine rows where available alleles will need to be flipped
  revcomp = list(A="T",C="G",G="C",T="A")
  flip = (pinfo[,"IlmnDesStrand"] == "+" & pinfo[,"dbSNPdesStrand"] == "R") | (pinfo[,"IlmnDesStrand"] == "-" & pinfo[,"dbSNPdesStrand"] == "F")
  alleles[flip] = lapply( alleles[flip], function(x) { unname(unlist(revcomp[ x ])) } )
  return(alleles)
}

##' Look up reference base.
##' Given a RangedData, look up refrence base at each position.
##' @param gs Genoset, RangedData, or RangesList
##' @return character, nucleotide at each position of rd in reference genome
##' @author Peter M. Haverty \email{phaverty@@gene.com}
##' @export
getRefBase <- function(gs) {
  refbase = getSeq(Hsapiens, start=start(gs), end=end(gs), names=paste("chr", sub("^chr","",as.character(space(gs))),sep=""))
  return(refbase)
}

##' Make ALT column
##' Take lookup of allele options from codeAlleleOptions() and output ALT column
##' string vector.  Basically take out ref option and use "." for only ref option.
##' This contains the sloweset loop ever. Someone please think of a way to do this
##' without mclapply.
##' 
##' @param allele.lookup list, like list( c("A","C"), c("G"), c("T","A","G"))
##' @return character, comma-separated list of alelle options, "." if alt is only ref
##' @author Peter M. Haverty email{phaverty@@gene.com}
##' @export
makeAltColumn <- function(allele.lookup) {
  allele.lengths = elementLengths(allele.lookup)
  alt = rep(".",length(allele.lookup))
  alt[ allele.lengths == 2 ] = unlist( lapply( allele.lookup[ allele.lengths == 2 ], function(x) { return(x[2]) } ) )
  alt[ allele.lengths == 3 ] = unlist( lapply( allele.lookup[ allele.lengths == 3 ], function(x) { sprintf("%s,%s", x[2], x[3]) } ) )
  alt[ allele.lengths == 4 ] = unlist( lapply( allele.lookup[ allele.lengths == 4 ], function(x) { sprintf("%s,%s,%s", x[2], x[3], x[4]) } ) )
  return(alt)
}

##' Write header part of VCF
##' @param filename character, name of file to create and write to
##' @param info.list list of lists. Sub-lists are name/value pairs for INFO rows.
##' @param filter.list list of lists. Sub-lists are name/value pairs for FILTER rows.
##' @param format.list list of lists. Sub-lists are name/value pairs for FORMAT rows.
##' @param sample.names character, names of sample to go in last header row
##' @return nothing
##' @author Peter M. Haverty \email{phaverty@@gene.com}
##' @export
writeVCFHeader <- function(filename, info.list=NULL, filter.list=NULL, format.list=NULL, sample.names) {
  con = file(filename, open="w")
  writeLines(c("##fileformat=VCFv4.0", paste("##fileDate=",format(Sys.time(), "%Y%m%d"),sep=""), "##source=gCGP Illumina Pipeline"),con=con)
  if (!is.null(info.list)) {
    writeLines(headerLines(info.list,"INFO"),con)
  }
  if (!is.null(filter.list)) {
    writeLines(headerLines(filter.list,"INFO"),con)
  }
  if (!is.null(format.list)) {
    writeLines(headerLines(format.list,"FORMAT"),con)
  }
  writeLines( paste("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT", sample.names, sep="\t"), con)
  return(con)
}

##' Format list as VCF header lines
##' @param x list of lists like info.list, format.list, etc. in writeVCFHeader
##' @param type character, like INFO or FORMAT or FILTER
##' @return list of characters
##' @author Peter M. Haverty \email{phaverty@@gene.com}
##' @export
headerLines <- function(x, type) {
  lines = sapply(x, function(y) {
    if ( "Description" %in% names(y) ) {
      y[["Description"]] = gsub("\"","",y[["Description"]])
      y[["Description"]] = paste("\"",y[["Description"]],"\"",sep="")
    }
    beginning = paste("##",type,"=<",sep="")
    end = paste(names(y),y,collapse=",",sep="=")
    new.line = paste(beginning, end, ">",sep="")
    return(new.line)
  })
  return(lines)
}

##' Code NA as "."
##' @param x vector
##' @return vector
##' @export 
##' @author Peter M. Haverty \email{phaverty@@gene.com}
fixNA <- function(x) {
  x = as.character(x)
  x[ is.na(x) ] = "."
  return(x)
}
