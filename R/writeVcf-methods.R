### =========================================================================
### writeVcf methods 
### =========================================================================

setMethod(writeVcf, c("SummarizedExperiment", "character"),
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
    REF <- as.character(values(rd)[["REF"]])
    ALT <- unlist(values(rd)[["ALT"]])
    QUAL <- values(rd)[["QUAL"]]
    FILTER <- values(rd)[["FILTER"]]
    INFO <- .makeVcfInfo(values(rd)[!colnames(values(rd)) %in% 
                         c("REF", "ALT", "QUAL", "FILTER")])
    GENO <- .makeVcfGeno(assays(obj))

    dat <- c(CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, GENO)
    dat <- gsub("NA", ".", dat)
    lst <- split(dat, seq_len(nrow(obj)))
    .pasteCollapse(CharacterList(lst), collapse="\t")
}

.makeVcfID <- function(id, ...)
{
    idx <- grep(":", id, fixed=TRUE)
    id[idx] <- "."
    id
}

.makeVcfFormat <- function(geno, ...)
{
    nms <- names(geno)
    idx <- lapply(seq_len(length(geno)), 
             function(i) rowSums(is.na(geno[[i]])) + i)
    cds <- do.call(c, idx)
    lst <- split(nms[cds], rep(seq_len(nrow(geno[[1]])), length(nms)))
    fmt <- lapply(lst, na.omit)
    .pasteCollapse(CharacterList(fmt), collapse=":")
}

.makeVcfGeno <- function(geno, ...)
{
    FORMAT <- .makeVcfFormat(geno)

    nsub <- ncol(geno[[1]])
    nrec <- nrow(geno[[1]])
    cls <- lapply(geno, class) 
    ary <- which(cls == "array") 
    arylst <- lapply(geno[ary], function(elt, geno) {
      ilst <- sapply(seq_len(nsub), function(i, elt) {
        d <- c(elt[,i,])
        s <- split(d, rep(seq_len(nrec), nsub)) 
        .pasteCollapse(CharacterList(s), collapse=",")
      }, elt, USE.NAMES=FALSE)
      data.frame(ilst, stringsAsFactors=FALSE)
    }, geno) 
    geno[ary] <- SimpleList(arylst)

    subj <- lapply(seq_len(nsub), function(i) {
      dat <- unlist(lapply(geno, function(fld) fld[,i]), use.names=FALSE)
      mat <- matrix(dat, ncol=length(geno))
      mat <- gsub("NA", NA, mat)
      lst <- split(mat, rep(seq_len(nrec), nsub))
      rmna <- lapply(lst, na.omit)
      .pasteCollapse(CharacterList(rmna), collapse=":") 
    })

    cbind(FORMAT, do.call(cbind, subj)) 
}

.makeVcfInfo <- function(info, ...)
{
    cls <- lapply(info, class) 
    cmp <- grep("Compressed", cls, fixed=TRUE)
    for (i in cmp) 
        info[,i] <- .pasteCollapse(info[,i], collapse=",")

    key <- rep(names(info), each=nrow(info))
    vlu <- unlist(info, use.names=FALSE)
    prs <- paste(key, "=", vlu, sep="")
    ## logical TRUE = present 
    idx <- which(vlu == TRUE)
    prs[idx] <- key[idx] 
    prs <- gsub("FALSE", NA, prs)
    ## missings
    prs <- gsub("NA", NA, prs)

    lst <- split(prs, seq_len(nrow(info)))
    rmna <- lapply(lst, na.omit) ## FIXME : better way to omit missings
    .pasteCollapse(CharacterList(rmna), collapse=";") 
}

.pasteCollapse <- rtracklayer:::pasteCollapse
 
.makeVcfHeader <- function(obj, ...)
{
    hdr <- exptData(obj)[["HEADER"]]
    header <- Map(.formatHeader, as.list(hdr), names(hdr))
    samples <- colnames(obj) 
    colnms <- paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",
                    "INFO", "FORMAT", samples[!is.null(samples)]), collapse="\t")

    unlist(c(header, colnms), use.names=FALSE) 
}

.formatHeader <- function(df, nms)
{
    if (nms == "META") {
        fd <- format(Sys.time(), "%Y%m%d")
        df[rownames(df) == "fileDate", ] <- fd
        paste("##", rownames(df), "=", df[,1], sep="")
    } else {
        if (any(colnames(df) %in% "Description")) 
            df$Description <- paste("\"", df$Description, "\"", sep="")
        prs <- paste(rep(colnames(df), each=nrow(df)), "=", unlist(df,
                     use.names=FALSE), sep="")
        lst <- split(prs, seq_len(nrow(df)))
        lns <- .pasteCollapse(CharacterList(lst), collapse=",") 
        paste("##", nms, "=<ID=", rownames(df), ",", lns, ">", sep="")
    }
}

## ------------------------------------------------------------
codeAlleleObservations <- function(ref, alt, observed, ...) {
    if (length(ref) != length(alt))
        stop("ref and alt must be the same length")

    len <- lapply(ref, length) 
    if (any(len > 1))
        stop("each element of ref must be a single nucleotide") 

    options <- .codeAlleleOptions(ref, alt)
    if (length(options) != length(observed))
        stop("ref, alt and observed must all be the same length")
    .Call(.code_allele_observations, options, observed)
}

.codeAlleleOptions <- function(ref, alt) {
    mapply(function(x, y) {
      unique(c(x, y))
    },ref, alt, USE.NAMES=FALSE)
}


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
