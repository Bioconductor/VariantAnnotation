##### Functions to create parts of a VCF file

##' @import BSgenome.Hsapiens.UCSC.hg19 IRanges
{}
##' @importFrom Biostrings getSeq
{}
##' @useDynLib Rvcftools

############
## Functions
############

{}  # Roxygen keeps adding junk to the import list if it doesn't hit some code

##' Read probeset info file
##'
##' Read file from Illumina with probe locs and possible nucleotides
##' 
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
##'
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
##'
##' Given a RangedData, look up refrence base at each position.
##' @param gs Genoset, RangedData, or RangesList
##' @return character, nucleotide at each position of rd in reference genome
##' @author Peter M. Haverty \email{phaverty@@gene.com}
##' @export
getRefBase <- function(gs) {
  refbase = getSeq(Hsapiens, start=start(gs), end=end(gs), names=paste("chr", sub("^chr","",as.character(space(gs))),sep=""))
  return(refbase)
}

##' Make lookup for alleles
##'
##' Take references bases (vector of A, C, G, or T) and alleles (list of
##' vectors with allowed letters (i.e list( c("A","C") )). Update list to make reference
##' first position. This lookup will be used to encode each new sample.
##' @param refbase character, one letter per position specifying base in reference genome
##' @param alleles list, one character vector per element in refbase, letters allowed at
##' this position.
##' @return character, like c("0/1","0/0")
##' @author Peter M. Haverty \email{phaverty@@gene.com}
##' @export
codeAlleleOptions <- function(refbase, alleles) {
  all.alleles = mapply(
    function(x,y) {
      return( unique( c(x,y) ) )
    },refbase, alleles, USE.NAMES=FALSE)
}


##' Make ALT column
##'
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

##' Look up integer vals for observed bases
##'
##' .Call version of codeAlleleObservations.
##' Take references bases (vector of A, C, G, or T) and alleles (list of
##' vectors with allowed letters (i.e list( c("A","C") )). Return column
##' for observed alleles in one sample.
##' @param allele.options list from codeAlleleOptions
##' @param observed character, like c("AA","AC","CA") from genotypes for a given sample
##' @return character, like c("0/1","0/0") for reference, alt option 1, then homref, etc.
##' @author Peter M. Haverty \email{phaverty@@gene.com}
##' @export
codeAlleleObservations <- function( allele.options, observed) {
  if (length(allele.options) != length(observed)) {
    stop("allele.options and observed must be the same length\n")
  }
  coded = .Call("code_allele_observations", allele.options, observed)
  return(coded)
}

##' Write header part of VCF
##'
##' Write header part of VCF
##' 
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
##'
##' Format list as VCF header lines
##' 
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
##'
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
