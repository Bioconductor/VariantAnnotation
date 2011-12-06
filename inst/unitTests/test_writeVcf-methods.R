library(BSgenome.Hsapiens.UCSC.hg19)

test_genotypesFromPinfo <- function() {
  pinfo.ds = data.frame(
    Chr=rep(1,6), MapInfo=seq(start=100,length=6,by=100),
    "SNP.design.strand."=c("[T/C]","[A/G]","[T/C]","[A/G]","[T/C]","[A/G]"),
    IlmnDesStrand=c("+","+","-","+","-","+"), dbSNPdesStrand=c("F","F","F","F","R","R"),
    TopGenomicSeq=c("AAAA[A/C]CCCC","CACA[A/G]CACA","TTTT[A/G]AATA","AAAA[A/G]CCCC","CACA[A/G]CACA","TTTT[A/G]AATA"),
    stringsAsFactors=FALSE)
  alleles = list(c("A","C"),c("A","G"),c("T","C"),c("A","G"),c("A","G"),c("T","C"))
  checkEquals( alleles, genotypesFromPinfo(pinfo.ds) )
}

test_getRefBase <- function() {
  rd = RangedData(ranges=IRanges(start=c(1e6,2e6,15e6),width=1,names=letters[1:3]),space=c("1","2","2"))
  nucleotides = c("T","A","T")
  checkEquals(nucleotides, as.character(suppressWarnings(getRefBase(rd))))
}

test_codeAlleleOptions <- function() {
  refbase = c("A","C","A","G","G")
  alleles = list( c("A","T"), c("A","G"), c("A"), c("C"), c("C","G"))
  checkEquals( list( c("A","T"), c("C","A","G"), c("A"), c("G","C"), c("G","C") ), codeAlleleOptions(refbase,alleles) )
}

test_makeAltColumn <- function() {
  allele.lookup = list( c("A","C"), c("A"), c("T","A","G"))
  checkEquals( c("C",".","A,G"), makeAltColumn(allele.lookup) )
}

test_codeAlleleObservations <- function() {
  allele.lookup = list( c("A","T"), c("C","A","G"), c("A"), c("G","C"), c("G","C"), c("G","A","C","TCAGT") )
  observed = list(c("A","A"),c("A","C"),c("A","A"),c("G","C"),c("-","G"),c("A","TCAGT","C"))
  checkEquals( c("0/0","1/0","0/0","0/1","./0","1/3/2"), codeAlleleObservations(allele.lookup, observed) )
}

test_headerLines <- function() {
  format.list = list(
    list(ID="GT",Number=1,Type="String",Description="Genotype"),
    list(ID="GQ",Number=1,Type="Integer",Description="Genotype Quality")
    )
  lines = c(
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">")
  checkEquals(lines, headerLines(format.list,"FORMAT"))
}
 
test_fixNA <- function() {
  checkEquals( c("1","3","1","4",".","4","."), fixNA(c(1,3,1,4,NA,4,NA)))
}

