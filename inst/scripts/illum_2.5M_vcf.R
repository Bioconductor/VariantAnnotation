#####  A script to write VCF files for Illumina 2.5M array data

library(genoset)
library(getopt)
library(Rvcftools)
library(multicore)

##Handle Options
opt.matrix = matrix(c(
  'pinfo.file',    'a', 1, 'character','File with probe annotations, e.g. HumanOmni2.5-4v1_hg19_StrandReport.txt',
  'gcscore.cutoff',   'c', 1, 'numeric',  'Min gcscore value for genotype to be used',
  'help',             'h', 0, 'logical',  'Print Usage',
  'output.prefix',    'o', 1, 'character','Prefix for output files',
  'sample.id.column', 'i', 1, 'character','pData column name with sample ID',
  'geno.file',        't', 1, 'character','RData with genoset containing genotype and gcscore'
  ),ncol=5,byrow=T)
opt = getopt( opt.matrix )

#Set defaults
if ( is.null(opt$output.prefix)) { opt$output.prefix = "gCGP.illumina" }
if ( is.null(opt$sample.id.column)) { opt$sample.id.column = "SAMPLE_ID" }
if ( is.null(opt$gcscore.cutoff)) { opt$gcscore.cutoff = 0.15 }
if ( is.null(opt$pinfo.file) ) { opt$pinfo.file = "/gne/research/data/isilon/dnaseq/external_data/SNP_Array/processed/Current/HumanOmni2.5-4v1_hg19_StrandReport.txt" }
if ( is.null(opt$geno.file) ) { opt$geno.file = "/gne/research/data/isilon/dnaseq/external_data/SNP_Array/processed/Current/gCGP.illumina.SNPArray.geno.RData" }


# Load genotypes
load(opt$geno.file)

# Load probe positions for hg19
pinfo = loadPinfo(opt$pinfo.file)  # Load genoset of genotypes and gcscores
pinfo = pinfo[ featureNames(geno), ]

# Assemble components of vcf common to all samples
pinfo.alleles = genotypesFromPinfo(pinfo)
refbase = getRefBase(locData(geno))
allele.options = codeAlleleOptions(refbase,pinfo.alleles)
alt.col = makeAltColumn(allele.options)
generic.body = cbind(chr(geno),pos(geno),featureNames(geno),refbase,alt.col,".",".",".","GT:GQ")

save(generic.body,allele.options,file="illumina.SNPArray.generic.VCF.body.RData")

# Write individual files
dir.create("VCF",showWarnings=FALSE)
format.list = list(
  list(ID="GT",Number=1,Type="String",Description="Genotype"),
  list(ID="GQ",Number=1,Type="Integer",Description="Genotype Quality 100X gcscore")
  )

lapply(sampleNames(geno), function(sample.name) {
  sample.id = paste("SAM",pData(geno)[sample.name,opt$sample.id.column],sep="")
  file.name = paste("VCF/gCGP.illumina.SNPArray",sample.name,sample.id,"vcf",sep=".")
  con = writeVCFHeader(file.name, format.list=format.list, sample.names=sample.id)
  observed = strsplit( assayDataElement(geno,"geno")[,sample.name], "" )
  coded.alleles = codeAlleleObservations2(allele.options, observed)
  qualities = fixNA(assayDataElement(geno,"gcscore")[,sample.name])
  sample.body = paste(coded.alleles,qualities,sep=":")
  full.body = cbind(generic.body, sample.body)
  write.table(full.body,file=con,sep="\t",row.names=FALSE,col.names=FALSE, quote=FALSE)
  close(con)
  return(invisible())
})
