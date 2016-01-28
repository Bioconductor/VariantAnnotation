library(snpStats) ## SnpMatrix class
quiet <- suppressWarnings

test_gSM_array_GT <- function() {
    mat <- matrix(c(".|.", "0|0", "0|1", "1|0", "1|1",
                    "./.", "0/0", "0/1", "1/0", "1/1"),
                  ncol=2, dimnames=list(1:5,1:2))
    sm <- new("SnpMatrix",
              matrix(as.raw(c(0, 1, 2, 2, 3,
                              0, 1, 2, 2, 3)),
                     nrow=2, byrow=TRUE, dimnames=list(1:2,1:5)))
    ref <- DNAStringSet(rep("A",5))
    alt <- DNAStringSetList("C", "G", "T", "C", "G")
    map <- DataFrame(snp.names=rownames(mat),
                             allele.1=ref, 
                             allele.2=alt,
                             ignore=rep(FALSE,5))
    gtsm <- genotypeToSnpMatrix(mat, ref, alt)
    checkIdentical(sm, gtsm$genotypes)
    checkIdentical(map, gtsm$map)
}

test_gSM_array_GT_2alt <- function() {
    mat <- matrix(c("0|1", "1|0", "1|1",
                    "1/2", "2/1", "2/2"),
                  ncol=2, dimnames=list(1:3,1:2))
    sm <- new("SnpMatrix",
              matrix(as.raw(rep(0,6)),
                     nrow=2, byrow=TRUE, dimnames=list(1:2,1:3)))
    ref <- DNAStringSet(rep("A",3))
    alt <- DNAStringSetList(c("C","G"), c("G","T"), c("T","C"))
    map <- DataFrame(snp.names=rownames(mat),
                             allele.1=ref, 
                             allele.2=alt,
                             ignore=rep(TRUE,3))
    gtsm <- quiet(genotypeToSnpMatrix(mat, ref, alt))
    checkIdentical(sm, gtsm$genotypes)
    checkIdentical(map, gtsm$map)
}

test_gSM_array_GT_nonsnv <- function() {
    mat <- matrix(c("0|0", "0|1", "1|0",
                    "0/0", "0/1", "1/0"),
                  ncol=2, dimnames=list(1:3,1:2))
    sm <- new("SnpMatrix",
              matrix(as.raw(rep(0,6)),
                     nrow=2, byrow=TRUE, dimnames=list(1:2,1:3)))
    ref <- DNAStringSet(c("A","ACG","ACG"))
    alt <- DNAStringSetList("CGT", "G", "GAC")
    map <- DataFrame(snp.names=rownames(mat),
                             allele.1=ref, 
                             allele.2=alt,
                             ignore=rep(TRUE,3))
    gtsm <- quiet(genotypeToSnpMatrix(mat, ref, alt))
    checkIdentical(sm, gtsm$genotypes)
    checkIdentical(map, gtsm$map)
}

test_gSM_VCF_GL <- function() {
    fl <- system.file("extdata", "gl_chr1.vcf", package="VariantAnnotation")
    vcf <- readVcf(fl, "hg19")
    gtsm <- quiet(genotypeToSnpMatrix(vcf, uncertain=TRUE))
    checkIdentical(colnames(vcf), rownames(gtsm$genotypes))
    checkIdentical(rownames(vcf), colnames(gtsm$genotypes))
    checkIdentical(rownames(vcf), gtsm$map$snp.names)
    checkIdentical(ref(vcf), gtsm$map$allele.1)
    checkIdentical(alt(vcf), gtsm$map$allele.2)
    checkEquals(unlist(GLtoGP(geno(vcf)$GL)[1,4]), 
                as.vector(g2post(gtsm$genotypes[4,1])))
}

test_gSM_VCF_PL <- function() {
    fl <- system.file("extdata", "hapmap_exome_chr22.vcf.gz", package="VariantAnnotation")
    vcf <- readVcf(fl, "hg19")
    gtsm <- quiet(genotypeToSnpMatrix(vcf, uncertain=TRUE))
    checkIdentical(colnames(vcf), rownames(gtsm$genotypes))
    checkIdentical(rownames(vcf), colnames(gtsm$genotypes))
    checkIdentical(rownames(vcf), gtsm$map$snp.names)
    checkIdentical(ref(vcf), gtsm$map$allele.1)
    checkIdentical(alt(vcf), gtsm$map$allele.2)
    checkEquals(unlist(PLtoGP(geno(vcf)$PL)[1,4]), 
                as.vector(g2post(gtsm$genotypes[4,1])), tolerance=0.01)
}

test_gSM_VCF_structural <- function() {
    fl <- system.file("extdata", "structural.vcf", package="VariantAnnotation")
    vcf <- readVcf(fl, "hg19")
    gtsn <- quiet(genotypeToSnpMatrix(vcf))
    checkTrue(nrow(gtsn$genotype) == 1L)
    checkTrue(rownames(gtsn$genotype) == "NA00001")
}

test_gSM_VCF_noSamples <- function() {
    fl <- system.file("unitTests", "cases", 
                      "FORMAT_header_no_SAMPLEs.vcf", 
                      package="VariantAnnotation")
    vcf <- readVcf(fl, "hg19")
    gtsm <- quiet(genotypeToSnpMatrix(vcf))
    checkEquals(0, nrow(gtsm$genotypes))
}

test_pSM_valid <- function() {
    probs <- matrix(c(1,0,0,
                      0,1,0,
                      0,0,1,
                      NA,NA,NA),
                    ncol=3, byrow=TRUE,
                    dimnames=list(1:4,c("RR","RA","AA")))
    sm <- new("SnpMatrix", matrix(as.raw(c(1,2,3,0)), nrow=1,
                                  dimnames=list(NULL,1:4)))
    checkIdentical(sm, probabilityToSnpMatrix(probs))
}

test_pSM_invalid <- function() {
    # invalid matrix - probs do not sum to 1
    probs <- matrix(c(1,1,0,
                      0,1,0,
                      0,0,1,
                      NA,NA,NA),
                    ncol=3, byrow=TRUE)
    checkException(probabilityToSnpMatrix(probs))
}

test_pSM_onerow <- function() {
    probs <- matrix(c(1,0,0,
                      NA,NA,NA),
                    ncol=3, byrow=TRUE,
                    dimnames=list(1:2,c("RR","RA","AA")))
    sm <- new("SnpMatrix", matrix(as.raw(c(1,0)), nrow=1,
                                  dimnames=list(NULL,1:2)))
    checkIdentical(sm, probabilityToSnpMatrix(probs))
}

test_GLtoGP_array <- function() {
    probs <- aperm(array(c(0.4,0.3,0.3,
                           0.5,0.1,0.4,
                           0.9,0.05,0.05,
                           0,1,0,
                           0,0,1,
                           1,NA,NA),
                         dim=c(3,3,2)),
                   c(2,3,1))
    gl <- probs
    for (i in 1:nrow(probs)) {
        for (j in 1:ncol(probs)) {
            gl[i,j,] <- log10(probs[i,j,])
        }
    }
    gp <- GLtoGP(gl)
    checkEquals(probs, gp)
}

test_PLtoGP_array <- function() {
    probs <- aperm(array(c(0.4,0.3,0.3,
                           0.5,0.1,0.4,
                           0.9,0.05,0.05,
                           0,1,0,
                           0,0,1,
                           1,NA,NA),
                         dim=c(3,3,2)),
                   c(2,3,1))
    pl <- probs
    for (i in 1:nrow(probs)) {
        for (j in 1:ncol(probs)) {
            pl[i,j,] <- -10*log10(probs[i,j,])
        }
    }
    gp <- PLtoGP(pl)
    checkEquals(probs, gp)
}

test_GLtoGP_matrix <- function() {
    probs <- matrix(c(list(c(0.4,0.3,0.3)),
                      list(c(0.5,0.1,0.4)),
                      list(c(0.9,0.05,0.05)),
                      list(c(0,1,0)),
                      list(c(0,0,1)),
                      list(c(1))),
                    ncol=2)
    gl <- probs
    for (i in 1:nrow(probs)) {
        for (j in 1:ncol(probs)) {
            gl[i,j] <- list(log10(unlist(probs[i,j])))
        }
    }
    gp <- GLtoGP(gl)
    checkEquals(probs, gp)
}

test_PLtoGP_matrix <- function() {
    probs <- matrix(c(list(c(0.4,0.3,0.3)),
                      list(c(0.5,0.1,0.4)),
                      list(c(0.9,0.05,0.05)),
                      list(c(0,1,0)),
                      list(c(0,0,1)),
                      list(c(1))),
                    ncol=2)
    pl <- probs
    for (i in 1:nrow(probs)) {
        for (j in 1:ncol(probs)) {
            pl[i,j] <- list(-10*log10(unlist(probs[i,j])))
        }
    }
    gp <- PLtoGP(pl)
    checkEquals(probs, gp)
}

test_matrixToArray <- function() {
    mat <- matrix(c(list(c(1,2,3)),
                    list(c(4,5,6)),
                    list(c(7,8,9)),
                    list(c(10,11,12)),
                    list(c(13,14)),
                    list(c(15))),
                  ncol=2)
    arr <- VariantAnnotation:::.matrixOfListsToArray(mat)
    for (i in 1:nrow(mat)) {
        for (j in 1:ncol(mat)) {
            n <- elementNROWS(mat[i,j])
            checkEquals(unlist(mat[i,j]), arr[i,j,1:n])
        } 
    }
    TRUE
}
