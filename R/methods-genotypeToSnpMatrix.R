### =========================================================================
### genotypeToSnpMatrix methods 
### =========================================================================

## Coding for snpMatrix :
## 0 = missing OR multiallelic OR multi-ALT values
## 1 = homozygous reference (0|0 or 0/0) 
## 2 = heterozygous (0|1 or 0/1 or 1|0 or 1/0) 
## 3 = homozygous alternate (risk) allele (1|1 or 1/1)

## empty matrix to return if conditions not met
.emptySnpMatrix <- function() {
    list(genotype=new("SnpMatrix"), 
         map=DataFrame(snp.names=character(), 
                       allele.1=DNAStringSet(), 
                       allele.2=DNAStringSetList(), 
                       ignore=character()))
}

setMethod("genotypeToSnpMatrix", "CollapsedVCF",
          function(x, uncertain=FALSE, ...)
{
    ok <- suppressWarnings(require("snpStats", quietly=TRUE, 
                                   character.only=TRUE))
    ok || stop("'snpStats' required; try biocLite('snpStats')", call.=FALSE) 

    alt <- alt(x)
    if (is(alt, "CompressedCharacterList")) {
        alt <- .toDNAStringSetList(alt)
        if (all(elementNROWS(alt) == 0L)) {
            warning("No nucleotide ALT values were detected.")
            return(.emptySnpMatrix())
        }
    }
    ref <- ref(x)

    if (ncol(x) == 0) {
        warning("no samples in VCF")
    }

    if (!uncertain) {
        gt <- geno(x)$GT
    } else {
        geno.cols <- row.names(geno(metadata(x)[["header"]]))
        if ("GP" %in% geno.cols) {
            gt <- geno(x)$GP
            if (storage.mode(gt) == "list") {
                gt <- .matrixOfListsToArray(gt)
            }
        } else if ("GL" %in% geno.cols) {
            gt <- geno(x)$GL
            if (storage.mode(gt) == "list") {
                gt <- .matrixOfListsToArray(gt)
            }
            gt <- GLtoGP(gt)
        } else if ("PL" %in% geno.cols) {
            gt <- geno(x)$PL
            if (mode(gt) == "list") {
                gt <- .matrixOfListsToArray(gt)
            }
            gt <- PLtoGP(gt)
        } else {
            warning("uncertain=TRUE requires GP, GL or PL; returning NULL")
            return(.emptySnpMatrix())
        }
    }

    callGeneric(gt, ref, alt)
})

setMethod("genotypeToSnpMatrix", "array",
          function(x, ref, alt, ...)
{
    if (!is(ref, "DNAStringSet"))
        stop("'ref' must be a DNAStringSet")
    if (!is(alt, "DNAStringSetList"))
        stop("'alt' must be a DNAStringSetList")
    # query ref and alt alleles for valid SNPs
    altelt <- elementNROWS(alt) == 1L 
    snv <- .testForSNV(ref, alt) 
 
    # if x is a matrix, we have GT with a single value for each snp
    if (is.matrix(x)) {
        if (!all(altelt)) {
            warning("variants with >1 ALT allele are set to NA")
            x[!altelt,] <- ".|."
        }

        if (!all(snv)) {
            message("non-single nucleotide variations are set to NA")
            x[!snv,] <- ".|."
        }
        map <- .genotypeToIntegerSNV(TRUE)
        diploid <- x %in% names(map)
        if (!all(diploid)) {
            warning("non-diploid variants are set to NA")
            x[!diploid] <- ".|."
        }

        mat <- matrix(map[x], nrow=ncol(x), ncol=nrow(x),
                      byrow=TRUE, dimnames=rev(dimnames(x)))
        genotypes <- new("SnpMatrix", mat)
    } else {
    # if x is a 3D array, we have GP with multiple values for each snp

        if (!all(altelt)) {
            warning("variants with >1 ALT allele are set to NA")
            x[!altelt,,] <- NA
        }

        if (!all(snv)) {
            message("non-single nucleotide variations are set to NA")
            x[!snv,,] <- NA
        }

        # if there is more than one ALT allele for any variant,
        # the 3rd dimension of the array will be too big
        # any values here should already have been set to NA above
        if (dim(x)[3] > 3) {
            x <- x[,,1:3]
        }
 
        # for each sample, call probabilityToSnpMatrix
        smlist <- list()
        for (s in 1:ncol(x)) {
            sm <- probabilityToSnpMatrix(x[,s,])
            rownames(sm) <- colnames(x)[s]
            smlist[[s]] <- sm
        }
        genotypes <- do.call(rbind, smlist)
    }
 
    flt <- !(snv & altelt)
    map <- .createMap(rownames(x), ref, alt, flt)

    list(genotypes = genotypes, map = map)
})

.createMap <- function(nms, ref, alt, flt)
{
    if (is.null(ref))
        DataFrame(snp.names=character(0), 
                  allele.1=DNAStringSet(), 
                  allele.2=DNAStringSetList(),
                  ignore=logical())
    else 
        DataFrame(snp.names=nms, 
                  allele.1=ref, 
                  allele.2=alt,
                  ignore=flt)
}

probabilityToSnpMatrix <- function(probs) {
    ok <- suppressWarnings(require("snpStats", quietly=TRUE, 
                                   character.only=TRUE))
    ok || stop("'snpStats' required; try biocLite('snpStats')", call.=FALSE) 

    if (ncol(probs) != 3)
        stop("input matrix should have 3 columns: P(A/A), P(A/B), P(B/B)")

    # skip missing values when checking for validity of probabilities
    missing <- rowSums(is.na(probs)) > 0
    if (!isTRUE(all.equal(rowSums(probs[!missing,,drop=FALSE]),
                          rep(1,sum(!missing)),
                          check.attributes=FALSE,
                          check.names=FALSE)))
        stop("sum of probabilities in each row of input matrix should = 1")

    # post2g can't handle missing data
    if (sum(missing) > 0) {
        probs[missing,] <- 0
        g <- post2g(probs)
        g[missing] <- as.raw(0)
    } else {
        g <- post2g(probs)
    }
    g <-  matrix(g, nrow=1, dimnames=list(NULL, rownames(probs)))
    new("SnpMatrix", g)
}

GLtoGP <- function(gl) {
    if (is.matrix(gl) && storage.mode(gl) == "list") {
        gp <- gl
        for (i in 1:length(gp)) {
            gp[[i]] <- 10^gl[[i]] / sum(10^gl[[i]], na.rm=TRUE)
        }
        gp
    } else if (is.array(gl) & length(dim(gl)) == 3) {
        aperm(apply(gl, c(1,2), function(x) 
                  10^x / sum(10^x, na.rm=TRUE)), 
              c(2,3,1))
    } else {
        stop("gl must be a matrix of lists or a 3D array")
    }
}

## PL is same as GL except for a factor of -10
## GL = log10(L), PL = -10*log10(L) (phred-scaled)
PLtoGP <- function(pl) {
    if (is.matrix(pl) && storage.mode(pl) == "list") {
        gl <- pl
        for (i in 1:length(gl)) {
            gl[[i]] <- pl[[i]]/(-10)
        }
    } else {
        gl <- pl/(-10)
    }
    GLtoGP(gl)
}

.listMatrixToArray <- function(x) {
  n <- elementNROWS(x)
  maxn <- max(n)
  v <- unlist(x, use.names=FALSE)
  a <- array(as(NA, class(v)), dim=c(maxn, nrow(x), ncol(x)),
             dimnames=list(NULL, rownames(x), colnames(x)))
  a[as.integer(IRanges(seq(1L, length(a), maxn), width=n))] <- v
  aperm(a, c(2,3,1))
}

.matrixOfListsToArray <- function(x) {
    # find number of elements of each cell of x
    n <- elementNROWS(x)
    maxn <- max(n)
 
    # for cells with less than the max number of elements, add NAs
    idx <- n < maxn
    x[idx] <- lapply(x[idx], function(a){c(a, rep(NA, maxn-length(a)))})
 
    # unlist and convert to array
    x <- array(unlist(x), dim=c(maxn, nrow(x), ncol(x)),
               dimnames=list(NULL, rownames(x), colnames(x)))
    x <- aperm(x, c(2,3,1))

    x
}
