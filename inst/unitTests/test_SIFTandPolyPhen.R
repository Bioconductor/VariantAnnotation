library(SIFT.Hsapiens.dbSNP132)
library(PolyPhen.Hsapiens.dbSNP131)
s <- SIFT.Hsapiens.dbSNP132
p <- PolyPhen.Hsapiens.dbSNP131


test_SIFT <- function()
{
    scol <- cols(s)
    checkIdentical(length(scol), 10L) 

    res <- select(s, "rs2142947")
    checkIdentical(nrow(res), 4L)

    res <- select(s, "rs2142947", cols="AACHANGE")
    checkIdentical(nrow(res), 1L)

    res <- 
      suppressWarnings(select(s, keys=c("rs17970171", "INVALID", "rs17970171")))
    checkIdentical(nrow(res), 9L)
    checkTrue(all(res$RSID %in% c("rs17970171", "INVALID"))) 
}

test_PolyPhen <- function()
{
    pcol <- cols(p)
    checkIdentical(length(pcol), 58L) 

    res <- select(p, "rs3026284")
    checkIdentical(nrow(res), 2L)

    res <- select(p, "rs3026284", cols="POS")
    checkIdentical(nrow(res), 1L)

    res <- 
      suppressWarnings(select(p, keys=c("rs3026284", "INVALID", "rs3026284")))
    checkIdentical(nrow(res), 5L)
    checkTrue(all(res$RSID %in% c("rs3026284", "INVALID"))) 
}

