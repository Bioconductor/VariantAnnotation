library(SIFT.Hsapiens.dbSNP132)
library(PolyPhen.Hsapiens.dbSNP131)
quiet <- suppressWarnings

test_SIFT_132 <- function()
{
    db <- SIFT.Hsapiens.dbSNP132
    scol <- columns(db)
    checkIdentical(length(scol), 10L) 

    res <- select(db, "rs2142947")
    checkIdentical(nrow(res), 4L)

    res <- select(db, "rs2142947", columns="AACHANGE")
    checkIdentical(nrow(res), 1L)

    res <- 
      quiet(select(db, keys=c("rs17970171", "INVALID", "rs17970171")))
    checkIdentical(nrow(res), 9L)
    checkTrue(all(res$RSID %in% c("rs17970171", "INVALID"))) 
}

test_PolyPhen <- function()
{
    db <- PolyPhen.Hsapiens.dbSNP131
    pcol <- columns(db)
    checkIdentical(length(pcol), 58L) 

    res <- select(db, "rs3026284")
    checkIdentical(nrow(res), 2L)

    res <- select(db, "rs3026284", columns="POS")
    checkIdentical(nrow(res), 1L)

    res <- 
      suppressWarnings(select(db, keys=c("rs3026284", "INVALID", "rs3026284")))
    checkIdentical(nrow(res), 5L)
    checkTrue(all(res$RSID %in% c("rs3026284", "INVALID"))) 
}
