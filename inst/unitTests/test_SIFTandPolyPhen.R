library(SIFT.Hsapiens.dbSNP132)
library(SIFT.Hsapiens.dbSNP137)
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

test_SIFT_137 <- function()
{
    db <- SIFT.Hsapiens.dbSNP137
    scol <- columns(db)
    checkIdentical(length(scol), 18L) 

    res <- select(db, keys=keys(db)[20:21])
    checkIdentical(nrow(res), 3L)
    checkIdentical(res$DBSNPID, c(665L, 665L, 666L))

    res <- quiet(select(db, keys=c("665", "foo", "666")))
    checkIdentical(res$DBSNPID, c("665", "665", "foo", "666"))

    res <- quiet(select(db, keys=c(665, 666), columns=c("LENGTH", "foo")))
    checkIdentical(res, data.frame())
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

