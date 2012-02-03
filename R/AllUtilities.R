### =========================================================================
### Helper functions used by more than one method 
### Not exported 
### =========================================================================


.listCumsum <- function(x) {
  x_unlisted <- unlist(x, use.names=FALSE)
  x_cumsum <- cumsum(x_unlisted)
  x_part <- PartitioningByWidth(elementLengths(x))
  x_cumsum - rep(x_cumsum[start(x_part)] - x_unlisted[start(x_part)],
                 width(x_part))
}

.listCumsumShifted <- function(x) {
  cs <- .listCumsum(x)
  shifted <- c(0L, head(cs, -1))
  shifted[start(PartitioningByWidth(elementLengths(x)))] <- 0L
  shifted
}

.getKcol <- function(conn)
{
    sql <- "SELECT value FROM metadata WHERE name='Key column'"
    as.character(dbGetQuery(conn, sql))
}

.sqlIn <- function(vals)
{
    if (length(vals) == 0L)
        return("")
    sql <-
      lapply(seq_len(length(vals)), 
        function(i) {
            v <- vals[[i]]
            if (!is.numeric(v))
            v <- paste("'", v, "'", sep="")
            v
        })
    paste(unlist(sql), collapse = ",")
}

.setActiveSubjectSeq <-
    function(query, subject)
    ## set active state of circular sequences of subject to FALSE;
    ## warn if query contains circular sequences
{
    queryseq <- seqlevels(query)
    circular <- isCircular(TxDb.Hsapiens.UCSC.hg19.knownGene)
    circNames <- intersect(queryseq, names(circular)[circular])
    if (0L != length(circNames))
        warning("circular sequence(s) in query '",
                paste(circNames, sep="' '"), "' ignored")

    isActiveSeq(subject)[] <- FALSE
    isActiveSeq(subject)[setdiff(queryseq, circNames)] <- TRUE
}
