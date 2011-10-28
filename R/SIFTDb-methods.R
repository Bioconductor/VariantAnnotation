### =========================================================================
### SIFTDb methods 
### =========================================================================

setMethod("keys", "SIFTDb",
    function(x)
    {
        sql <- paste("SELECT rsid FROM siftdata ", sep="") 
        unqsnp <- unique(dbGetQuery(x$conn, sql))[,1] 
        paste("rs", unqsnp, sep="") 
    }
) 

setMethod("cols", "SIFTDb",
    function(x)
    {
        c("RSID", "PROTEINID", "AACHANGE", "METHOD", "AA", "PREDICTION",
          "SCORE", "MEDIAN", "POSTIONSEQS", "TOTALSEQS") 
    }
) 

setMethod("select", "SIFTDb",
    function(x, keys, cols, ...)
    {
        if (missing(cols))
            cols <- NULL
        if (missing(keys)) {
            keys <- NULL
            sql <- paste("SELECT * FROM siftdata", sep="")
        } else {
            fmtkeys <- .sqlIn(.formatRSID(keys))
            sql <- paste("SELECT * FROM siftdata WHERE RSID IN (",
                         fmtkeys, ")", sep="")
        }
        raw <- dbGetQuery(x$conn, sql)
        .formatSIFTDbSelect(raw, keys=keys, cols=cols)
    }
)

.formatRSID <- function(rsid) 
{
    gsub("rs", "", rsid, fixed=TRUE)
}

.formatSIFTDbSelect <- function(raw, keys = NULL, cols = NULL)
{
    ## restore order
    if (!is.null(keys)) {
        fmtkeys <- gsub("rs", "", keys, fixed=TRUE)
        missing <- (!fmtkeys %in% as.character(raw$RSID))
        if (any(missing))
            warning(paste("keys not found in database : ", keys[missing], 
                          "\n", sep="")) 
        reorder <- match(fmtkeys[!missing], as.character(raw$RSID))
        raw <- raw[reorder, ]
    }

    ## format common columns 
    repfactor <- nrow(raw)*4
    id <- lapply(as.list(raw$RSID), function(x) rep(x, 4))
    rsid <- paste("rs", unlist(id, use.names=FALSE), sep="") 
    protein_id <- lapply(as.list(raw$PROTEINID), function(x) rep(x, 4))
    aa_change <- lapply(as.list(raw$AACHANGE), function(x) rep(x, 4))
    method <- rep(c("BEST HITS", "BEST HITS", "ALL HITS", "ALL HITS"),
        length(raw$RSID)) 

    ## create single aa column, snp followed by ref
    aa <- lapply(raw$AACHANGE, 
          function(x) {
              ref_aa <- substr(as.character(x), 1, 1)
              snp_aa <- substr(as.character(x), nchar(as.character(x)), 
                  nchar(as.character(x)))
              c(snp_aa, ref_aa, snp_aa, ref_aa)
          })

    ## partition and stack rows 
    dat <- raw[,-c(1:3)]
    grp <- list(1:5, 6:10, 11:15, 16:20)
    lst <- lapply(as.list(1:nrow(dat)), 
           function(i, dat) {
               do.call(rbind, split(dat[i,], sort(rep(1:4, 5))))
           }, dat=as.matrix(dat))
    res <- data.frame(do.call(rbind, lst), row.names=NULL)
    colnames(res) <- c("PREDICTION", "SCORE", "MEDIAN", "POSITIONSEQS",
                       "TOTALSEQS")
 
    df <- data.frame(RSID=unlist(rsid), PROTEINID=unlist(protein_id), 
                     AACHANGE=unlist(aa_change), METHOD=method, AA=unlist(aa),
                     res, row.names=NULL)
    if (!is.null(cols)) {
      if (!"RSID" %in% cols) cols <- c("RSID", cols)
          df <- df[,colnames(df) %in% cols] 
    }
    df
}
