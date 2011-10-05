### =========================================================================
### PolyPhenDb methods 
### =========================================================================

setMethod("keys", "PolyPhenDb",
    function(x)
    {
        sql <- paste("SELECT rsid FROM ppdata ", sep="") 
        unique(dbGetQuery(x$conn, sql))[,1] 
    }
) 

setMethod("cols", "PolyPhenDb",
    function(x)
    {
        dbListFields(conn=x$conn, "ppdata")
    }
) 

setMethod("select", c("PolyPhenDb", "character", "character"),
    function(db, keys, cols, ...)
    {
        fmtkeys <- .sqlIn(keys)
        if (!"rsid" %in% cols)
            cols <- c("rsid", cols) 
        fmtcols <- paste(cols, collapse=",") 
        sql <- paste("SELECT ", fmtcols, " FROM ppdata WHERE rsid
          IN (", fmtkeys, ")", sep="")
        raw <- dbGetQuery(db$conn, sql)
        .formatPPDbSelect(raw, keys=keys) 
    }
)

setMethod("select", c("PolyPhenDb", "missing", "character"),
    function(db, keys = NULL, cols, ...)
    {
        if (!"rsid" %in% cols)
            cols <- c("rsid", cols) 
        fmtcols <- paste(cols, collapse=",") 
        sql <- paste("SELECT ", fmtcols, " FROM ppdata", sep="")
        dbGetQuery(db$conn, sql)
    }
)

setMethod("select", c("PolyPhenDb", "character", "missing"),
    function(db, keys, cols = NULL, ...)
    {
        fmtkeys <- .sqlIn(keys) 
        sql <- paste("SELECT * FROM ppdata WHERE rsid IN (",
            fmtkeys, ")", sep="")
        raw <- dbGetQuery(db$conn, sql)
        .formatPPDbSelect(raw, keys=keys) 
    }
)

setMethod("select", c("PolyPhenDb", "missing", "missing"),
    function(db, keys = NULL, cols = NULL, ...)
    {
        dbGetQuery(db$conn, "SELECT * FROM ppdata")
    }
)

.formatPPDbSelect <- function(raw, keys)
{
    ## restore order
    missing <- (!keys %in% as.character(raw$rsid))
    if (any(missing))
        warning(paste("keys not found in database : ", keys[missing],
            "\n", sep=""))

    lst <- lapply(keys[!missing], function(x, raw) {
        raw[raw$rsid %in% x, ]
    }, raw=raw)
    do.call(rbind, lst)
}



