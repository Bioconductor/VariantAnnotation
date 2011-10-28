### =========================================================================
### PolyPhenDb methods 
### =========================================================================

setMethod("keys", "PolyPhenDb",
    function(x)
    {
        sql <- paste("SELECT RSID FROM ppdata ", sep="") 
        unique(dbGetQuery(x$conn, sql))[,1] 
    }
) 

setMethod("cols", "PolyPhenDb",
    function(x)
    {
        dbListFields(conn=x$conn, "ppdata")
    }
) 

setMethod("select", "PolyPhenDb",
    function(x, keys, cols, ...)
    {
        if (missing(keys) & missing(cols)) {
            sql <- "SELECT * FROM ppdata"
        }
        if (!missing(keys) & !missing(cols)) {
            if (!"RSID" %in% cols)
                cols <- c("RSID", cols) 
            fmtcols <- paste(cols, collapse=",") 
            fmtkeys <- .sqlIn(keys)
            sql <- paste("SELECT ", fmtcols, " FROM ppdata WHERE RSID 
                         IN (", fmtkeys, ")", sep="")
        } 
        if (!missing(keys) & missing(cols)) {
            fmtkeys <- .sqlIn(keys)
            sql <- paste("SELECT * FROM ppdata WHERE RSID IN (",
                         fmtkeys, ")", sep="")
        }
        if (missing(keys) & !missing(cols)) {
            if (!"RSID" %in% cols)
                cols <- c("RSID", cols) 
            fmtcols <- paste(cols, collapse=",") 
            sql <- paste("SELECT ", fmtcols, " FROM ppdata", sep="")
        }
        if (missing(keys)) {
            raw <- dbGetQuery(x$conn, sql)
        } else { 
            raw <- dbGetQuery(x$conn, sql)
            .formatPPDbSelect(raw, keys=keys) 
        }
    }
)

.formatPPDbSelect <- function(raw, keys)
{
    ## restore order
    missing <- (!keys %in% as.character(raw$RSID))
    if (any(missing))
        warning(paste("keys not found in database : ", keys[missing],
                      "\n", sep=""))

    lst <- lapply(keys[!missing], 
             function(x, raw) {
               raw[raw$RSID %in% x, ]
             }, raw=raw)
    do.call(rbind, lst)
}



