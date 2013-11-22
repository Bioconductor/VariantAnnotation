### ========================================================================
### PROVEANDb methods 
### =========================================================================

## FIXME: make index reference generic

setMethod("keys", "PROVEANDb",
    function(x)
    {
        sql <- paste("SELECT DISTINCT DBSNPID FROM siftdata ", sep="") 
        dbGetQuery(x$conn, sql)[,1]
    }
) 

setMethod("columns", "PROVEANDb",
    function(x)
        DBI:::dbListFields(x$conn, "siftdata") 
) 

setMethod("select", "PROVEANDb",
    function(x, keys, columns)
    {
        sql <- .createPROVEANDbQuery(x, keys, columns)
        if (length(sql)) {
            raw <- dbGetQuery(x$conn, sql)
            .formatPROVEANDbSelect(raw, keys, columns)
        } else {
            data.frame()
        }
    }
)

.createPROVEANDbQuery <- function(x, keys, columns)
{
    if (missing(keys) && missing(columns)) {
        sql <- "SELECT * FROM siftdata"
    }
    if (!missing(keys)) {
        if(.missingKeys(x, keys, "PROVEAN"))
            return(character())
        if (!missing(columns)) {
            if (.missingCols(x, columns, "PROVEAN"))
                return(character())
            columns <- union("DBSNPID", columns)
            fmtcols <- paste(columns, collapse=",")
            fmtkeys <- .sqlIn(keys)
            sql <- paste("SELECT ", fmtcols, " FROM siftdata WHERE ",
                         "DBSNPID IN (", fmtkeys, ")", sep="")
        } else {
            fmtkeys <- .sqlIn(keys)
            sql <- paste("SELECT * FROM siftdata WHERE DBSNPID ",
                         "IN (", fmtkeys, ")", sep="")
        }
    } else {
        if (.missingCols(x, columns, "PROVEAN"))
            return(character())
        columns <- union("DBSNPID", columns)
        fmtcols <- paste(columns, collapse=",")
        sql <- paste("SELECT ", fmtcols, " FROM siftdata", sep="")
    }
    sql
}

.formatPROVEANDbSelect <- function(raw, keys, columns)
{
    ## no data
    if (!nrow(raw))
        return(data.frame())

    ## remove duplicate rows
    if (any(dup <- duplicated(raw)))
        raw <- raw[!dup, ]

    ## reorder columns 
    if (!missing(columns)) {
      if (!"DBSNPID" %in% columns) columns <- c("DBSNPID", columns)
          raw <- raw[,colnames(raw) %in% columns] 
    }

    ## return keys not found 
    index <- unique(raw$DBSNPID)
    missing <- rep(FALSE, length(index))
    if (!missing(keys)) 
        missing <- (!keys %in% as.character(index))

    lst <- as.list(rep(NA_character_, length(keys)))
    for (i in which(missing == FALSE))
        lst[[i]] <- raw[raw$DBSNPID %in% keys[i], ]

    df <- do.call(rbind, lst)
    df$DBSNPID[is.na(df$DBSNPID)] <- keys[missing]
    rownames(df) <- NULL
    df
}



