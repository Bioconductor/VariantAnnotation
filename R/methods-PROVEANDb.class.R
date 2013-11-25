### ========================================================================
### PROVEANDb methods 
### =========================================================================

setMethod("keys", "PROVEANDb",
    function(x, keytype, ...)
    {
        if (missing(keytype))
            keytype <- "DBSNPID"
        sql <- paste("SELECT DISTINCT ", keytype, " FROM proveandata ", sep="") 
        dbGetQuery(x$conn, sql)[,1]
    }
) 

setMethod("columns", "PROVEANDb",
    function(x)
        DBI:::dbListFields(x$conn, "proveandata") 
) 

setMethod("keytypes", "PROVEANDb",
    function(x) "DBSNPID"
) 

setMethod("select", "PROVEANDb",
    function(x, keys, columns, keytype, ...)
    {
        if (missing(keytype)) keytype <- "DBSNPID"
        sql <- .createPROVEANDbQuery(x, keys, columns, keytype)
        if (length(sql)) {
            raw <- dbGetQuery(x$conn, sql)
            .formatPROVEANDbSelect(raw, keys, columns, keytype)
        } else {
            data.frame()
        }
    }
)

.createPROVEANDbQuery <- function(x, keys, columns, keytype)
{
    if (missing(keys) && missing(columns)) {
        sql <- "SELECT * FROM proveandata"
    }
    if (!missing(keys)) {
        if(.missingKeys(x, keys, "PROVEAN"))
            return(character())
        if (!missing(columns)) {
            if (.missingCols(x, columns, "PROVEAN"))
                return(character())
            columns <- union(keytype, columns)
            fmtcols <- paste(columns, collapse=",")
            fmtkeys <- .sqlIn(keys)
            sql <- paste("SELECT ", fmtcols, " FROM proveandata WHERE ",
                         keytype, " IN (", fmtkeys, ")", sep="")
        } else {
            fmtkeys <- .sqlIn(keys)
            sql <- paste("SELECT * FROM proveandata WHERE ", keytype,
                         " IN (", fmtkeys, ")", sep="")
        }
    } else {
        if (.missingCols(x, columns, db="PROVEAN"))
            return(character())
        columns <- union(keytype, columns)
        fmtcols <- paste(columns, collapse=",")
        sql <- paste("SELECT ", fmtcols, " FROM proveandata", sep="")
    }
    sql
}

.formatPROVEANDbSelect <- function(raw, keys, columns, keytype)
{
    ## no data
    if (!nrow(raw))
        return(data.frame())

    ## remove duplicate rows
    if (any(dup <- duplicated(raw)))
        raw <- raw[!dup, ]

    ## reorder columns 
    if (!missing(columns)) {
      if (!keytypes %in% columns) columns <- c(keytypes, columns)
          raw <- raw[,colnames(raw) %in% columns] 
    }

    ## return keys not found 
    index <- unique(raw[[keytype]])
    missing <- rep(FALSE, length(index))
    if (!missing(keys)) 
        missing <- (!keys %in% as.character(index))

    lst <- as.list(rep(NA_character_, length(keys)))
    for (i in which(missing == FALSE))
        lst[[i]] <- raw[raw[[keytype]] %in% keys[i], ]

    df <- do.call(rbind, lst)
    df[[keytype]][is.na(df[[keytype]])] <- keys[missing]
    rownames(df) <- NULL
    df
}



