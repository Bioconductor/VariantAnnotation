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
    function(x, keys, cols, keytype, ...)
    {
        sql <- .createPPDbQuery(keys, cols) 
        raw <- dbGetQuery(x$conn, sql)
        .formatPPDbSelect(raw, keys) 
    }
)

.createPPDbQuery <- function(keys, cols)
{
    if (missing(keys) && missing(cols)) {
        sql <- "SELECT * FROM ppdata"
    }
    if (!missing(keys) && !missing(cols)) {
        if (!"RSID" %in% cols)
            cols <- c("RSID", cols)
        fmtcols <- paste(cols, collapse=",")
        fmtkeys <- .sqlIn(keys)
        sql <- paste("SELECT ", fmtcols, " FROM ppdata WHERE RSID 
            IN (", fmtkeys, ")", sep="")
    }
    if (!missing(keys) && missing(cols)) {
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
    sql
}

.formatPPDbSelect <- function(raw, keys)
{
    if (missing(keys)) {
        df <- data.frame(raw)
        rownames(df) <- NULL
        df
    } else {
    ## restore key order
        missing <- (!keys %in% as.character(raw$RSID))
        if (any(missing)) {
            msg <- paste(IRanges:::selectSome(keys[missing], 5), collapse=" ")
            warning(sprintf(msg, "keys not found in database : %s"))
        }
        lst <- as.list(rep(NA_character_, length(keys)))
        raw <- raw[!duplicated(raw), ]
        for (i in which(missing == FALSE))
            lst[[i]] <- raw[raw$RSID %in% keys[i], ]

        df <- do.call(rbind, lst)
        df$RSID[is.na(df$RSID)] <- keys[missing]
        rownames(df) <- NULL
        df
    }
}

duplicateRSID <- function(db, keys, ...)
{
    fmtrsid <- .sqlIn(keys)
    sql <- paste("SELECT * FROM duplicates WHERE RSID IN (",
                 fmtrsid, ")", sep="")
    q1 <- dbGetQuery(db$conn, sql)

    fmtgp <- .sqlIn(unique(q1$DUPLICATEGROUP))
    gpsql <- paste("SELECT * FROM duplicates WHERE DUPLICATEGROUP IN (",
                   fmtgp, ")", sep="")
    q2 <- dbGetQuery(db$conn, gpsql)

    matched <- q2[!q2$RSID %in% keys, ]
    matchedlst <- split(matched$RSID, matched$DUPLICATEGROUP)
    names(matchedlst) <- q1$RSID[match(names(matchedlst), q1$DUPLICATEGROUP)]

    missing <- !keys %in% q2$RSID
    if (any(missing)) {
        warning(paste("keys not found in database : ", keys[missing],
                      sep=""))
        missinglst <- list(rep(NA, sum(missing)))
        names(missinglst) <- keys[missing]
        matchedlst <- c(matchedlst, missinglst)
    }

    matchedlst[order(match(names(matchedlst), keys))]
}

