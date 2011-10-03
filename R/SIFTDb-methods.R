### =========================================================================
### SIFTDb methods 
### =========================================================================

setMethod("keys", "SIFTDb",
    function(x)
    {
        Kcol <- .getKcol(x$conn)
        sql <- paste("SELECT ", Kcol, " FROM siftdata ", sep="") 
        unqsnp <- unique(dbGetQuery(x$conn, sql))[,1] 
        paste("rs", unqsnp, sep="") 
    }
) 

setMethod("cols", "SIFTDb",
    function(x)
    {
        c("rsID", "proteinID", "aaChange", "method", "aa", "prediction",
          "score", "median", "positionSeqs", "totalSeqs") 
    }
) 

setMethod("select", c("SIFTDb", "character", "character"),
    function(db, keys, cols, ...)
    { 
        fmtkeys <- .sqlIn(.formatRSID(keys))
        sql <- paste("SELECT * FROM siftdata WHERE rsid
          IN (", fmtkeys, ")", sep="")
        raw <- dbGetQuery(db$conn, sql)
        .formatSelect(raw, keys=keys, cols=cols)
    }
)

setMethod("select", c("SIFTDb", "missing", "character"),
    function(db, keys = NULL, cols, ...)
    {
        sql <- paste("SELECT * FROM siftdata", sep="")
        raw <- dbGetQuery(db$conn, sql)
        .formatSelect(raw, cols=cols)
    }
)

setMethod("select", c("SIFTDb", "character", "missing"),
    function(db, keys, cols = NULL, ...)
    {
        fmtkeys <- .sqlIn(.formatRSID(keys)) 
        sql <- paste("SELECT * FROM siftdata WHERE rsid IN (",
            fmtkeys, ")", sep="")
        raw <- dbGetQuery(db$conn, sql)
        .formatSelect(raw, keys=keys)
    }
)

setMethod("select", c("SIFTDb", "missing", "missing"),
    function(db, keys = NULL, cols = NULL, ...)
    {
        raw <- dbGetQuery(db$conn, "SELECT * FROM siftdata")
        .formatSelect(raw)
    }
)

.formatRSID <- function(rsid) 
{
    gsub("rs", "", rsid, fixed=TRUE)
}

.formatSelect <- function(raw, keys = NULL, cols = NULL)
{
    repfactor <- nrow(raw)*4
    id <- lapply(as.list(raw$rsid), function(x) rep(x, 4))
    rsid <- paste("rs", unlist(id, use.names=FALSE), sep="") 
    protein_id <- lapply(as.list(raw$protein_id), function(x) rep(x, 4))
    aa_change <- lapply(as.list(raw$aa_change), function(x) rep(x, 4))
    method <- rep(c("BEST HITS", "BEST HITS", "ALL HITS", "ALL HITS"),
        length(raw$rsid)) 

    ## Create single aa column, snp followed by ref
    aa <- lapply(raw$aa_change, function(x) {
        ref_aa <- substr(as.character(x), 1, 1)
        snp_aa <- substr(as.character(x), nchar(as.character(x)), 
            nchar(as.character(x)))
        c(snp_aa, ref_aa, snp_aa, ref_aa)})

    ## Replace missing values in 'prediction' column with 'not scored'

    dat <- raw[,-c(1:3)]
    #grp <- as.factor(sort(rep(1:4, 5)))
    grp <- list(1:5, 6:10, 11:15, 16:20)
    lst <- lapply(as.list(1:nrow(dat)), function(i, dat) {
            do.call(rbind, split(dat[i,], sort(rep(1:4, 5))))}, 
            dat=as.matrix(dat))
    res <- data.frame(do.call(rbind, lst), row.names=NULL)
    colnames(res) <- c("prediction", "score", "median", "positionSeqs",
        "totalSeqs")
 
    df <- data.frame(rsID=unlist(rsid), proteinID=unlist(protein_id), 
            aaChange=unlist(aa_change), method=method, aa=unlist(aa),
            res, row.names=NULL)
    if (!is.null(cols))
        df <- df[,colnames(df) %in% cols] 
    if (!is.null(keys)) {
        missing <- !keys %in% df$rsID
        if (any(missing))
            warning(paste("keys not found in database : ", keys[missing], 
                sep="")) 
    #    reorder <- match(df$rsID, keys[!missing])
    #    df <- df[reorder, ]
    }
    df
}

