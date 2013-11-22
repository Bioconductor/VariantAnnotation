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

## NOTE: The database table has 23 columns. There are two
## 'METHODS' (BEST HITS and ALL HITS) of obtaining related 
## sequences using PSI-BLAST. Each variable outcome was listed
## once for each method which created a very wide non-coehesive
## table. We decided to consolidate the data by these methods 
## to return a more compact result. (Not sure if this was wise 
## or not.) This 'consolidating' of the data could have been 
## done before creating the database or here before returning 
## results to the user. We decided to keep the database table
## table as close to original as possible so curating is done here. 
setMethod("columns", "SIFTDb",
    function(x)
    {
        c("RSID", "PROTEINID", "AACHANGE", "METHOD", "AA", "PREDICTION",
          "SCORE", "MEDIAN", "POSTIONSEQS", "TOTALSEQS") 
    }
) 

setMethod("select", "SIFTDb",
    function(x, keys, columns, ...)
    {
        sql <- .createSIFTDbQuery(x, keys, columns)
        if (length(sql)) {
            raw <- dbGetQuery(x$conn, sql)
            .formatSIFTDbSelect(raw, keys, columns)
        } else {
            data.frame()
        }
    }
)

.createSIFTDbQuery <- function(x, keys, cols)
{
   if (.missingKeys(x, keys, "SIFT"))
       return(character())
   if (.missingCols(x, cols, "SIFT"))
       return(character())
   if (missing(keys)) {
       sql <- paste("SELECT * FROM siftdata", sep="")
   } else {
       fmtkeys <- .sqlIn(.formatRSID(keys))
       sql <- paste("SELECT * FROM siftdata WHERE RSID IN (",
                    fmtkeys, ")", sep="")
   }
}

.formatRSID <- function(rsid) 
{
    gsub("rs", "", rsid, fixed=TRUE)
}

.formatSIFTDbSelect <- function(raw, keys, cols)
{
    missing <- rep(FALSE, length(unique(raw$RSID)))
    if (!missing(keys)) {
        fmtkeys <- .formatRSID(keys) 
        missing <- (!fmtkeys %in% as.character(raw$RSID))
    }

    ## no data
    if (!nrow(raw))
        return(data.frame())

    ## create variable columns 
    dat <- raw[,-c(1:3)]
    res <- matrix(t(dat), ncol=5, byrow=TRUE)
    colnames(res) <- c("PREDICTION", "SCORE", "MEDIAN", "POSITIONSEQS",
        "TOTALSEQS")

    ## create common columns 
    id <- lapply(as.list(raw$RSID), function(x) rep(x, 4L))
    rsid <- paste("rs", unlist(id, use.names=FALSE), sep="") 
    protein_id <- lapply(as.list(raw$PROTEINID), function(x) rep(x, 4L))
    aachange <- lapply(as.list(raw$AACHANGE), function(x) rep(x, 4L))
    aa <- .createAAColumn(raw$AACHANGE) 
    method <- rep(c("BEST HITS", "BEST HITS", "ALL HITS", "ALL HITS"),
        length(raw$RSID)) 
    dd <- data.frame(RSID=unlist(rsid), PROTEINID=unlist(protein_id), 
        AACHANGE=unlist(aachange), METHOD=method, AA=unlist(aa),
        res, row.names=NULL, stringsAsFactors=FALSE)

    ## reorder columns 
    if (!missing(cols)) {
      if (!"RSID" %in% cols) cols <- c("RSID", cols)
          dd <- dd[,colnames(dd) %in% cols] 
    }

    ## return keys not found 
    df <- dd[!duplicated(dd), ]
    lst <- as.list(rep(NA_character_, length(keys)))
    for (i in which(missing == FALSE))
        lst[[i]] <- df[df$RSID %in% keys[i], ]

    df <- do.call(rbind, lst)
    df$RSID[is.na(df$RSID)] <- keys[missing]
    rownames(df) <- NULL
    df
}

.createAAColumn <- function(x)
{
    ## create single aa column, snp followed by ref
    lapply(x, function(elt) {
        ref_aa <- substr(as.character(elt), 1, 1)
        snp_aa <- substr(as.character(elt), nchar(as.character(elt)), 
            nchar(as.character(elt)))
        c(snp_aa, ref_aa, snp_aa, ref_aa)
    })
}
