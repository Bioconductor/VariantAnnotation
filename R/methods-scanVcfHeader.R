### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Quote-aware VCF header line parser
###
### scanBcfHeader() (htslib) splits structured header fields on '=' without
### respecting double-quoted values, so a Description or CommandLineOptions
### containing '=' gets silently truncated.  We re-parse the raw header text
### here and patch each DataFrame back to the correct values.
###

## Parse a single structured header body  "<ID=foo,Description="a=b",..."
## into a named character vector, correctly handling "=" inside double quotes.
.parseVcfHeaderBody <- function(body) {
    ## Strip surrounding < >
    body <- sub("^<", "", sub(">$", "", body))

    ## Walk character by character, splitting on ',' outside quotes,
    ## then splitting each token on the FIRST '=' outside quotes.
    chars  <- strsplit(body, "")[[1]]
    n      <- length(chars)
    in_q   <- FALSE
    tokens <- character(0)
    start  <- 1L

    for (i in seq_len(n)) {
        ch <- chars[i]
        if (ch == '"') {
            in_q <- !in_q
        } else if (ch == ',' && !in_q) {
            tokens <- c(tokens, paste(chars[start:(i - 1L)], collapse = ""))
            start  <- i + 1L
        }
    }
    tokens <- c(tokens, paste(chars[start:n], collapse = ""))

    ## For each token split on the first '=' that is NOT inside quotes
    result <- character(0)
    for (tok in tokens) {
        tc   <- strsplit(tok, "")[[1]]
        tn   <- length(tc)
        tq   <- FALSE
        split_at <- NA_integer_
        for (j in seq_len(tn)) {
            if (tc[j] == '"') {
                tq <- !tq
            } else if (tc[j] == '=' && !tq) {
                split_at <- j
                break
            }
        }
        if (!is.na(split_at)) {
            key <- paste(tc[seq_len(split_at - 1L)], collapse = "")
            val <- paste(tc[(split_at + 1L):tn],     collapse = "")
            ## strip surrounding quotes from value
            val <- sub('^"(.*)"$', "\\1", val)
            result[key] <- val
        }
        ## tokens without '=' (e.g. bare flags) are skipped — consistent with
        ## what scanBcfHeader produces
    }
    result
}

## Parse all ##TYPE=<...> meta-lines from raw header text.
## Returns a named list of data.frames (one per TYPE), each row being one
## record, columns being the union of keys seen.  Row names are the ID field.
.parseRawVcfHeader <- function(raw_lines) {
    ## only structured lines:  ##KEY=<...>
    struct <- grep("^##[A-Za-z0-9_]+=<", raw_lines, value = TRUE)
    if (!length(struct))
        return(list())

    ## extract TYPE and BODY
    type  <- sub("^##([A-Za-z0-9_]+)=<.*$", "\\1", struct)
    body  <- sub("^##[A-Za-z0-9_]+=(<.*>)$", "\\1", struct)

    ## parse each line into a named character vector
    parsed <- lapply(body, .parseVcfHeaderBody)

    ## group by TYPE
    types_unique <- unique(type)
    result <- vector("list", length(types_unique))
    names(result) <- types_unique

    for (tp in types_unique) {
        rows <- parsed[type == tp]
        ## union of all keys
        all_keys <- unique(unlist(lapply(rows, names)))
        mat <- matrix(NA_character_, nrow = length(rows), ncol = length(all_keys),
                      dimnames = list(NULL, all_keys))
        for (i in seq_along(rows)) {
            k <- names(rows[[i]])
            mat[i, k] <- rows[[i]][k]
        }
        df <- as.data.frame(mat, stringsAsFactors = FALSE)
        rownames(df) <- if ("ID" %in% colnames(df)) df$ID else seq_len(nrow(df))
        result[[tp]] <- df
    }
    result
}

## Patch the DataFrames produced by scanBcfHeader with correctly-parsed values
## from the raw header lines.
.patchVcfHeader <- function(hdr_list, raw_lines) {
    parsed <- .parseRawVcfHeader(raw_lines)
    if (!length(parsed))
        return(hdr_list)

    for (tp in names(parsed)) {
        ref_df  <- parsed[[tp]]           # correctly-parsed data.frame
        curr_df <- hdr_list[[tp]]         # possibly-truncated DataFrame (or NULL)

        if (is.null(curr_df) || !is(curr_df, "DataFrame"))
            next

        ## For each column that exists in both, overwrite with re-parsed values
        ## matched by row name (ID).
        common_ids <- intersect(rownames(curr_df), rownames(ref_df))
        if (!length(common_ids))
            next

        for (col in intersect(colnames(curr_df), colnames(ref_df))) {
            curr_df[common_ids, col] <- ref_df[common_ids, col]
        }
        hdr_list[[tp]] <- curr_df
    }
    hdr_list
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### scanVcfHeader methods
###

setMethod(scanVcfHeader, "missing",
    function(file, ...)
{
    VCFHeader()
})

setMethod(scanVcfHeader, "character",
    function(file, ...)
{
    if (length(file)) {
        f1       <- .ensure_bgzf(file[[1]])
        hdr      <- scanBcfHeader(f1, ...)[[1]]
        ## Read raw header lines; use the BGZF path for gzip, plain path otherwise
        raw_lines <- if (grepl("\\.gz$", f1, ignore.case = TRUE))
                         readLines(gzcon(file(f1, "rb")))
                     else
                         readLines(file[[1]])
        raw_lines <- raw_lines[startsWith(raw_lines, "##")]
        patched  <- .patchVcfHeader(hdr$Header, raw_lines)
        VCFHeader(hdr$Reference, hdr$Sample, patched)
    } else {
        VCFHeader()
    }
})

setMethod(scanVcfHeader, "TabixFile",
    function(file, ...)
{
    if (isOpen(file)) {
        ## already open: fall back to path-based method which reads raw lines
        scanVcfHeader(path(file), ...)
    } else {
        hdr       <- scanBcfHeader(path(file), ...)[[1]]
        raw_lines <- headerTabix(file)$header
        raw_lines <- raw_lines[startsWith(raw_lines, "##")]
        patched   <- .patchVcfHeader(hdr$Header, raw_lines)
        VCFHeader(hdr$Reference, hdr$Sample, patched)
    }
})
