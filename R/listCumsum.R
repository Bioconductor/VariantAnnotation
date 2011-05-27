listCumsum <- function(x) 
{
    x_unlisted <- unlist(x, use.names=FALSE)
    x_cumsum <- cumsum(x_unlisted)
    x_part <- PartitioningByWidth(elementLengths(x))
    x_cumsum - rep(x_cumsum[start(x_part)] - 
        x_unlisted[start(x_part)], width(x_part))
}

