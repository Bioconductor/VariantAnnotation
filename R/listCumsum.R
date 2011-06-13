listCumsum <- function(x) 
{
    xUnlisted <- unlist(x, use.names=FALSE)
    xCumsum <- cumsum(xUnlisted)
    xPart <- PartitioningByWidth(elementLengths(x))
    xCumsum - rep(xCumsum[start(xPart)] - 
        xUnlisted[start(xPart)], width(xPart))
}

