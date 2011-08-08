
.cumsumShifted <- function(x) 
{
    cs <- cumsum(x)
    shifted <- c(0L, head(cs, -1))
    shifted[start(PartitioningByWidth(elementLengths(x)))] <- 0L
    shifted
}

