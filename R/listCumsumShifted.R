
.listCumsumShifted <- function(x) 
{
    cs <- listCumsum(x)
    shifted <- c(0L, head(cs, -1))
    shifted[start(PartitioningByWidth(elementLengths(x)))] <- 0L
    shifted
}

