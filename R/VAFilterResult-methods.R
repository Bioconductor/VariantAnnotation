
VAFilterResult <-
    function(x=logical(), name=NA_character_,
             input=length(x), passing=sum(x), op=NA_character_)
{
    new("VAFilterResult", x, name=mkScalar(as.character(name)[length(name)]),
        stats=data.frame(Name=as.character(name), Input=input, Passing=passing,
          Op=op, stringsAsFactors=FALSE))
}

setMethod(name, "VAFilterResult", function(x, ...) slot(x, "name"))

setMethod(stats, "VAFilterResult", function(x, ...) slot(x, "stats"))

setMethod("Logic", c("VAFilterResult", "VAFilterResult"),
    function(e1, e2)
{
    x <- callNextMethod()
    s1 <- stats(e1); s2 <- stats(e2)
    op <- as.character(.Generic)
    name <- sprintf("(%s %s %s)", name(e1), op, name(e2))
    s <- rbind(stats(e1), stats(e2),
               data.frame(Name=name, Input=length(x), Passing=sum(x),
                          Op=op, stringsAsFactors=FALSE))
    VAFilterResult(x, s$Name, s$Input, s$Passing, s$Op)
})

setMethod("!", "VAFilterResult",
    function(x)
{
    name <- sprintf("!(%s)", name(x))
    y <- callNextMethod()
    s <- rbind(stats(x),
               data.frame(Name=name, Input=length(y), Passing=sum(y),
                          Op="!", stringsAsFactors=FALSE))
    VAFilterResult(y, s$Name, s$Input, s$Passing, s$Op)
})

setMethod(show, "VAFilterResult",
    function(object) 
{
    cat("class:", class(object), "\n")
    cat("name:", name(object), "\n")
    cat("pass filter:", selectSome(object), "\n")
    cat("stats:\n")
    print(stats(object))
})

