
.vaValidity <- function(object) TRUE

setGeneric(".vaValidity")

setClass(".VAUtil",
         representation=representation("VIRTUAL"))

setClass("VAFilter",
         contains=c("function", ".VAUtil"),
         representation=representation(
           name="ScalarCharacter"),
         validity=.vaValidity)

setClass("VAFilterResult",
         contains=c("logical", ".VAUtil"),
         representation=representation(
           name="ScalarCharacter",
           stats="data.frame"))

.SIFTDb <- setRefClass("SIFTDb",
    contains = "AnnotationDb",
)


