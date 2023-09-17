# Accessor methods ----------------------------------------------------------

#' @rdname getGenomicVariants
#' @exportMethod getGenomicVariants
setMethod("getGenomicVariants", "KatDetect", function(x) {
    return(x@genomicVariants)
})

#' @rdname getSegments
#' @exportMethod getSegments
setMethod("getSegments", "KatDetect", function(x) {
    return(x@segments)
})

#' @rdname getKataegisFoci
#' @exportMethod getKataegisFoci
setMethod("getKataegisFoci", "KatDetect", function(x) {
    return(x@kataegisFoci)
})

#' @rdname getInfo
#' @exportMethod getInfo
setMethod("getInfo", "KatDetect", function(x) {
    return(x@info)
})


# summary and show functions ---------------

#' @export
summary.KatDetect <- function(object, ...) {
    base::cat("Sample name:                                ", getInfo(object)$sampleName, "\n")
    base::cat("Total number of genomic variants:           ", getInfo(object)$totalGenomicVariants, "\n")
    base::cat("Total number of putative Kataegis foci:     ", getInfo(object)$totalKataegisFoci, "\n")
    base::cat("Total number of variants in a Kataegis foci:", getInfo(object)$totalVariantsInKataegisFoci)
}

#' @export
print.KatDetect <- function(x, ...) {
    return(summary(x, ...))
}

setMethod("show", "KatDetect", function(object) {
    base::cat("Class \"KatDetect\" : KatDetect Object\n")
    base::cat("                  : S4 class containing", base::length(base::attributes(object)) - 1, "slots with names:\n")
    base::cat("                   ", base::names(base::attributes(object))[seq_len(base::length(base::attributes(object))) - 1], "\n\n")
    base::cat("Created on:        ", getInfo(object)$date, "\n")
    base::cat("katdetectr version:", getInfo(object)$version, "\n\n")
    base::cat("summary: \n--------------------------------------------------------\n")
    summary(object)
    base::cat("\n--------------------------------------------------------")
})
