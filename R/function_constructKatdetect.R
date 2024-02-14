#' @title constructKatdetect
#' @description Constructor function for a KatDetect object.
#'
#' @param genomicVariants (\link[VariantAnnotation]{VRanges})
#' @param segments (\link[GenomicRanges]{GRanges})
#' @param kataegisFoci (\link[GenomicRanges]{GRanges})
#' @param info (list)
#'
#' @examples
#' constructKatdetect()
#'
#' @return (KatDetect): Returns a KatDetect object.
#' @author Daan Hazelaar
#'
#' @export
constructKatdetect <- function(genomicVariants = VariantAnnotation::VRanges(), segments = GenomicRanges::GRanges(), kataegisFoci = GenomicRanges::GRanges(), info = list()) {
    # create new KatDetect object
    kd <- methods::new(
        "KatDetect",
        genomicVariants = genomicVariants,
        segments = segments,
        kataegisFoci = kataegisFoci,
        info = info
    )

    return(kd)
}
