
# Class - KatDetect -------------------------------------------------------

#' @title KatDetect-class: KatDetect objects
#' @description The katdetectr package introduces a new S4 object which stores all relevant information regarding kataegis detection.
#'
#'
#' @slot kataegisFoci (\link[GenomicRanges]{GRanges}): Contains all
#' annotated putative kataegis foci.
#' @slot genomicVariants (\link[VariantAnnotation]{VRanges}): Contains all
#' processed and annotated genomic variants.
#' @slot segments (\link[GenomicRanges]{GRanges}): Contains all segments
#' detected after changepoint analysis.
#' @slot info (list): Contains some general information and model parameters used for kataegis detection.
#'
#' @examples
#'
#' syntheticData <- generateSyntheticData()
#' kd <- detectKataegis(syntheticData)
#'
#' getKataegisFoci(kd)
#' getGenomicVariants(kd)
#' getSegments(kd)
#' getInfo(kd)
#'
#' @rdname KatDetect
#' @exportClass KatDetect
#'
#' @author Daan Hazelaar
#' @author Job van Riet
#' @export
setClass(
    Class = 'KatDetect',
    slots = methods::representation(
        kataegisFoci = 'GRanges',
        genomicVariants = 'VRanges',
        segments = 'GRanges',
        info = 'list'
    ),
    prototype = base::list(
        kataegisFoci = GenomicRanges::GRanges(),
        genomicVariants = VariantAnnotation::VRanges(),
        segments = GenomicRanges::GRanges(),
        info = base::list()
    )
)


.initkatdetect <- function(genomicVariantsAnnotatedKat, segments, kataegisFoci, info){

    # create new KatDetect object
    kd <- methods::new(
        'KatDetect',
        genomicVariants = genomicVariantsAnnotatedKat,
        segments = segments,
        kataegisFoci = kataegisFoci,
        info = info
    )

    return(kd)
}
