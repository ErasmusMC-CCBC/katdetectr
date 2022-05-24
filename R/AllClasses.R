
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
#' new('KatDetect', genomicVariants = VariantAnnotation::VRanges(), info = list())
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
#' @importFrom rlang .data
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


.initkatdetect <- function(genomicVariantsAnnotatedKat, segments, kataegisFoci){

    # store some information of the kataegis detection analysis
    info <- base::list(
        sampleName = base::as.character(base::unique(genomicVariantsAnnotatedKat@sampleNames)),
        totalGenomicVariants = base::length(genomicVariantsAnnotatedKat),
        totalKataegisFoci = base::length(kataegisFoci),
        totalVariantsInKataegisFoci = base::sum(genomicVariantsAnnotatedKat$putativeKataegis),
        version = base::as.character(utils::packageVersion('katdetectr')),
        date = base::date()
    )

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
