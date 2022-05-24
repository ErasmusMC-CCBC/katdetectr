# Accessor methods --------------------------------------------------------

#' @title Retrieve kateagis foci from a KatDetect object.
#' @inheritParams getGenomicVariants
#'
#' @examples
#' x <- new('KatDetect', genomicVariants = VariantAnnotation::VRanges(), info = list())
#' getKataegisFoci(x)
#'
#' @return (GRanges): Returns a GRanges with annotated kataegis foci.
#' @rdname getKataegisFoci
#' @export
#'
setGeneric('getKataegisFoci', function(x) standardGeneric("getKataegisFoci"))

#' @title Retrieve genomic variants from KatDetect object.
#' @param x (KatDetect): KatDetect object.
#'
#' @examples
#' x <- new('KatDetect', genomicVariants = VariantAnnotation::VRanges(), info = list())
#' getGenomicVariants(x)
#'
#' @return (VRanges): Returns a VRanges with annotated genomic variants.
#' @rdname getGenomicVariants
#' @export
#'
setGeneric('getGenomicVariants', function(x) standardGeneric("getGenomicVariants"))

#' @title Retrieve segments from KatDetect a object.
#' @inheritParams getGenomicVariants
#'
#' @examples
#' x <- new('KatDetect', genomicVariants = VariantAnnotation::VRanges(), info = list())
#' getSegments(x)
#'
#' @return (GRanges): Returns a GRanges with annotated segments.
#' @rdname getSegments
#' @export
#'
setGeneric('getSegments', function(x) standardGeneric("getSegments"))

#' @title Retrieve model parameters from a KatDetect object.
#' @inheritParams getGenomicVariants
#'
#' @examples
#' x <- new('KatDetect', genomicVariants = VariantAnnotation::VRanges(), info = list())
#' getInfo(x)
#'
#' @return (list): Returns a list with all model parameters used for kataegis detection.
#' @rdname getInfo
#' @export
#'
setGeneric('getInfo', function(x) standardGeneric("getInfo"))
