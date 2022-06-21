# Accessor methods --------------------------------------------------------

#' @title Retrieve kateagis foci from a KatDetect object.
#' @inheritParams getGenomicVariants
#'
#' @examples
#' kd <- constructKatdetect()
#' getKataegisFoci(kd)
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
#' kd <- constructKatdetect()
#' getGenomicVariants(kd)
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
#' kd <- constructKatdetect()
#' getSegments(kd)
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
#' kd <- constructKatdetect()
#' getInfo(kd)
#'
#' @return (list): Returns a list with all model parameters used for kataegis detection.
#' @rdname getInfo
#' @export
#'
setGeneric('getInfo', function(x) standardGeneric("getInfo"))
