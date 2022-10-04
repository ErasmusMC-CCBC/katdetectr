#' @title detectKataegis
#' @description Detection of kataegis foci using changepoint detection.
#'
#' Changepoint detection is performed on the intermutation distance (IMD)
#' of the variants using the changepoint (killick 2014) package.
#'
#' Note that we recommend using the default parameters for the detection of kataegis.
#'
#' @param genomicVariants (\link[VariantAnnotation]{VRanges}, VCF or MAF): VRanges, or path to VCF or MAF file containing genomic variants.
#' @param minSizeKataegis (integer): Minimal number of variants required within a segment for classification as a kataegis foci.
#' @param maxMeanIMD (integer): Max. mean IMD within a segment for classification as a kataegis foci.
#' @param test.stat (character): Distribution that is fitted to the data (Exponential or Empirical). See \link[changepoint]{cpt.meanvar}.
#' @param penalty (character): Penalty used to guard against overfitting (BIC or Manual). See \link[changepoint]{cpt.meanvar}.
#' @param pen.value (integer): Only needed for manual penalty. See \link[changepoint]{cpt.meanvar}.
#' @param minseglen (integer): Min. size of segments (no. of variants).
#' @param BPPARAM (\link[BiocParallel]{BiocParallelParam}): Can be used for parallelization purposes.
#' @param aggregateRecords (logical): Aggregate multiple samples and treat as-if all records originate from a single sample.
#'
#' @examples
#' syntheticData <- generateSyntheticData()
#' kd <- detectKataegis(syntheticData)
#'
#' @return (KatDetect): Returns a KatDetect object including putative kataegis foci.
#' @author Daan Hazelaar
#' @author Job van Riet
#'
#' @references
#' Killick R, Eckley I (2014). “changepoint: An R package for changepoint analysis.” Journal of statistical software, 58(3), 1–19.
#'
#'
#' @export
detectKataegis <- function(genomicVariants, minSizeKataegis = 6, maxMeanIMD = 1000, test.stat = 'Exponential', penalty = 'BIC', pen.value = 0, minseglen = 2, BPPARAM = BiocParallel::SerialParam(), aggregateRecords = FALSE){

    validateInputdetectKataegis(genomicVariants, minSizeKataegis, maxMeanIMD, test.stat, penalty, pen.value , minseglen, BPPARAM, aggregateRecords)

    # Import, pre-process and annotate genomic variants. ----
    genomicVariantsImported <- .importGenomicVariants(x = genomicVariants, aggregateRecords)
    genomicVariantsProcessed <- .processGenomicVariants(genomicVariantsImported)
    genomicVariantsAnnotated <- .annotateGenomicVariants(genomicVariantsProcessed)

    # Changepoint detection. ----
    changepointsPerChromosome <- .performChangepointDetection(genomicVariantsAnnotated, test.stat, penalty, pen.value, minseglen, BPPARAM)

    # Annotate segments --------------------------------------------------------
    segments <- .annotateSegments(changepointsPerChromosome, genomicVariantsAnnotated)

    # Determine kataegis foci --------------------------------------------------
    kataegisFoci <- .determineKataegisFoci(segments, genomicVariantsAnnotated, minSizeKataegis, maxMeanIMD)

    # Add kataegis foci ID to genomicVariants
    genomicVariantsAnnotatedKat <- .addIDsToVariants(kataegisFoci, genomicVariantsAnnotated, changepointsPerChromosome)

    # Obtain relevant info for info slot
    info <- .getInfo(genomicVariantsAnnotatedKat, segments, kataegisFoci, minSizeKataegis, maxMeanIMD, test.stat, penalty, pen.value, minseglen, aggregateRecords)

    # Initialize KatDetect ---------------------------------------------------
    kd <- .initkatdetect(genomicVariantsAnnotatedKat, segments, kataegisFoci, info)

    return(kd)
}

