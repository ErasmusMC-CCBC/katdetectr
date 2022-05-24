#' @title detectKataegis
#' @description Detection of kataegis foci using changepoint detection.
#'
#' Changepoint detection is performed on the intermutation distance (IMD)
#' of the variants using the changepoint \insertCite{killick2014changepoint}{katdetectr} package. 
#' 
#' Note that we recommend using the default parameters for the detection of kataegis.
#'
#' @param genomicVariants (\link[VariantAnnotation]{VRangesList}, VCF or MAF): VRanges, or path to VCF or MAF file containing genomic variants.
#' @param minSizeKataegis (integer): Minimal no. of variants required within a segment for classification as a kataegis foci.
#' @param maxMeanIMD (integer): Max. mean IMD within a segment for classification as a kataegis foci.
#' @param test.stat (character): Distribution that is fitted to the data (Exponential or Empirical). See \link[changepoint]{cpt.meanvar}.
#' @param penalty (character): Penalty used to guard against overfitting (BIC or Manual). See \link[changepoint]{cpt.meanvar}.
#' @param pen.value (integer): Only needed for manual penalty. See \link[changepoint]{cpt.meanvar}.
#' @param minseglen (integer): Min. size of segments (no. of variants).
#' @param bpworkers (integer): Number of cores to be used in parallelisation.
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
#' \insertRef{killick2014changepoint}{katdetectr}
#'
#' @importFrom Rdpack reprompt
#' @import BiocStyle
#' @export
detectKataegis <- function(genomicVariants, minSizeKataegis = 6, maxMeanIMD = 1000, test.stat = 'Exponential', penalty = 'BIC', pen.value = 0, minseglen = 2, bpworkers = 1, aggregateRecords = FALSE){

    validateInputdetectKataegis(genomicVariants, minSizeKataegis, maxMeanIMD, test.stat, penalty, pen.value , minseglen, bpworkers, aggregateRecords)

    # Import, pre-process and annotate genomic variants. ----
    genomicVariantsImported <- .importGenomicVariants(x = genomicVariants, aggregateRecords = aggregateRecords)
    genomicVariantsProcessed <- .processGenomicVariants(genomicVariantsImported)
    genomicVariantsAnnotated <- .annotateGenomicVariants(genomicVariantsProcessed)

    # Changepoint detection. ----
    changepointsPerChromosome <- .performChangepointDetection(genomicVariantsAnnotated, test.stat, penalty, pen.value, minseglen, bpworkers = bpworkers)

    # Annotate segments --------------------------------------------------------
    segments <- .annotateSegments(changepointsPerChromosome, genomicVariantsAnnotated)

    # Determine kataegis foci --------------------------------------------------
    kataegisFoci <- .determineKataegisFoci(segments, genomicVariantsAnnotated, minSizeKataegis, maxMeanIMD)

    # Add kataegis foci ID to genomicVariants
    genomicVariantsAnnotatedKat <- .addKataegisIDtoVariants(kataegisFoci, genomicVariantsAnnotated)

    # Initialize KatDetect ---------------------------------------------------
    kd <- .initkatdetect(genomicVariantsAnnotatedKat, segments, kataegisFoci)

    return(kd)
}

