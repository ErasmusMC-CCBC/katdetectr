.addIDsToVariants <- function(kataegisFoci, genomicVariantsAnnotated, changepointsPerChromosome){

    # obtain results from changepoint analysis
    changepoints <- .getChangepoints(changepointsPerChromosome)
    # determine for each variant in which segment it belongs based on the detected changepoints in the IMD
    segmentIDs <- .determineSegmentID(changepoints)

    # add column which specifies the segment ID each variant belong to
    genomicVariantsAnnotatedSeg <- genomicVariantsAnnotated |>
        .addSegmentsIDtovariants(segmentIDs) |>
        GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE) |>
        VariantAnnotation::makeVRangesFromGRanges()

    # add coumn which specifies if a variant lies within a kataegis foci
    genomicVariantsAnnotatedSegKat <- .addKataegisIDtoVariantsPerSample(kataegisFoci, genomicVariantsAnnotatedSeg)

    return(genomicVariantsAnnotatedSegKat)
}

.addKataegisIDtoVariantsPerSample <- function(kataegisFoci, genomicVariantsAnnotated){

    if(base::length(kataegisFoci) > 0){
        # label variants that are in a detected kataegis foc
        kataegisVariants <- IRanges::subsetByOverlaps(genomicVariantsAnnotated, kataegisFoci)
        kataegisVariants$putativeKataegis <- TRUE

        # label variants that are not in a detected kataegis foci
        noKataegisVariants <- genomicVariantsAnnotated[-kataegisVariants$variantID]
        noKataegisVariants$putativeKataegis <- FALSE

        # update variants
        genomicVariantsAnnotated <- GenomicRanges::sort(c(kataegisVariants, noKataegisVariants))
    } else {
        genomicVariantsAnnotated$putativeKataegis <- FALSE
    }

    return(genomicVariantsAnnotated)
}
