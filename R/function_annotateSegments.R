

.annotateSegments <- function(changepointsPerChromosome, genomicVariantsAnnotated){

    # obtain results from changepoint analysis
    changepoints <- .getChangepoints(changepointsPerChromosome)
    rates <- .getRates(changepointsPerChromosome)

    # determine for each variant in which segment it belongs based on the detected changepoints in the IMD
    segmentIDs <- .determineSegmentID(changepoints)

    # determine the segments and calculate additional info for each segment
    segments <- determineSegments(genomicVariantsAnnotated, segmentIDs, rates)

    return(segments)
}


.getChangepoints <- function(changepointsPerChromosome){

    changepoints <- lapply(changepointsPerChromosome, function(chromosome){chromosome$changepointsChromosome})

    return(changepoints)
}

.getRates <- function(changepointsPerChromosome){

    rates <- lapply(changepointsPerChromosome, function(chromosome){chromosome$rateChromosome}) |>
        unlist(use.names = FALSE)

    return(rates)
}

.determineSegmentID <- function(changepoints){

    segmentIDs <- base::lapply(changepoints, .determineSegmentIDperChr) |>
        base::unlist(use.names = FALSE)

    return(segmentIDs)
}

.determineSegmentIDperChr <- function(changepointsPerChr){

    nSegments <- base::length(changepointsPerChr) - 1
    difference <- base::diff(changepointsPerChr)
    segmentID <- base::unname(base::rep(base::seq_len(nSegments), difference))

    return(segmentID)
}

.addSegmentsIDtovariants <- function(variants, segmentIDs){

    variants <- variants |>
        tibble::as_tibble() |>
        dplyr::mutate(segmentID = base::unlist(segmentIDs, use.names = FALSE))

    return(variants)
}

determineSegments <- function(genomicVariantsAnnotated, segmentIDs, rates){

    segments <- genomicVariantsAnnotated |>
        .addSegmentsIDtovariants(segmentIDs = segmentIDs) |>
        dplyr::group_by(.data$sampleNames, .data$seqnames, .data$segmentID) |>
        dplyr::summarise(
            .groups = 'keep',
            totalVariants = base::sum(dplyr::n()),
            firstVariantID = base::min(.data$variantID),
            lastVariantID = base::max(.data$variantID),
            start = base::min(.data$start),
            end = base::max(.data$end),
            meanIMD = base::mean(.data$IMD, na.rm = TRUE),
            sampleNames = base::as.character(base::unique(.data$sampleNames))
        ) |>
        # replace mean IMD to NA if NAN
        dplyr::na_if('NaN') |>
        dplyr::group_by(.data$sampleNames, .data$seqnames) |>
        # update start column to make sure segments border each other.
        dplyr::mutate(
            diff = base::ifelse(!base::is.na(dplyr::lag(.data$end)), .data$start - dplyr::lag(.data$end), 0),
            start = .data$start - .data$diff + 1
        ) |>
        dplyr::ungroup() |>
        dplyr::mutate(mutationRate = rates) |>
        dplyr::select(!diff) |>
        GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)

    return(segments)
}
