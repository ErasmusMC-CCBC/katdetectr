
.annotateSegments <- function(changepoints, genomicVariantsAnnotated){

    # determine for each variant in which segment it belongs based on the detected changepoints in the IMD
    segmentIDs <- .determineSegmentID(changepoints)
    # determine the segments and calculate additional info for each segment
    segments <- determineSegments(genomicVariantsAnnotated, segmentIDs)

    return(segments)
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

determineSegments <- function(genomicVariantsAnnotated, segmentIDs){

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
        dplyr::select(!diff) |>
        GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)

    return(segments)
}
