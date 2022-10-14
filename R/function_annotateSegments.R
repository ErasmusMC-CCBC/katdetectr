

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

.addEmptySegments <- function(segments){

    chromosomeNames <- levels(segments$seqnames)
    chromosomesWithMutations <- unique(segments$seqnames)

    # select chromosomes without any mutations in the sample
    emptyChromosomes <- chromosomeNames[!chromosomeNames %in% chromosomesWithMutations]

    if(length(emptyChromosomes) != 0){

        # construct tibble for the empty chromosomes and fill in columns accordingly
        segmentsInclEmpty <- tibble::tibble(
            seqnames = emptyChromosomes,
            segmentID = 1,
            totalVariants = 0,
            firstVariantID = NA,
            lastVariantID = NA,
            start = 1,
            meanIMD = NA,
            mutationRate = 0,
            sampleNames = base::as.character(base::unique(segments$sampleNames))
        ) |>
            dplyr::rowwise() |>
            dplyr::mutate(
                end = getChromosomeLength(chromosome = base::unique(.data$seqnames))
            ) |>
            dplyr::ungroup() |>
            dplyr::bind_rows(segments)

    } else {
        segmentsInclEmpty <- segments
    }

    return(segmentsInclEmpty)
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
            end = base::max(.data$start),
            start = base::min(.data$start),
            meanIMD = base::mean(.data$IMD, na.rm = TRUE),
            sampleNames = base::as.character(base::unique(.data$sampleNames))
        ) |>
        dplyr::group_by(.data$sampleNames, .data$seqnames) |>
        dplyr::mutate(
            # update start column to make sure segments border each other.
            diff = base::ifelse(!base::is.na(dplyr::lag(.data$end)), .data$start - dplyr::lag(.data$end), 0),
            start = .data$start - .data$diff + 1,
            # make sure the first segment starts at the beginning of the sequence
            start = c(1, .data$start[-1])
        ) |>
        dplyr::ungroup() |>
        dplyr::mutate(mutationRate = rates) |>
        dplyr::select(!diff) |>
        .addEmptySegments() |>
        GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE) |>
        GenomicRanges::sort()

    return(segments)
}
