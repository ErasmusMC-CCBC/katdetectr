.annotateSegments <- function(changepointsPerChromosome, genomicVariantsAnnotated, refSeq) {
    # obtain results from changepoint analysis
    changepoints <- .getChangepoints(changepointsPerChromosome)
    rates <- .getRates(changepointsPerChromosome)

    # determine for each variant in which segment it belongs based on the detected changepoints in the IMD
    segmentIDs <- .determineSegmentID(changepoints)

    # determine the segments and calculate additional info for each segment
    segments <- determineSegments(genomicVariantsAnnotated, segmentIDs, rates, refSeq)

    return(segments)
}


.getChangepoints <- function(changepointsPerChromosome) {
    changepoints <- lapply(changepointsPerChromosome, function(chromosome) {
        chromosome$changepointsChromosome
    })

    return(changepoints)
}

.getRates <- function(changepointsPerChromosome) {
    rates <- lapply(changepointsPerChromosome, function(chromosome) {
        chromosome$rateChromosome
    }) |>
        unlist(use.names = FALSE)

    return(rates)
}

.determineSegmentID <- function(changepoints) {
    segmentIDs <- base::lapply(changepoints, .determineSegmentIDperChr) |>
        base::unlist(use.names = FALSE)

    return(segmentIDs)
}

.determineSegmentIDperChr <- function(changepointsPerChr) {
    nSegments <- base::length(changepointsPerChr) - 1
    difference <- base::diff(changepointsPerChr)
    segmentID <- base::unname(base::rep(base::seq_len(nSegments), difference))

    return(segmentID)
}

.addSegmentsIDtovariants <- function(variants, segmentIDs) {
    variants <- variants |>
        tibble::as_tibble() |>
        dplyr::mutate(segmentID = base::unlist(segmentIDs, use.names = FALSE))

    return(variants)
}

.addEmptySegments <- function(segments, refSeq) {
    if (all(refSeq %in% c("hg19", "hg38"))) {
        chromosomeNames <- paste0("chr", c(1:22, "X", "Y", "M"))
    } else {
        chromosomeNames <- levels(segments$seqnames)
    }
    chromosomesWithMutations <- unique(segments$seqnames)

    # select chromosomes without any mutations in the sample
    emptyChromosomes <- chromosomeNames[!chromosomeNames %in% chromosomesWithMutations]

    if (length(emptyChromosomes) != 0) {
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
                end = getChromosomeLength(chromosome = base::unique(seqnames), refSeq)
            ) |>
            dplyr::ungroup() |>
            dplyr::bind_rows(segments)
    } else {
        segmentsInclEmpty <- segments
    }

    return(segmentsInclEmpty)
}

determineSegments <- function(genomicVariantsAnnotated, segmentIDs, rates, refSeq) {
    segments <- genomicVariantsAnnotated |>
        .addSegmentsIDtovariants(segmentIDs = segmentIDs) |>
        dplyr::group_by(sampleNames, seqnames, segmentID) |>
        dplyr::summarise(
            .groups = "keep",
            totalVariants = base::sum(dplyr::n()),
            firstVariantID = base::min(variantID),
            lastVariantID = base::max(variantID),
            end = base::max(start),
            start = base::min(start),
            meanIMD = base::mean(IMD, na.rm = TRUE),
            sampleNames = base::as.character(base::unique(sampleNames))
        ) |>
        dplyr::group_by(sampleNames, seqnames) |>
        dplyr::mutate(
            # update start column to make sure segments border each other.
            diff = base::ifelse(!base::is.na(dplyr::lag(end)), start - dplyr::lag(end), 0),
            start = start - diff + 1,
            # make sure the first segment starts at the beginning of the sequence
            start = c(1, start[-1]),
            # make sure the last segment ends at the end of the sequence
            end = c(end[-base::length(end)], getChromosomeLength(chromosome = as.character(base::unique(seqnames)), refSeq))
        ) |>
        dplyr::ungroup() |>
        # add the mutation rate of the segments
        dplyr::mutate(
            mutationRate = rates
        ) |>
        dplyr::group_by(sampleNames, seqnames) |>
        # recalculate the mean IMD and rate of the last segment of each sequence (due to the added pseudo count during changepoint analysis)
        dplyr::mutate(
            meanIMD = c(
                meanIMD[-base::length(meanIMD)],
                (end[base::length(end)] - start[base::length(start)]) / totalVariants[base::length(totalVariants)]
            ),
            mutationRate = c(
                mutationRate[-base::length(mutationRate)],
                1 / meanIMD[base::length(meanIMD)]
            )
        ) |>
        dplyr::select(!diff) |>
        .addEmptySegments(refSeq) |>
        GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE) |>
        GenomicRanges::sort()

    GenomeInfoDb::seqlevels(segments) <- GenomeInfoDb::seqlevelsInUse(segments)

    return(segments)
}
