.determineKataegisFoci <- function(segments, genomicVariantsAnnotated, minSizeKataegis, IMDcutoffValues){

    kataegisFoci <- segments |>
        .determineKataegisSegments(IMDcutoffValues = IMDcutoffValues) |>
        .mergeKataegisSegments(minSizeKataegis = minSizeKataegis) |>
        .annotateKataegisSegments(genomicVariantsAnnotated = genomicVariantsAnnotated)

    return(kataegisFoci)
}

.determineKataegisSegments <- function(segments, IMDcutoffValues){

    selectedSegments <- segments |>
        tibble::as_tibble() |>
        dplyr::mutate(IMDcutoff = {{IMDcutoffValues}}) |>
        dplyr::filter(.data$meanIMD <= .data$IMDcutoff) |>
        dplyr::group_by(.data$seqnames)

    return(selectedSegments)
}

# function for annotating kataegis foci also merges segments that are in a single kataegis foci.
.determinefociID <- function(segmentIDs){

    nFoci <- base::sum(c(1, base::diff(segmentIDs)) != 1) + 1
    nSegmentsInFoci <- base::diff(c(1, base::which(c(1, base::diff(segmentIDs)) != 1), base::length(segmentIDs) + 1))
    fociID <- base::rep(base::seq_len(nFoci), nSegmentsInFoci)

    return(fociID)
}

.mergeKataegisSegments <- function(kataegisSegments, minSizeKataegis){

    # when no kataegis foci are present in the segments
    if(base::nrow(kataegisSegments) == 0){
        kataegisFoci <- tibble::tibble()
    } else {
        kataegisFoci <- kataegisSegments |>
            dplyr::mutate(
                fociID = .determinefociID(.data$segmentID)
            ) |>
            dplyr::group_by(.data$seqnames, .data$fociID) |>
            dplyr::summarise(
                .groups = "keep",
                seqnames = base::unique(.data$seqnames),
                start = base::min(.data$start),
                end = base::max(.data$end),
                sampleNames = base::unique(.data$sampleNames),
                totalVariants = base::sum(.data$totalVariants),
                firstVariantID = base::min(.data$firstVariantID),
                lastVariantID = base::max(.data$lastVariantID),
                meanIMD = base::mean(.data$meanIMD),
                IMDcutoff = base::unique(.data$IMDcutoff)
            ) |>
            dplyr::ungroup() |>
            dplyr::filter(.data$totalVariants >= minSizeKataegis - 1)
    }

    return(kataegisFoci)
}

.annotateKataegisSegments <- function(kataegisFoci, genomicVariantsAnnotated){

    # when still no kataegis foci are detected after merging of segments
    if(base::nrow(kataegisFoci) == 0){
        kataegisFociAnnotated <- GenomicRanges::GRanges()
    }else{
        kataegisFociAnnotated <- kataegisFoci |>
            # remove old fociID
            dplyr::select(!.data$fociID) |>
            # add new correct fociID number based on rowindex
            tibble::rowid_to_column("fociID") |>
            dplyr::rowwise() |>
            # manually add the last variants to the detected kataegis foci
            dplyr::mutate(
                firstVariantOfSeqname = base::min(genomicVariantsAnnotated$variantID[as.logical(GenomicRanges::seqnames(genomicVariantsAnnotated) == .data$seqnames)]),
                firstVariantID = base::ifelse(.data$firstVariantID != .data$firstVariantOfSeqname, .data$firstVariantID - 1, .data$firstVariantID),
                totalVariants = .data$lastVariantID - .data$firstVariantID + 1
            ) |>
            dplyr::ungroup() |>
            dplyr::select(!.data$firstVariantOfSeqname) |>
            # update start of kataegis foci and add column with sample name
            dplyr::mutate(
                start = GenomicRanges::start(genomicVariantsAnnotated[.data$firstVariantID])
            ) |>
            # convert to granges
            GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)
    }

    return(kataegisFociAnnotated)
}
