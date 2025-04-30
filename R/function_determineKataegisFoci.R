.determineKataegisFoci <- function(segments, genomicVariantsAnnotated, minSizeKataegis, IMDcutoffValues) {
    kataegisFoci <- segments |>
        .determineKataegisSegments(IMDcutoffValues = IMDcutoffValues) |>
        .mergeKataegisSegments(minSizeKataegis = minSizeKataegis) |>
        .annotateKataegisSegments(genomicVariantsAnnotated = genomicVariantsAnnotated)

    return(kataegisFoci)
}

.determineKataegisSegments <- function(segments, IMDcutoffValues) {
    selectedSegments <- segments |>
        tibble::as_tibble() |>
        dplyr::mutate(IMDcutoff = {{ IMDcutoffValues }}) |>
        dplyr::filter(meanIMD <= IMDcutoff) |>
        dplyr::group_by(seqnames)

    return(selectedSegments)
}

# function for annotating kataegis foci also merges segments that are in a single kataegis foci.
.determinefociID <- function(segmentIDs) {
    nFoci <- base::sum(c(1, base::diff(segmentIDs)) != 1) + 1
    nSegmentsInFoci <- base::diff(c(1, base::which(c(1, base::diff(segmentIDs)) != 1), base::length(segmentIDs) + 1))
    fociID <- base::rep(base::seq_len(nFoci), nSegmentsInFoci)

    return(fociID)
}

.mergeKataegisSegments <- function(kataegisSegments, minSizeKataegis) {
    # when no kataegis foci are present in the segments
    if (base::nrow(kataegisSegments) == 0) {
        kataegisFoci <- tibble::tibble()
    } else {
        kataegisFoci <- kataegisSegments |>
            dplyr::mutate(
                fociID = .determinefociID(segmentID)
            ) |>
            dplyr::group_by(seqnames, fociID) |>
            dplyr::summarise(
                .groups = "keep",
                seqnames = base::unique(seqnames),
                start = base::min(start),
                end = base::max(end),
                sampleNames = base::unique(sampleNames),
                totalVariants = base::sum(totalVariants),
                firstVariantID = base::min(firstVariantID),
                lastVariantID = base::max(lastVariantID),
                meanIMD = base::mean(meanIMD),
                IMDcutoff = base::unique(IMDcutoff)
            ) |>
            dplyr::ungroup() |>
            dplyr::filter(totalVariants >= minSizeKataegis - 1)
    }

    return(kataegisFoci)
}

.annotateKataegisSegments <- function(kataegisFoci, genomicVariantsAnnotated) {
    # when still no kataegis foci are detected after merging of segments
    if (base::nrow(kataegisFoci) == 0) {
        kataegisFociAnnotated <- GenomicRanges::GRanges()
    } else {
        kataegisFociAnnotated <- kataegisFoci |>
            # remove old fociID
            dplyr::select(!fociID) |>
            # add new correct fociID number based on rowindex
            tibble::rowid_to_column("fociID") |>
            dplyr::rowwise() |>
            # manually add the last variants to the detected kataegis foci
            dplyr::mutate(
                firstVariantOfSeqname = base::min(genomicVariantsAnnotated$variantID[as.logical(as.character(GenomicRanges::seqnames(genomicVariantsAnnotated)) == as.character(seqnames))]),
                firstVariantID = base::ifelse(firstVariantID != firstVariantOfSeqname, firstVariantID - 1, firstVariantID),
                totalVariants = lastVariantID - firstVariantID + 1
            ) |>
            dplyr::ungroup() |>
            dplyr::select(!firstVariantOfSeqname) |>
            # update start of kataegis foci and add column with sample name
            dplyr::mutate(
                start = GenomicRanges::start(genomicVariantsAnnotated[firstVariantID])
            ) |>
            # convert to granges
            GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)
    }

    return(kataegisFociAnnotated)
}
