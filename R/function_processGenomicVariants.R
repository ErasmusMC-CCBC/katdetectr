# Internal - Clean-up genomic variants prior to analysis. ----
.processGenomicVariants <- function(genomicVariantsImported) {
    # Reduce overlapping variants.
    genomicVariantsProcessed <- .reduceOverlappingVariants(genomicVariantsImported)

    return(genomicVariantsProcessed)
}


.reduceOverlappingVariants <- function(genomicVariantsImported) {
    # Reduce overlapping mutations to a single range.
    variantsReduced <- GenomicRanges::reduce(genomicVariantsImported, min.gapwidth = 0, ignore.strand = TRUE, with.revmap = TRUE)

    # Select variants that (partially) overlap with another variant.
    overlappingVariants <- variantsReduced[S4Vectors::elementNROWS(variantsReduced$revmap) != 1, ]

    # Select variants that do not overlap with other variants.
    nonOverlappingVariants <- genomicVariantsImported[-base::unlist(overlappingVariants$revmap)]

    # If there are no overlapping variants, return original.
    # Otherwise, add a column that specifies which variants overlapped and were reduced.
    if (base::length(overlappingVariants) == 0) {
        variantsUnique <- GenomicRanges::sort(genomicVariantsImported)
    } else {
        # Add necessary columns for VRanges
        variantsUnique <- overlappingVariants |>
            tibble::as_tibble() |>
            dplyr::mutate(
                ref = "X",
                alt = "XX",
                sampleNames = base::unique(VariantAnnotation::sampleNames(genomicVariantsImported))
            ) |>
            # join non with  overlapping variants
            dplyr::bind_rows(tibble::as_tibble(nonOverlappingVariants)) |>
            # coerce to VRanges
            GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE) |>
            GenomicRanges::sort() |>
            VariantAnnotation::makeVRangesFromGRanges()
    }

    # Add a unique id number for each variant.
    variantsUnique$variantID <- base::seq_len(base::length(variantsUnique))

    # Store which variants from the original data belong to each variant in the reduced data (i.e. the mutation events in a range)
    variantsUnique$revmap <- variantsReduced$revmap

    variantsUnique <- GenomicRanges::sort(variantsUnique)

    # drop unused levels
    GenomeInfoDb::seqlevels(variantsUnique) <- GenomeInfoDb::seqlevelsInUse(variantsUnique)

    return(variantsUnique)
}
