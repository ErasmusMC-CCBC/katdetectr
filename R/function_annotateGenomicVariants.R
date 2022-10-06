# Internal - Annotate genomic variants. ----
.annotateGenomicVariants <- function(genomicVariants){

    genomicVariantsSplitonChr <- base::split(genomicVariants, GenomeInfoDb::seqnames(genomicVariants))

    genomicVariantsAnnotatedList <- base::lapply(genomicVariantsSplitonChr, .determineIMD)

    genomicVariantsAnnotated <- plyranges::bind_ranges(genomicVariantsAnnotatedList)

    return(genomicVariantsAnnotated)
}

# Helper - Calculate the 5' intermutation distance (IMD) ----
.determineIMD <- function(genomicVariants){

    genomicVariants$IMD <- c(GenomicRanges::start(genomicVariants)[1],
                             GenomicRanges::start(genomicVariants)[-1] - GenomicRanges::start(genomicVariants)[-base::length(genomicVariants)]
                             )

    return(genomicVariants)
}
