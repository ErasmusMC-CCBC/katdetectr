# Internal - Annotate genomic variants. ----
.annotateGenomicVariants <- function(genomicVariants){

    genomicVariantsAnnotated <- .determineIMD(genomicVariants)

    return(genomicVariantsAnnotated)
}

# Helper - Calculate the 5' intermutation distance (IMD) ----
.determineIMD <- function(genomicVariants){

    genomicVariants$IMD <- c(base::as.numeric(NA), GenomicRanges::distance(genomicVariants[-base::length(genomicVariants)], genomicVariants[-1]))

    return(genomicVariants)
}
