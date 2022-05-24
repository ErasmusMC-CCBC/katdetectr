.addKataegisIDtoVariants <- function(kataegisFoci, genomicVariantsAnnotated){

    genomicVariantsAnnotatedKat <- .addKataegisIDtoVariantsPerSample(kataegisFoci, genomicVariantsAnnotated)

    return(genomicVariantsAnnotatedKat)
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
