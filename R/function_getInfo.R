.getInfo <- function(genomicVariantsAnnotatedKat, segments, kataegisFoci, minSizeKataegis, maxMeanIMD, test.stat, penalty, pen.value, minseglen, aggregateRecords){

    parametersList <- base::list(
        minSizeKataegis = minSizeKataegis,
        maxMeanIMD = maxMeanIMD,
        test.stat = test.stat,
        penalty = penalty,
        pen.value = pen.value,
        minseglen = minseglen,
        aggregateRecords = aggregateRecords

    )

    # store some information of the kataegis detection analysis
    info <- base::list(
        sampleName = base::as.character(base::unique(VariantAnnotation::sampleNames(genomicVariantsAnnotatedKat))),
        totalGenomicVariants = base::length(genomicVariantsAnnotatedKat),
        totalKataegisFoci = base::length(kataegisFoci),
        totalVariantsInKataegisFoci = base::sum(genomicVariantsAnnotatedKat$putativeKataegis),
        version = base::as.character(utils::packageVersion('katdetectr')),
        date = base::date(),
        parameters = parametersList
    )

    return(info)
}
