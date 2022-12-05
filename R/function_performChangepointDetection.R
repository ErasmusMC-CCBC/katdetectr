
# this function fits a model to the intermutations distances of each chromosome of each sample. The detected changepoints are returned
.performChangepointDetection <- function(genomicVariantsAnnotated, refSeq, test.stat, penalty, pen.value, method, minseglen, BPPARAM){

    # Split on chromosome and obtain IMD.
    perChromosomeIMD <- .getIMD(genomicVariantsAnnotated, refSeq)

    # Fit model to the IMD data and return the changepoints
    changepoints <- .runCD(refSeq, perChromosomeIMD, test.stat, penalty, pen.value, method, minseglen, BPPARAM)

    return(changepoints)
}

# Helper - Retrieve IMD per chromosome. ----
.getIMD <- function(genomicVariantsAnnotated, refSeq){

    genomicVariantsPerChromosome <- base::split(genomicVariantsAnnotated, GenomeInfoDb::seqnames(genomicVariantsAnnotated))

    IMDs <- base::lapply(genomicVariantsPerChromosome, function(chromosome){

        IMD <- chromosome$IMD

        # In order to consider the whole (finite) DNA sequence I add one pseudo IMD which is the distance from the last variant to the end of the DNA sequence
        # This is necessary to make sure that the sum of all the rates detected in changepoint analysis equal the mutation rate of the entire chromosome
        chromosomeLength <- chromosome@seqnames |>
            as.character() |>
            unique() |>
            getChromosomeLength(refSeq)

        distanceToEndChr <- chromosomeLength - sum(IMD)

        if(distanceToEndChr < 0){base::stop("Make sure to provide the correct reference genome using the refSeq argument")}

        IMDfullSequence <- c(IMD, distanceToEndChr)

        return(IMDfullSequence)
    })

    return(IMDs)
}



# Helper - Perform changepoint detection. ----
.runCD <- function(refSeq, perChromosomeIMD, test.stat, penalty, pen.value, method, minseglen, BPPARAM){

    # loop over the elements of perChromosomeIMD. Normal map is not possible as the names of the list must be available in the function body
    changepointsPerChromosome <- BiocParallel::bplapply(seq(perChromosomeIMD), function(i, .refSeq = refSeq, .test.stat = test.stat, .penalty = penalty, .pen.value = pen.value, .method = method, .minseglen = minseglen, .BPPARAM = BPPARAM){

        # At least 4 observations (IMDs) are needed for changepoint analysis.
        if(base::length(perChromosomeIMD[[i]]) >= 4){
            if(.test.stat == 'Exponential'){

                cptChromosome <- changepoint::cpt.meanvar(
                    data = perChromosomeIMD[[i]],
                    penalty = .penalty,
                    pen.value = .pen.value,
                    method = .method,
                    test.stat = 'Exponential',
                    class = TRUE,
                    minseglen = .minseglen
                )
            }

            if(.test.stat == 'Empirical'){
                cptChromosome <- changepoint.np::cpt.np(
                    data = perChromosomeIMD[[i]],
                    penalty = .penalty,
                    pen.value = .pen.value,
                    method = 'PELT',
                    test.stat = 'empirical_distribution',
                    class = TRUE,
                    minseglen = .minseglen
                )
            }

            # Add 0 as the first changepoint. And substract 1 from the last changepoint due to the added pseudo IMD
            changepointsChromosome <- c(0,
                                        cptChromosome@cpts[-base::length(cptChromosome@cpts)],
                                        cptChromosome@cpts[base::length(cptChromosome@cpts)] - 1)
            rateChromosome <- changepoint::param.est(cptChromosome)$rate
        }else{
            # For <4 observations, return first and last variant as changepoints to set a single segment.
            # ubstract 1 from the last changepoint due to the added pseudo IMD
            changepointsChromosome <- c(0, base::length(perChromosomeIMD[[i]]) - 1)
            # For <4 observations, the rate must be calculated manually. rate = nVariants / length of chromosome
            rateChromosome <- length(perChromosomeIMD[[i]]) / getChromosomeLength(chromosome = names(perChromosomeIMD[i]), .refSeq)
        }

        resultsChromosome <- list(
            changepointsChromosome = changepointsChromosome,
            rateChromosome = rateChromosome
        )
        return(resultsChromosome)
    },
    BPPARAM = BPPARAM)

    # add back chromosome names as the names of the list
    names(changepointsPerChromosome) <- names(perChromosomeIMD)

    return(changepointsPerChromosome)
}
