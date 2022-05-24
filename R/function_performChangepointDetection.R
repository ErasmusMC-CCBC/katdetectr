
# this function fits a PCF model to the intermutations distances of each chromosome of each sample. The detected changepoints are returned
.performChangepointDetection <- function(genomicVariantsAnnotated, test.stat, penalty, pen.value, minseglen, bpworkers){

    # Split on chromosome and obtain IMD.
    perChromosomeIMD <- .getIMD(genomicVariantsAnnotated)

    # Fit a PCF model to the IMD data and return the changepoints
    changepoints <- .runCD(perChromosomeIMD, test.stat, penalty, pen.value, minseglen, bpworkers)

    return(changepoints)
}

# Helper - Retrieve IMD per chromosome. ----
.getIMD <- function(x){

    genomicVariantsPerChromosome <- base::split(x, GenomeInfoDb::seqnames(x))

    IMDs <- base::lapply(genomicVariantsPerChromosome, function(chromosome){

        # Remove the first variant of the sequence (which has no IMD).
        IMDraw <- stats::na.omit(chromosome$IMD)

        # Change 0 to 0.001 as 0 cannot be sampled from an exponential distribution
        IMD <- base::ifelse(IMDraw == 0, 0.001, IMDraw)

        return(IMD)
    })

    return(IMDs)
}

# Helper - Perform changepoint detection. ----
.runCD <- function(perChromosomeIMD, test.stat, penalty, pen.value, minseglen, bpworkers){

    changepointsPerChromosome <- BiocParallel::bplapply(perChromosomeIMD, function(IMD, .test.stat = test.stat, .penalty = penalty, .pen.value = pen.value, .minseglen = minseglen, .bpworkers = bpworkers){

        # At least 4 observations (variants) are needed for changepoint analysis.
        if(base::length(IMD) >= 4){
            if(.test.stat == 'Exponential'){
                changepointsChromosome <- changepoint::cpt.meanvar(
                    data = IMD,
                    penalty = .penalty,
                    pen.value = .pen.value,
                    method = 'PELT',
                    test.stat = 'Exponential',
                    class = FALSE,
                    minseglen = .minseglen
                )
            }

            if(.test.stat == 'Empirical'){
                changepointsChromosome <- changepoint.np::cpt.np(
                    data = IMD,
                    penalty = .penalty,
                    pen.value = .pen.value,
                    method = 'PELT',
                    test.stat = 'empirical_distribution',
                    class = FALSE,
                    minseglen = .minseglen
                )
            }

            # Add pseudo-count of 1 to all changepoints as the first variant was not included as it had no 5' IMD.
            # Add 0 as the first changepoint.
            changepointsChromosome <- c(0, (changepointsChromosome + 1))

        }else{
            # For <4 observations, return first and last variant to set a single segment.
            changepointsChromosome <- c(0, (base::length(IMD) + 1))
        }

        return(changepointsChromosome)
    },
    BPPARAM = BiocParallel::MulticoreParam(workers = bpworkers))

    return(changepointsPerChromosome)
}
