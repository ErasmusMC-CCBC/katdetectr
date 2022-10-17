# TODO make sure that the final segments span the entire genome. Now the segments only span from the first to the last variant

# this function fits a PCF model to the intermutations distances of each chromosome of each sample. The detected changepoints are returned
.performChangepointDetection <- function(genomicVariantsAnnotated, test.stat, penalty, pen.value, minseglen, BPPARAM){

    # Split on chromosome and obtain IMD.
    perChromosomeIMD <- .getIMD(genomicVariantsAnnotated)

    # Fit a PCF model to the IMD data and return the changepoints
    changepoints <- .runCD(perChromosomeIMD, test.stat, penalty, pen.value, minseglen, BPPARAM)

    return(changepoints)
}

# Helper - Retrieve IMD per chromosome. ----
.getIMD <- function(x){

    genomicVariantsPerChromosome <- base::split(x, GenomeInfoDb::seqnames(x))

    IMDs <- base::lapply(genomicVariantsPerChromosome, function(chromosome){

        IMD <- chromosome$IMD

        return(IMD)
    })

    return(IMDs)
}

#TODO add lookup table for hg38 and add logic for sequences that are not present in the lookup table

getChromosomeLength <- function(chromosome){

    # lookup table for calculating mutation rate of segments with very few variants
    # length op chromosome in bp according to BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
    chromosomeLengthTib <- tibble::tibble(
        chr1 = 249250621,
        chr2 = 243199373,
        chr3 = 198022430,
        chr4 = 191154276,
        chr5 = 180915260,
        chr6 = 171115067,
        chr7 = 159138663,
        chr8 = 146364022,
        chr9 = 141213431,
        chr10 = 135534747,
        chr11 = 135006516,
        chr12 = 133851895,
        chr13 = 115169878,
        chr14 = 107349540,
        chr15 = 102531392,
        chr16 = 90354753,
        chr17 = 81195210,
        chr18 = 78077248,
        chr19 = 59128983,
        chr20 = 63025520,
        chr21 = 48129895,
        chr22 = 51304566,
        chrX = 155270560,
        chrY = 59373566,
        chrM = 16571
    )

    chromosomeLength <- chromosomeLengthTib[[chromosome]]

    return(chromosomeLength)
}

# Helper - Perform changepoint detection. ----
.runCD <- function(perChromosomeIMD, test.stat, penalty, pen.value, minseglen, BPPARAM){

    # loop over the elements of perChromosomeIMD. Normal map is not possible as the names of the list must be available in the function body
    changepointsPerChromosome <- BiocParallel::bplapply(seq(perChromosomeIMD), function(i, .test.stat = test.stat, .penalty = penalty, .pen.value = pen.value, .minseglen = minseglen, .BPPARAM = BPPARAM){

        # At least 4 observations (IMDs) are needed for changepoint analysis.
        if(base::length(perChromosomeIMD[[i]]) >= 4){
            if(.test.stat == 'Exponential'){

                cptChromosome <- changepoint::cpt.meanvar(
                    data = perChromosomeIMD[[i]],
                    penalty = .penalty,
                    pen.value = .pen.value,
                    method = 'PELT',
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

            # Add pseudo-count of 1 to all changepoints as the first variant was not included as it had no 5' IMD.
            # Add 0 as the first changepoint.
            changepointsChromosome <- c(0, cptChromosome@cpts)
            rateChromosome <- changepoint::param.est(cptChromosome)$rate
        }else{
            # For <4 observations, return first and last variant to set a single segment.
            changepointsChromosome <- c(0, (base::length(perChromosomeIMD[[i]])))
            # For <4 observations, the rate must be calculated manually. rate = nVariants / length of chromosome
            rateChromosome <- length(perChromosomeIMD[[i]]) / getChromosomeLength(chromosome = names(perChromosomeIMD[i]))
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
