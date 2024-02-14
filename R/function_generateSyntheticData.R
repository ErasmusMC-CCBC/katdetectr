#' @title Generate genomic variants with pre-defined kataegis foci.
#' @description This function generates a synthetic dataset (VRanges) containing background genomic variants
#' and a no. of desired interjected kataegis foci.
#'
#' @param genome (\link[BSgenome]{BSgenome}): The genome (DNA) which will be sampled for genomic variants.
#' @param nBackgroundVariants (integer): The no. of generated background genomic variants.
#' @param seqnames (character): The sequences on which genomic variants will be sampled. If NULL, then all human autosomes and sex chromosomes will be used.
#' @param probMutationType (numeric): The probability of a generated variant being an SNV, MNV, Deletion or Insertion, respectively.
#' @param nKataegisFoci (integer): The no. of generated and interjected kataegis foci.
#' @param nKataegisVariants (integer): The no. of genomic variants within each kataegis foci.
#' @param expectedIMD (integer): The expected mean intermutational distance (IMD) of the generated kataegis foci.
#' @param sampleName (character): The name of the sample
#' @param removeValidationColumns (logical): Include columns with extra information regarding mutation sampling?
#'
#' @examples
#' syntheticData1 <- generateSyntheticData()
#'
#' syntheticData2 <- generateSyntheticData(
#'     genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
#'     nBackgroundVariants = 75,
#'     seqnames = c("chr1", "chrX"),
#'     nKataegisFoci = 1,
#'     nKataegisVariants = 25,
#'     sampleName = "testSample",
#'     removeValidationColumns = FALSE
#' )
#'
#' @return (\link[VariantAnnotation]{VRanges}): VRanges containing background genomic variants and pre-defined kataegis foci.
#'
#' @author Daan Hazelaar
#' @author Job van Riet
#' @export
generateSyntheticData <- function(genome = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19,
                                  nBackgroundVariants = 100,
                                  seqnames = NULL,
                                  probMutationType = c(.8, .1, .1),
                                  nKataegisFoci = 5,
                                  nKataegisVariants = 20,
                                  expectedIMD = 100,
                                  sampleName = "syntheticData",
                                  removeValidationColumns = TRUE) {
    .validateInputGenerateSyntheticData(genome, nBackgroundVariants, seqnames, probMutationType, nKataegisFoci, nKataegisVariants, expectedIMD, sampleName, removeValidationColumns)

    selectedSequences <- .selectSequencesgenerateSyntheticData(seqnames, genome)

    dataInit <- .determineInitialProbabilities(selectedSequences, genome)

    kataegisData <- getKataegisData(nKataegisFoci, dataInit, nKataegisVariants, expectedIMD, genome)

    dataKataegisFoci <- getDataKataegisFoci(kataegisData)

    dataKataegisVariants <- getDataKataegisVariants(kataegisData)

    dataBackgroundVariants <- .getBackgroundData(dataInit, nBackgroundVariants, probMutationType, genome)

    dataCombined <- .combineData(dataKataegisVariants, dataBackgroundVariants, dataKataegisFoci, sampleName)

    dataFinal <- .removeValidationColumns(dataCombined, removeValidationColumns)

    return(dataFinal)
}

.validateInputGenerateSyntheticData <- function(genome, nBackgroundVariants, seqnames, probMutationType, nKataegisFoci, nKataegisVariants, expectedIMD, sampleName, removeValidationColumns) {
    checkmate::assertClass(genome, "BSgenome")
    checkmate::assertInt(nBackgroundVariants, lower = 0, upper = 1E6)
    checkmate::assertCharacter(seqnames, null.ok = TRUE, pattern = "^chr")
    checkmate::assertNumeric(probMutationType, lower = 0, upper = 10, len = 3, any.missing = FALSE)
    checkmate::assertInt(nKataegisFoci, lower = 0, upper = 1000)
    checkmate::assertInt(nKataegisVariants, lower = 0, upper = 1E6)
    checkmate::assertInt(expectedIMD, lower = 0, null.ok = FALSE)
    checkmate::assertCharacter(sampleName)
    checkmate::assertLogical(removeValidationColumns)

    # Check if all necessary parameters were provided
    if (nBackgroundVariants == 0 & any(c(nKataegisFoci, nKataegisVariants) == 0)) {
        base::stop("Make sure more than 0 (background or kataegis) mutations are generated")
    }
}

.selectSequencesgenerateSyntheticData <- function(seqnames, genome) {
    # If no sequence names are specified select the 22 autosomes and the sex chromosomes from the reference genome
    if (base::is.null(seqnames)) {
        seqnames <- base::names(genome)[base::seq_len(24)]
    }

    # check if all requested sequence names are present in the provided reference genome
    sequenceInGenome <- seqnames[!seqnames %in% base::names(genome)]
    if (!S4Vectors::isEmpty(sequenceInGenome)) {
        base::stop("The following sequence name(s) are not present in the provided reference genome: ", sequenceInGenome)
    }

    return(seqnames)
}

.determineInitialProbabilities <- function(selectedSequences, genome) {
    initData <- tibble::tibble(.rows = base::length(selectedSequences)) |>
        dplyr::mutate(
            # Name of sequence/chromosome.
            seqnames = selectedSequences,
            # Length of sequence/chromosome.
            length = GenomeInfoDb::seqlengths(genome)[.data$seqnames],
            # Probability of a mutation occurring in each seqname (assuming equal probability for all basepairs).
            probMut = .data$length / base::sum(.data$length)
        )
    return(initData)
}

getKataegisData <- function(nKataegisFoci, dataInit, nKataegisVariants, expectedIMD, genome) {
    kataegisData <- lapply(seq_len(nKataegisFoci), function(synthFociID) {
        generateKataegisData(synthFociID, dataInit, nKataegisVariants, expectedIMD, genome)
    })

    return(kataegisData)
}

generateKataegisData <- function(synthFociID, dataInit, nKataegisVariants, expectedIMD, genome) {
    dataKataegisFoci <- .generateKataegisFoci(synthFociID, dataInit, nKataegisFoci = 1, nKataegisVariants, expectedIMD)
    dataKataegisVariants <- .generateKataegisVariants(dataKataegisFoci, nKataegisVariants, genome)

    if (any(dataKataegisVariants$REF == "N")) {
        generateKataegisData(synthFociID, dataInit, nKataegisVariants, expectedIMD, genome)
    } else {
        kataegisData <- list(dataKataegisFoci = dataKataegisFoci, dataKataegisVariants = dataKataegisVariants)
        return(kataegisData)
    }
}

getDataKataegisFoci <- function(kataegisData) {
    dataKataegisFoci <- lapply(kataegisData, function(x) {
        return(x$dataKataegisFoci)
    }) |> dplyr::bind_rows()

    return(dataKataegisFoci)
}

getDataKataegisVariants <- function(kataegisData) {
    dataKataegisVariants <- lapply(kataegisData, function(x) {
        return(x$dataKataegisVariants)
    }) |> dplyr::bind_rows()

    return(dataKataegisVariants)
}

.generateKataegisFoci <- function(synthFociID, dataInit, nKataegisFoci, nKataegisVariants, expectedIMD) {
    if (nKataegisFoci != 0 & nKataegisVariants != 0) {
        # Retrieve the specified number of kataegis foci from random places on the genome.
        dataKataegisFoci <- tibble::tibble(
            # Determine the chromosomes on which kataegis foci are present.
            seqnames = base::sample(dataInit$seqnames, prob = dataInit$probMut, size = nKataegisFoci, replace = TRUE)
        ) |>
            # Add column with length of each the seqname.
            dplyr::inner_join(dplyr::select(dataInit, .data$seqnames, .data$length), by = "seqnames") |>
            dplyr::rowwise() |>
            dplyr::mutate(
                # Determine genomic start and end location of each kataegis region
                startKataegisFoci = base::sample(base::seq_len(.data$length - ((nKataegisVariants + 1) * expectedIMD)), size = 1, replace = TRUE),
                endKataegisFoci = .data$startKataegisFoci + ((nKataegisVariants + 1) * expectedIMD) - 1,
                synthFociID = synthFociID
            ) |>
            dplyr::ungroup()
    } else {
        dataKataegisFoci <- NULL
    }

    return(dataKataegisFoci)
}

# Generate kataegis variants
.generateKataegisVariants <- function(dataKataegisFoci, nKataegisVariants, genome) {
    if (!base::is.null(dataKataegisFoci)) {
        dataKataegisVariants <- dataKataegisFoci |>
            # Make a separate row for all mutations
            tidyr::uncount(nKataegisVariants) |>
            # Sample start and end point within each kataegis foci
            dplyr::group_by(.data$synthFociID) |>
            dplyr::mutate(
                start = base::sample(base::seq(base::unique(.data$startKataegisFoci), base::unique(.data$endKataegisFoci)), size = nKataegisVariants),
                end = .data$start,
                mutType = "SNV"
            ) |>
            dplyr::ungroup()

        # Retrieve the reference base from the reference genome.
        dataKataegisVariants$REF <- BSgenome::getSeq(
            genome,
            names = dataKataegisVariants$seqnames,
            start = dataKataegisVariants$start,
            end = dataKataegisVariants$end,
            as.character = TRUE
        ) |>
            base::unname()

        # The alternative base for a SNV is different base than REF
        dataKataegisVariants$ALT <- base::vapply(
            dataKataegisVariants$REF,
            FUN = function(x) {
                nucleotides <- c("A", "T", "C", "G")
                base::sample(nucleotides[nucleotides != x], size = 1)
            }, FUN.VALUE = "s"
        ) |>
            base::unname()
    } else {
        dataKataegisVariants <- NULL
    }

    return(dataKataegisVariants)
}

.getBackgroundData <- function(dataInit, nBackgroundVariants, probMutationType, genome) {
    if (nBackgroundVariants != 0) {
        backgroundData <- .generateBackgroundVariants(dataInit, nBackgroundVariants, probMutationType, genome)
    } else {
        backgroundData <- NULL
    }

    return(backgroundData)
}

.generateBackgroundVariants <- function(dataInit, nBackgroundVariants, probMutationType, genome, previous = NULL) {
    # Construct a tibble which contains a row for each mutation event
    dataBackgroundVariants <- tibble::tibble(
        # Determine on which seqname the mutations occur according to previously determined probabilities
        seqnames = base::sample(dataInit$seqnames, prob = dataInit$probMut, size = nBackgroundVariants, replace = TRUE),
        # Determine the type of mutation (SNV, Deletion or Insertion) for each mutation event.
        mutType = base::sample(c("SNV", "Deletion", "Insertion"), size = nBackgroundVariants, prob = probMutationType, replace = TRUE)
    ) |>
        # Add column with length of each the seqname.
        dplyr::inner_join(dplyr::select(dataInit, .data$seqnames, .data$length), by = "seqnames") |>
        # Determine the number a basepairs that are involved in each mutation event.
        dplyr::mutate(
            mutLength = dplyr::case_when(
                .data$mutType == "SNV" ~ 1,
                .data$mutType == "Deletion" ~ base::as.numeric(sample(2:10, size = nBackgroundVariants, replace = TRUE)),
                .data$mutType == "Insertion" ~ base::as.numeric(sample(2:9, size = nBackgroundVariants, replace = TRUE))
            )
        ) |>
        # Perform per row.
        dplyr::rowwise() |>
        # determine genomic start location for each mutation event
        dplyr::mutate(
            start = base::sample(base::seq_len(.data$length), size = 1, replace = TRUE)
        ) |>
        dplyr::group_by(.data$mutType) |>
        # determine genomic end location for each mutation event
        dplyr::mutate(
            end = dplyr::if_else(.data$mutType == "SNV" | .data$mutType == "Insertion", base::as.numeric(.data$start), base::as.numeric(.data$start + .data$mutLength - 1))
        ) |>
        dplyr::ungroup()

    # Retrieve the reference base(s) from the reference genome.
    dataBackgroundVariants$REF <- BSgenome::getSeq(
        genome,
        names = dataBackgroundVariants$seqnames,
        start = dataBackgroundVariants$start,
        end = dataBackgroundVariants$end,
        as.character = TRUE
    ) |>
        base::unname()

    # Determine alternate base(s) for the different mutation types
    # The alternative base for a SNV is different base than REF
    dataSNV <- dataBackgroundVariants[dataBackgroundVariants$mutType == "SNV", ]
    dataSNV$ALT <- base::vapply(dataSNV$REF, FUN = function(x) {
        nucleotides <- c("A", "T", "C", "G")
        base::sample(nucleotides[nucleotides != x], size = 1)
    }, FUN.VALUE = "s") |>
        base::unname()

    # The alternative base for an insertion is the first base of REF
    dataInsertion <- dataBackgroundVariants[dataBackgroundVariants$mutType == "Insertion", ]
    dataInsertion$ALT <- base::paste0(dataInsertion$REF, base::vapply(
        dataInsertion$mutLength,
        FUN = function(x) {
            nucleotides <- c("A", "T", "C", "G")
            base::paste(base::sample(nucleotides, size = x, replace = TRUE), collapse = "")
        }, FUN.VALUE = "s"
    )) |>
        base::unname()

    # The alternative bases for a deletion start with the reference base followed by a sequence of random bases according to mutation length
    dataDeletion <- dataBackgroundVariants[dataBackgroundVariants$mutType == "Deletion", ]
    dataDeletion$ALT <- base::substring(dataDeletion$REF, 1, 1)

    # combine data of different mutation types
    dataBackgroundVariants <- dplyr::bind_rows(dataSNV, dataInsertion, dataDeletion)

    # remove variants that contain an N base
    dataBackgroundVariantsNoN <- dataBackgroundVariants |> dplyr::filter(!base::grepl("N", .data$REF))
    nResample <- nBackgroundVariants - base::nrow(dataBackgroundVariantsNoN)

    # combine all sampled variants in single tibble
    dataBackgroundVariantsCombined <- dplyr::bind_rows(dataBackgroundVariantsNoN, previous)

    # recursion is used for resampling mutations that occured in N regions
    if (nResample != 0) {
        .generateBackgroundVariants(dataInit = dataInit, nBackgroundVariants = nResample, probMutationType = probMutationType, genome = genome, previous = dataBackgroundVariantsCombined)
    } else {
        return(dataBackgroundVariantsCombined)
    }
}

.combineData <- function(dataKataegisVariants, dataBackgroundVariants, dataKataegisFoci, sampleName) {
    # Combine background and kataegis mutations.
    dataCombined <- dplyr::bind_rows(dataKataegisVariants, dataBackgroundVariants) |>
        dplyr::mutate(sampleNames = sampleName) |>
        GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE) |>
        VariantAnnotation::makeVRangesFromGRanges() |>
        GenomicRanges::sort()

    # add column which specifies if a variants lies within a kataegis foci
    if (!nrow(dataKataegisFoci) == 0) {
        # convert to GRanges
        colnames(dataKataegisFoci)[c(3, 4)] <- c("start", "end")
        dataKataegisFocigr <- GenomicRanges::makeGRangesFromDataFrame(dataKataegisFoci)

        # find all variants that lie within a kataegis foci
        kataegisVariants <- S4Vectors::queryHits(IRanges::findOverlaps(dataCombined, dataKataegisFocigr))

        # add column to the combined data
        dataCombined$kataegis <- FALSE
        dataCombined[kataegisVariants]$kataegis <- TRUE
    } else {
        dataCombined$kataegis <- FALSE
    }

    return(dataCombined)
}

.removeValidationColumns <- function(dataCombined, removeValidationColumns) {
    # remove validation columns if so specified
    if (removeValidationColumns) {
        dataCombined <- dataCombined |>
            tibble::as_tibble() |>
            dplyr::select(c(.data$seqnames, .data$start, .data$end, .data$width, .data$strand, .data$ref, .data$alt, .data$totalDepth, .data$refDepth, .data$altDepth, .data$sampleNames)) |>
            GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE) |>
            VariantAnnotation::makeVRangesFromGRanges()
    } else {
        dataCombined <- dataCombined
    }

    return(dataCombined)
}
