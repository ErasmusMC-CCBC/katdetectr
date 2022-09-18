testthat::test_that("test .annotateSegments():", {

    # test on maf file
    genomicVariantsAnnotatedCPTAC <- system.file('extdata', 'CPTAC_Breast.vcf', package = 'katdetectr') |>
        .importGenomicVariants() |>
        .processGenomicVariants() |>
        .annotateGenomicVariants()

    changepointsCPTAC <- .performChangepointDetection(genomicVariantsAnnotated = genomicVariantsAnnotatedCPTAC, test.stat = "Exponential", penalty = "BIC", pen.value = 0, minseglen = 2, BPPARAM = BiocParallel::SerialParam())
    segmentsCPTAC <- .annotateSegments(changepoints = changepointsCPTAC, genomicVariantsAnnotated = genomicVariantsAnnotatedCPTAC)

    testthat::expect_equal(base::length(segmentsCPTAC), 454)
    testthat::expect_equal(base::unique(segmentsCPTAC$sampleNames), "CPTAC")
    testthat::expect_equal(base::as.character(base::unique(GenomicRanges::seqnames(segmentsCPTAC)[1])), "chr1")
    testthat::expect_equal(segmentsCPTAC$segmentID[5], 5)
    testthat::expect_equal(GenomicRanges::end(GenomicRanges::ranges(segmentsCPTAC)[446]), 50728062)
    testthat::expect_equal(segmentsCPTAC$lastVariantID[446], 3622)
    testthat::expect_equal(segmentsCPTAC$totalVariants[443], 18)
    testthat::expect_equal(segmentsCPTAC$meanIMD[4], 1081.25)
    testthat::expect_equal(segmentsCPTAC$mutationRate[4], 0.0009248555)
    testthat::expect_equal(segmentsCPTAC@ranges@start[1], 1)
    testthat::expect_equal(segmentsCPTAC$totalVariants[1] / segmentsCPTAC@ranges@width[1], segmentsCPTAC$mutationRate[1])
    testthat::expect_equal(segmentsCPTAC@ranges@start[449], 1)
    testthat::expect_equal(segmentsCPTAC@ranges@start[450], 3228145)
    testthat::expect_equal(segmentsCPTAC$totalVariants[454] / segmentsCPTAC@ranges@width[454], segmentsCPTAC$mutationRate[454])

    # The weighted mean of the mutation rate of all segments in a chromosome must equal the rate of the entire chromosome
    testTotalRateChrx <- segmentsCPTAC |>
        dplyr::as_tibble() |>
        dplyr::filter(seqnames == "chrX")

    MutationRateChrX <- base::sum(testTotalRateChrx$totalVariants) / base::sum(testTotalRateChrx$width)
    meanRateSegments <- stats::weighted.mean(testTotalRateChrx$mutationRate, testTotalRateChrx$width)

    testthat::expect_equal(MutationRateChrX, meanRateSegments)
})


testthat::test_that("test .determineSegmentID:", {

    changepoints <- list(testSample1 = c(0 , 5, 43, 46, 52, 56), testSample2 = c(0, 5))

    segmentIDs <- .determineSegmentID(changepoints)

    testthat::expect_equal(base::length(segmentIDs), 61)
    testthat::expect_equal(segmentIDs[1], 1)
    testthat::expect_equal(segmentIDs[56], 5)
    testthat::expect_equal(segmentIDs[61], 1)
})

testthat::test_that("test .determineSegmentIDperChr:", {

    changepoints <- c(0 , 5, 43, 46, 52, 56, 100, 10200)

    segmentIDs <- .determineSegmentIDperChr(changepoints)

    testthat::expect_equal(base::length(segmentIDs), 10200)
    testthat::expect_equal(segmentIDs[1], 1)
    testthat::expect_equal(segmentIDs[56], 5)
    testthat::expect_equal(segmentIDs[61], 6)
})

