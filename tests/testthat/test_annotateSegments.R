testthat::test_that("test .annotateSegments():", {

    # test on maf file
    genomicVariantsAnnotatedCPTAC <- system.file('extdata', 'CPTAC_Breast.vcf', package = 'katdetectr') |>
    .importGenomicVariants() |>
    .processGenomicVariants() |>
    .annotateGenomicVariants()
    changepointsCPTAC <- .performChangepointDetection(genomicVariantsAnnotated = genomicVariantsAnnotatedCPTAC, test.stat = "Exponential", penalty = "BIC", pen.value = 0, minseglen = 2, bpworkers = 1)
    segmentsCPTAC <- .annotateSegments(changepoints = changepointsCPTAC, genomicVariantsAnnotated = genomicVariantsAnnotatedCPTAC)

    testthat::expect_equal(base::length(segmentsCPTAC), 446)
    testthat::expect_equal(base::unique(segmentsCPTAC$sampleNames), "CPTAC")
    testthat::expect_equal(base::as.character(base::unique(segmentsCPTAC@seqnames[1])), "chr1")
    testthat::expect_equal(segmentsCPTAC$segmentID[5], 5)
    testthat::expect_equal(GenomicRanges::end(segmentsCPTAC@ranges[446]), 153764217)
    testthat::expect_equal(segmentsCPTAC$lastVariantID[446], 3684)
    testthat::expect_equal(segmentsCPTAC$lastVariantID[446], 3684)
    testthat::expect_equal(segmentsCPTAC$totalVariants[443], 38)
    testthat::expect_equal(segmentsCPTAC$meanIMD[4], 1080.25)
    testthat::expect_equal(segmentsCPTAC$mutationRate[4], 0.00092571164)
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
