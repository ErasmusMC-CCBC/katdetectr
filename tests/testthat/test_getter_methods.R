
testthat::test_that("Check accessor methods:", {

    set.seed(1)

    syntheticData <- generateSyntheticData(nBackgroundVariants = 20, nKataegisFoci = 1)
    kd <- detectKataegis(genomicVariants = syntheticData)

    testthat::expect_equal(base::length(getGenomicVariants(kd)), 40)
    testthat::expect_equal(base::length(getKataegisFoci(kd)), 1)
    testthat::expect_equal(getSegments(kd)$firstVariantID[1], 1)
    testthat::expect_equal(getInfo(kd)$sampleName, "syntheticData")
    testthat::expect_equal(getInfo(kd)$parameters$IMDcutoff, 1000)

    set.seed(NULL)
})
