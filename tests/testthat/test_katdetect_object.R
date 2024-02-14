testthat::test_that("test KatDetect object:", {
    kd <- detectKataegis(genomicVariants = system.file("extdata", "CPTAC_Breast.vcf", package = "katdetectr"))

    testthat::expect_equal(katdetectr::getGenomicVariants(kd)$IMD[1], 935222)
    testthat::expect_equal(katdetectr::getKataegisFoci(kd)$meanIMD[1], 435.16667)
    testthat::expect_equal(katdetectr::getSegments(kd)$firstVariantID[1], 1)
    testthat::expect_equal(katdetectr::getInfo(kd)$sampleName, "CPTAC")
})
