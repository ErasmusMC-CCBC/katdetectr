
testthat::test_that("test KatDetect object:", {

    kd <- detectKataegis(genomicVariants = system.file('extdata', 'CPTAC_Breast.vcf', package = 'katdetectr'))

    testthat::expect_equal(kd@genomicVariants$IMD[1], base::as.numeric(NA))
    testthat::expect_equal(kd@kataegisFoci$meanIMD[1], 434.16667)
    testthat::expect_equal(kd@segments$firstVariantID[1], 1)
    testthat::expect_equal(kd@info$sampleName, "CPTAC")
})

