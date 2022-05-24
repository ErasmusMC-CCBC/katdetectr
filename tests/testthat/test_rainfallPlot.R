
testthat::test_that("test rainfallPlot():", {

    kd <- detectKataegis(genomicVariants = system.file('extdata', 'CPTAC_Breast.vcf', package = 'katdetectr'))

    testthat::expect_error(rainfallPlot(kd = kd, showSequence = c("chr1", "ch")), regexp = "Requested sequences are not present in the data.")
    testthat::expect_error(rainfallPlot(kd = kd, showSequence = "chr300"), regexp = "Requested sequences are not present in the data.")

    p1 <- rainfallPlot(kd = kd, showSequence = "All")
    testthat::expect_equal(base::nrow(p1$data), 3684)

    p2 <- rainfallPlot(kd = kd, showSequence = "Kataegis")
    testthat::expect_equal(base::nrow(p2$data), 1342)

    p3 <- rainfallPlot(kd = kd, showSequence = c("chr1", "chr2"))
    testthat::expect_equal(base::nrow(p3$data), 682)
})


