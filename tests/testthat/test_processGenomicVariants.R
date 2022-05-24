
testthat::test_that("test .reduceOverlappingVariants()", {

    vrReduced <- system.file('extdata', 'CPTAC_Breast.vcf', package = 'katdetectr') |>
        .importGenomicVariants() |>
        .processGenomicVariants()

    testthat::expect_equal(base::length(vrReduced), 3684)
    testthat::expect_equal(vrReduced@ref[33], "X")
    testthat::expect_equal(vrReduced@alt[1254], "XX")
    testthat::expect_equal(vrReduced$revmap[[33]], c(33, 34))
    testthat::expect_equal(vrReduced$variantID[33], 33)
})

