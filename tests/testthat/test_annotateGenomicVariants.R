
testthat::test_that("test .annotateGenomicVariants:", {

    genomicVariantsAnnotated <- system.file('extdata', 'CPTAC_Breast.vcf', package = 'katdetectr') |>
        .importGenomicVariants() |>
        .processGenomicVariants() |>
        .annotateGenomicVariants()

    testthat::expect_equal(genomicVariantsAnnotated$IMD[1], base::as.numeric(NA))
    testthat::expect_equal(genomicVariantsAnnotated$IMD[2], 14385)
})

testthat::test_that("test .determineIMD():", {

    vr <- VariantAnnotation::VRanges(
        seqnames = "test",
        ranges = IRanges::IRanges(start = c(1, 6, 20, 25, 31), end = c(5, 6, 20, 27, 33)),
        ref = c("A", "A", "T", "A", "T"),
        alt = c("T", "C", "A", "C", "A")
        )

    vrIMD <- .determineIMD(vr)

    testthat::expect_equal(vrIMD$IMD[1], base::as.numeric(NA))
    testthat::expect_equal(vrIMD$IMD[2], 0)
    testthat::expect_equal(vrIMD$IMD[3], 13)
    testthat::expect_equal(vrIMD$IMD[4], 4)
    testthat::expect_equal(vrIMD$IMD[5], 3)
})
