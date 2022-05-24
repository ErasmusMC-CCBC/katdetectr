testthat::test_that("test .importGenomicVariants():", {

    testthat::expect_error(.importGenomicVariants(system.file('extdata', 'APL_primary.maf', package = 'katdetectr')))

    vrCPTAC <- .importGenomicVariants(system.file('extdata', 'CPTAC_Breast.vcf', package = 'katdetectr'))

    testthat::expect_s4_class(vrCPTAC, "VRanges")
    testthat::expect_equal(base::length(vrCPTAC), 3687)
    testthat::expect_equal(vrCPTAC@ref[1], "C")
})

testthat::test_that("test coerceMAFtoVRanges():", {

    vr <- .coerceMAFtoVRanges(path = system.file('extdata', 'APL_primary.maf', package = 'katdetectr'))

    testthat::expect_s4_class(vr, "VRanges")
    testthat::expect_equal(base::length(vr), 224)
    testthat::expect_equal(base::length(base::levels(base::unique(vr@sampleNames))), 97)
    testthat::expect_equal(vr@ref[2], "G")
    testthat::expect_equal(base::as.character(vr@seqnames[1]), "chr17")
})

testthat::test_that("test coerceVCFtoVRanges():", {

    pathToVCF <- system.file('extdata', 'CPTAC_Breast.vcf', package = 'katdetectr')
    vr <- .coerceVCFtoVRanges(path = pathToVCF)

    testthat::expect_s4_class(vr, "VRanges")
    testthat::expect_equal(base::length(vr), 3687)
    testthat::expect_equal(base::length(base::levels(base::unique(vr@sampleNames))), 1)
    testthat::expect_equal(vr@ref[21], "G")
    testthat::expect_equal(base::as.character(vr@seqnames[1]), "chr1")
})
