testthat::test_that("test .performChangepointDetection:", {

    # test on larger breast cancer sample
    genomicVariantsAnnotatedCPTAC <- system.file('extdata', 'CPTAC_Breast.vcf', package = 'katdetectr') |>
        .importGenomicVariants() |>
        .processGenomicVariants() |>
        .annotateGenomicVariants()

    resCpCPTAC <- .performChangepointDetection(genomicVariantsAnnotated = genomicVariantsAnnotatedCPTAC, test.stat = "Exponential", penalty = "BIC", pen.value = 0, method = "PELT", minseglen = 2, BPPARAM = BiocParallel::SerialParam())

    testthat::expect_equal(base::length(resCpCPTAC), 23)
    testthat::expect_equal(resCpCPTAC$chr1$changepointsChromosome[1], 0)
    testthat::expect_equal(resCpCPTAC$chr1$changepointsChromosome[2], 11)
    testthat::expect_equal(resCpCPTAC$chr1$changepointsChromosome[48], 435)

    testthat::expect_equal(round(resCpCPTAC$chr1$rateChromosome[1], 4), 0)
    testthat::expect_equal(round(resCpCPTAC$chr1$rateChromosome[2], 4), 0.0001)
    testthat::expect_equal(round(resCpCPTAC$chr1$rateChromosome[27], 4), 0.0150)
    testthat::expect_equal(round(resCpCPTAC$chr1$rateChromosome[47], 4), 0)

    testthat::expect_equal(resCpCPTAC$chrX$changepointsChromosome[1], 0)
    testthat::expect_equal(resCpCPTAC$chrX$changepointsChromosome[2], 2)
    testthat::expect_equal(resCpCPTAC$chrX$changepointsChromosome[7], 56)
})
