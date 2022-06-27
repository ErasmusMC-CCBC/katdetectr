

testthat::test_that("test .performChangepointDetection:", {

    # test on larger breast cancer sample
    genomicVariantsAnnotatedCPTAC <- system.file('extdata', 'CPTAC_Breast.vcf', package = 'katdetectr') |>
        .importGenomicVariants() |>
        .processGenomicVariants() |>
        .annotateGenomicVariants()

    resPCFCPTAC <- .performChangepointDetection(genomicVariantsAnnotated = genomicVariantsAnnotatedCPTAC, test.stat = "Exponential", penalty = "BIC", pen.value = 0, minseglen = 2, BPPARAM = BiocParallel::SerialParam())



    testthat::expect_equal(base::length(resPCFCPTAC), 23)
    testthat::expect_equal(resPCFCPTAC$chr1$changepointsChromosome[1], 0)
    testthat::expect_equal(resPCFCPTAC$chr1$changepointsChromosome[2], 11)
    testthat::expect_equal(resPCFCPTAC$chr1$changepointsChromosome[48], 435)

    testthat::expect_equal(round(resPCFCPTAC$chr1$rateChromosome[1], 4), 0)
    testthat::expect_equal(round(resPCFCPTAC$chr1$rateChromosome[2], 4), 0.0001)
    testthat::expect_equal(round(resPCFCPTAC$chr1$rateChromosome[27], 4), 0.0153)
    testthat::expect_equal(round(resPCFCPTAC$chr1$rateChromosome[47], 4), 0)

    testthat::expect_equal(resPCFCPTAC$chrX$changepointsChromosome[1], 0)
    testthat::expect_equal(resPCFCPTAC$chrX$changepointsChromosome[2], 5)
    testthat::expect_equal(resPCFCPTAC$chrX$changepointsChromosome[6], 56)
})
