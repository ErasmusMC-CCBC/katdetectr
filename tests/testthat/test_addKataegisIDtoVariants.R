testthat::test_that("test .addIDsToVariants()", {

    # test on larger breast cancer sample
    genomicVariantsAnnotatedCPTAC <- system.file('extdata', 'CPTAC_Breast.vcf', package = 'katdetectr') |>
        .importGenomicVariants() |>
        .processGenomicVariants() |>
        .annotateGenomicVariants()
    changepointsCPTAC <- .performChangepointDetection(genomicVariantsAnnotated = genomicVariantsAnnotatedCPTAC, test.stat = "Exponential", penalty = "BIC", pen.value = 0, method = "PELT", minseglen = 2, BPPARAM = BiocParallel::SerialParam())
    segmentsCPTAC <- .annotateSegments(changepoints = changepointsCPTAC, genomicVariantsAnnotated = genomicVariantsAnnotatedCPTAC)
    kataegisFociCPTAC <- .determineKataegisFoci(segments = segmentsCPTAC, genomicVariantsAnnotated = genomicVariantsAnnotatedCPTAC, minSizeKataegis = 5, IMDcutoff = 1000)
    genomicVariantsAnnotatedCPTACkat <- .addIDsToVariants(kataegisFoci = kataegisFociCPTAC, genomicVariantsAnnotated = genomicVariantsAnnotatedCPTAC, changepointsPerChromosome = changepointsCPTAC)

    testthat::expect_equal(base::length(genomicVariantsAnnotatedCPTACkat), 3684)
    testthat::expect_equal(genomicVariantsAnnotatedCPTACkat$putativeKataegis[1], FALSE)
    testthat::expect_equal(genomicVariantsAnnotatedCPTACkat$putativeKataegis[363], TRUE)
    testthat::expect_equal(genomicVariantsAnnotatedCPTACkat$putativeKataegis[367], TRUE)
    testthat::expect_equal(genomicVariantsAnnotatedCPTACkat$putativeKataegis[368], FALSE)
    testthat::expect_equal(genomicVariantsAnnotatedCPTACkat$putativeKataegis[3138], FALSE)
    testthat::expect_equal(genomicVariantsAnnotatedCPTACkat$putativeKataegis[3139], TRUE)
    testthat::expect_equal(genomicVariantsAnnotatedCPTACkat$putativeKataegis[3157], TRUE)
    testthat::expect_equal(genomicVariantsAnnotatedCPTACkat$putativeKataegis[3158], FALSE)


    testthat::expect_equal(genomicVariantsAnnotatedCPTACkat$segmentID[1], 1)
    testthat::expect_equal(genomicVariantsAnnotatedCPTACkat$segmentID[3634], 3)
    testthat::expect_equal(genomicVariantsAnnotatedCPTACkat$segmentID[3671], 3)
    testthat::expect_equal(genomicVariantsAnnotatedCPTACkat$segmentID[3684], 7)

})
