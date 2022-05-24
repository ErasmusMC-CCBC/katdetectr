

testthat::test_that("test .performChangepointDetection:", {

    # test on larger breast cancer sample
    genomicVariantsAnnotatedCPTAC <- system.file('extdata', 'CPTAC_Breast.vcf', package = 'katdetectr') |>
        .importGenomicVariants() |>
        .processGenomicVariants() |>
        .annotateGenomicVariants()

    resPCFCPTAC <- .performChangepointDetection(genomicVariantsAnnotated = genomicVariantsAnnotatedCPTAC, test.stat = "Exponential", penalty = "BIC", pen.value = 0, minseglen = 2, bpworkers = 1)

    testthat::expect_equal(base::length(resPCFCPTAC), 23)
    testthat::expect_equal(resPCFCPTAC$chr1[1], 0)
    testthat::expect_equal(resPCFCPTAC$chr1[2], 11)
    testthat::expect_equal(resPCFCPTAC$chr1[48], 435)

    testthat::expect_equal(resPCFCPTAC$chrX[1], 0)
    testthat::expect_equal(resPCFCPTAC$chrX[2], 5)
    testthat::expect_equal(resPCFCPTAC$chrX[6], 56)
})
#
#
# testthat::test_that("test .splitPerChromosome:", {
#
#     genomicVariantsTest <- VariantAnnotation::VRanges(
#         seqnames = c("chr1", "chr2", "chr2", "chr1", "chr3"),
#         ranges = IRanges::IRanges(start = c(1, 60, 20, 25, 31), end = c(5, 60, 20, 27, 33)),
#         ref = c("A", "A", "T", "A", "T"),
#         alt = c("T", "C", "A", "C", "A")
#     )
#     vrSplit <- .splitPerChromosome(genomicVariantsTest)
#
#     testthat::expect_equal(base::length(vrSplit), 3)
#     testthat::expect_equal(base::length(vrSplit$chr1), 2)
# })


# testthat::test_that("test .getIMD:", {
#
#     testData <- base::list(
#         chr1 = tibble::tibble(IMD = c(0, 1, 5, 2, 8, 19, 0)),
#         chr2 = tibble::tibble(IMD = c(550, 100019, 10)))
#
#     IMD <- .getIMD(x = testData)
#
#     testthat::expect_equal(IMD$chr1[1], 0.001)
#     testthat::expect_equal(IMD$chr1[2], 1)
#     testthat::expect_equal(IMD$chr1[7], 0.001)
#     testthat::expect_equal(IMD$chr2[2], 100019)
# })



# testthat::test_that("test .runPCF:", {
#
#     pcfBIC <- .runCD(IMDsample = base::list(IMD = c(1,1,4,2,6,65,76,87,98)),
#                       test.stat = "Exponential",
#                       penalty = "BIC",
#                       pen.value = 0,
#                       minseglen = 2,
#                       bpworkers = 1)
#
#     pcfManual <- .runCD(IMDsample = base::list(IMD = c(1,1,4,2,6,65,76,87,98)),
#                          test.stat = "Exponential",
#                          penalty = "Manual",
#                          pen.value = 100,
#                          minseglen = 2,
#                          bpworkers = 1)
#
#     pcfEmpirical <- .runCD(IMDsample = base::list(IMD = c(1,1,4,2,6,65,76,87,98)),
#                             test.stat = "Empirical",
#                             penalty = "BIC",
#                             pen.value = 0,
#                             minseglen = 3,
#                             bpworkers = 1)
#
#     testthat::expect_equal(pcfBIC$IMD, c(0, 6, 10))
#     testthat::expect_equal(pcfManual$IMD, c(0, 10))
#     testthat::expect_equal(pcfEmpirical$IMD, c(0, 5, 10))
# })
