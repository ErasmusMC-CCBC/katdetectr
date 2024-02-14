testthat::test_that("test non standard sequences:", {
    set.seed(1)

    syndata1 <- generateSyntheticData(seqnames = c("chr1_gl000191_random", "chr4_ctg9_hap1"))
    syndata2 <- generateSyntheticData(seqnames = "chr1", nKataegisFoci = 0)
    syndata <- suppressWarnings(c(syndata1, syndata2))

    set.seed(NULL)

    sequenceLength1 <- data.frame(
        chr1 = 249250621,
        chr1_gl000191_random = 106433,
        chr4_ctg9_hap1 = 590426
    )

    kd <- detectKataegis(genomicVariants = syndata, refSeq = sequenceLength1)

    sequenceLength2 <- data.frame(
        chr1 = 249250621,
        chr1_gl000191_random = 106433
    )

    testthat::expect_error(detectKataegis(genomicVariants = syndata, refSeq = sequenceLength2))
})

testthat::test_that("test reference genome:", {
    testthat::expect_error(detectKataegis(genomicVariants = system.file("extdata", "CPTAC_Breast.vcf", package = "katdetectr"), refSeq = "hg38"))
})
