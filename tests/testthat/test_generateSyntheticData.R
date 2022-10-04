testthat::test_that("test .generateBackgroundVariants():", {

    base::set.seed(1)

    syntheticDataTest <- generateSyntheticData(
        genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
        nBackgroundVariants = 75,
        seqnames = c("chr1", "chrX"),
        nKataegisFoci = 1,
        nKataegisVariants = 25,
        sampleName = "testSample",
        removeValidationColumns = FALSE
    )

    syntheticDataTest2 <- generateSyntheticData(
        genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
        nBackgroundVariants = 100,
        seqnames = c("chr1", "chrX"),
        nKataegisFoci = 0,
        nKataegisVariants = 25,
        sampleName = "testSample",
        removeValidationColumns = FALSE
    )

    base::set.seed(NULL)


    testthat::expect_equal(base::length(syntheticDataTest), 100)
    testthat::expect_equal(base::levels(VariantAnnotation::sampleNames(syntheticDataTest)[1]), "testSample")
    testthat::expect_equal(VariantAnnotation::ref(syntheticDataTest)[3], "G")
    testthat::expect_equal(VariantAnnotation::ref(syntheticDataTest)[6], "C")
    testthat::expect_equal(VariantAnnotation::alt(syntheticDataTest)[9], "T")
    testthat::expect_equal(syntheticDataTest$kataegis[1], TRUE)
    testthat::expect_equal(syntheticDataTest$kataegis[25], TRUE)
    testthat::expect_equal(syntheticDataTest$kataegis[26], FALSE)
    testthat::expect_equal(VariantAnnotation::ref(syntheticDataTest)[30], "G")
    testthat::expect_equal(VariantAnnotation::ref(syntheticDataTest)[31], "A")
    testthat::expect_equal(VariantAnnotation::alt(syntheticDataTest)[31], "G")
    testthat::expect_equal(sum(grepl("N", VariantAnnotation::ref(syntheticDataTest))), 0)
    testthat::expect_equal(sum(grepl("N", VariantAnnotation::alt(syntheticDataTest))), 0)

    testthat::expect_equal(base::length(syntheticDataTest2), 100)
    testthat::expect_equal(sum(syntheticDataTest2$kataegis), 0)
})

testthat::test_that("test .generateBackgroundVariants():", {

    base::set.seed(1)

    hg38 <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
    testDataInit <- tibble::tibble(
        seqnames= c("chr1", "chr2"),
        length = c(248956422, 242193529),
        probMut = c(.8, .2)
    )

    testVariants <- .generateBackgroundVariants(dataInit = testDataInit, nBackgroundVariants = 100, probMutationType = c(.1, .1, .1), genome = hg38)

    base::set.seed(NULL)

    testthat::expect_equal(base::nrow(testVariants), 100)
    testthat::expect_equal(testVariants$seqnames[5], "chr2")
    testthat::expect_equal(testVariants$mutType[5], "SNV")
    testthat::expect_equal(testVariants$REF[5], "A")
    testthat::expect_equal(testVariants$ALT[5], "C")
    testthat::expect_equal(testVariants$mutType[50], "Insertion")
    testthat::expect_equal(testVariants$REF[50], "A")
    testthat::expect_equal(testVariants$ALT[50], "AGCA")
    testthat::expect_equal(testVariants$mutType[90], "Deletion")
    testthat::expect_equal(testVariants$REF[90], "GG")
    testthat::expect_equal(testVariants$ALT[90], "G")


})


testthat::test_that("test .generateKataegisVariants():", {

    base::set.seed(1)

    hg38 <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
    dataKataegisFociTest <- tibble::tibble(
        synthFociID = c(1, 2),
        seqnames = c("chr1", "chr1"),
        startKataegisFoci = c(1000000, 1),
        endKataegisFoci = c(1001000, 100)
    )
    testKataegisVariants <- .generateKataegisVariants(dataKataegisFoci = dataKataegisFociTest, nKataegisVariants = 15, genome = hg38)

    base::set.seed(NULL)

    testthat::expect_equal(testKataegisVariants$REF[1], "C")
    testthat::expect_equal(testKataegisVariants$ALT[1], "T")
    testthat::expect_equal(testKataegisVariants$start[5], 1000508)
    testthat::expect_equal(testKataegisVariants$REF[30], "N")
    testthat::expect_equal(testKataegisVariants$ALT[30], "C")
})

testthat::test_that("test .generateKataegisFoci():", {

    base::set.seed(1)

    hg38 <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
    initDataTest1 <- .determineInitialProbabilities(selectedSequences = c("chr1", "chrX"), genome = hg38)
    testFoci <- .generateKataegisFoci(synthFociID = 1 ,dataInit = initDataTest1, nKataegisFoci = 5, nKataegisVariants = 10, expectedIMD = 500)

    base::set.seed(NULL)

    testthat::expect_equal(base::nrow(testFoci), 5)
    testthat::expect_equal(testFoci$synthFociID[1], 1)
    testthat::expect_equal(testFoci$seqnames[1], "chr1")
    testthat::expect_equal(testFoci$startKataegisFoci[3], 221425633)
    testthat::expect_equal(testFoci$endKataegisFoci[5], 233666245)
})

testthat::test_that("test .selectSequencesgenerateSyntheticData():", {

    hg38 <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38

    selectedSequences1 <- .selectSequencesgenerateSyntheticData(seqnames = c("chr1", "chrX"), genome = hg38)
    testthat::expect_equal(selectedSequences1, c("chr1", "chrX"))

    selectedSequences2 <- .selectSequencesgenerateSyntheticData(seqnames = NULL, genome = hg38)
    testthat::expect_equal(base::length(selectedSequences2), 24)
    testthat::expect_equal(selectedSequences2[5], "chr5")

    selectedSequences3 <- .selectSequencesgenerateSyntheticData(seqnames = NULL, genome = hg38)
    testthat::expect_equal(base::length(selectedSequences3), 24)

    testthat::expect_error(.selectSequencesgenerateSyntheticData(seqnames = "test", genome = hg38))
})

testthat::test_that("test .determineInitProbabilities():", {

    hg38 <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38

    initDataTest1 <- .determineInitialProbabilities(selectedSequences = c("chr1", "chrX"), genome = hg38)
    testthat::expect_equal(initDataTest1$seqnames, c("chr1", "chrX"))
    testthat::expect_equal(base::sum(initDataTest1$probMut), 1)

    initDataTest2 <- .determineInitialProbabilities(selectedSequences = c("chr1", "chrX", "chr2", "chr18"), genome = hg38)
    testthat::expect_equal(initDataTest2$seqnames, c("chr1", "chrX", "chr2", "chr18"))
    testthat::expect_equal(base::sum(initDataTest2$probMut), 1)
})

