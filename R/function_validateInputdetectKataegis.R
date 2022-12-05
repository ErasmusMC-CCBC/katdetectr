validateInputdetectKataegis <- function(genomicVariants, refSeq, minSizeKataegis, IMDcutoff, test.stat, penalty, pen.value, method, minseglen, BPPARAM, aggregateRecords){

    checkmate::assert(
        checkmate::checkClass(genomicVariants, 'VRanges'),
        checkmate::checkAccess(genomicVariants, access = 'r')
    )

    checkmate::assert(
        checkmate::checkClass(refSeq, "data.frame"),
        checkmate::checkCharacter(refSeq, pattern = "hg19|hg38")
    )

    checkmate::assertInt(minSizeKataegis, lower = 1)

    checkmate::assert(
        checkmate::checkNumeric(IMDcutoff, 0),
        checkmate::checkClass(IMDcutoff, "function")
    )

    checkmate::assert_character(test.stat, pattern = 'Exponential')
    checkmate::assert_character(penalty, pattern = 'BIC|Manual')
    checkmate::assertInt(pen.value, lower = 0)
    checkmate::assert_character(method, pattern = 'PELT|AMOC|SegNeigh|BinSeg')
    checkmate::assertInt(minseglen, lower = 2)
    checkmate::check_class(BPPARAM, "BiocParallelParam")
    checkmate::assertLogical(aggregateRecords)
}
