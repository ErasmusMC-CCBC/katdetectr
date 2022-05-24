validateInputdetectKataegis <- function(genomicVariants, minSizeKataegis, maxMeanIMD, test.stat, penalty, pen.value, minseglen, bpworkers, aggregateRecords){

    checkmate::assert(
        checkmate::checkClass(genomicVariants, 'VRanges'),
        checkmate::checkAccess(genomicVariants, access = 'r')
    )

    # PCF-specific parameters.
    checkmate::assertInt(minSizeKataegis, lower = 1)
    checkmate::assertNumeric(maxMeanIMD, 0)
    checkmate::assert_character(test.stat, pattern = 'Exponential|Empirical')
    checkmate::assert_character(penalty, pattern = 'BIC|Manual')
    checkmate::assertInt(pen.value, lower = 0)
    checkmate::assertInt(minseglen, lower = 2)
    checkmate::assertInt(bpworkers, lower = 1)
    checkmate::assertLogical(aggregateRecords)
}
