.determineIMDcutoffValues <- function(IMDcutoff, genomicVariantsAnnotated, segments) {
    if (is.numeric(IMDcutoff)) {
        IMDcutoffValues <- IMDcutoff
    } else if (is.function(IMDcutoff)) {
        IMDcutoffValues <- IMDcutoff(genomicVariantsAnnotated, segments)
    } else {
        IMDcutoffValues <- 1000
    }

    return(IMDcutoffValues)
}
