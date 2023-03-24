# Internal - Import supplied genomic variants. ----

.importGenomicVariants <- function(x, aggregateRecords, refSeq){

    # Check type of provided data and import accordingly. ----
    if(methods::is(x, 'VRanges')){
        genomicVariants <- x
    }else{
        # Coerce based on file extension.
        genomicVariants <- base::switch(
            ifelse(tools::file_ext(x) == 'gz', tools::file_ext(base::gsub('.gz', '', x)), tools::file_ext(x)),
            vcf = .coerceVCFtoVRanges(x),
            vcf.gz = .coerceVCFtoVRanges(x),
            maf = .coerceMAFtoVRanges(x),
            maf.gz = .coerceMAFtoVRanges(x),
            base::stop('Provide a VRanges object or path to MAF or VCF file (vcf, vcf.gz, maf, maf.gz)')
        )}

    # Check multiple samples
    if(base::length(base::unique(Biobase::sampleNames(genomicVariants))) != 1){
        if(!aggregateRecords){
            base::stop('Your genomic variant data contains multiple sample names. Provide data from a single sample or set aggregateRecords = TRUE')
        }else{
            Biobase::sampleNames(genomicVariants) <- base::paste(base::unique(Biobase::sampleNames(genomicVariants)), collapse = ', ')
            genomicVariants <- IRanges::unique(genomicVariants)
        }
    }

    # set seqlevel style to USCS
    GenomeInfoDb::seqlevelsStyle(genomicVariants) <- "UCSC"

    # check for non standard sequences
    seqLev <- GenomeInfoDb::seqlevelsInUse(genomicVariants)
    allStandardSeqlev <- seqLev %in% paste0("chr", c(1:22, "X", "Y", "M"))
    seqLevNonStandard <- seqLev[!allStandardSeqlev]
    seqLevInTib <- seqLev %in% names(refSeq)

    if(!all(allStandardSeqlev) & !all(seqLevInTib)){
        base::stop(base::paste0(c("Your genomic variant data contains the following non standard sequences:",
                                  seqLevNonStandard,
                                  "Please check the vignette on how to deal with non standard sequences."
        ), collapse = " "))
    }

    return(genomicVariants)
}


# Helper - Coerce VCF into VRanges. ----
.coerceVCFtoVRanges <- function(path = NULL){

    # Read vcf as VRanges.
    vr <- VariantAnnotation::readVcfAsVRanges(path)

    # Remove unwanted columns.
    GenomicRanges::mcols(vr) <- NULL

    # Change seqlevel style.
    GenomeInfoDb::seqlevelsStyle(vr) <- 'UCSC'

    return(vr)
}

# Helper - Coerce MAF into VRanges. ----
.coerceMAFtoVRanges <- function(path = NULL){

    # use read.maf function from maftools for import
    maf <- maftools::read.maf(maf = path, verbose = FALSE)

    # specify the correct column names in order to coerce maf to granges
    columnNames <- c(
        seqnames = 'Chromosome',
        start = 'Start_Position',
        end = 'End_Position',
        strand = 'Strand',
        ref = 'Reference_Allele',
        alt = 'Tumor_Seq_Allele2',
        sampleNames = 'Tumor_Sample_Barcode'
    )

    # a MAF object contains separate slots for the synonymous (silent) and the non-synonymous variants
    mafSilent <-  maf@maf.silent |>
        tibble::as.tibble()

    # convert maf to VRanges
    vr <- maf@data |>
        tibble::as_tibble() |>
        # combine synonymous and non-synonymous variants
        dplyr::bind_rows(mafSilent) |>
        # rename columns
        dplyr::rename(dplyr::any_of(columnNames)) |>
        # remove rows of which start or end is NA
        dplyr::filter(!is.na(.data$start) | !is.na(.data$end)) |>
        # coerce to GRanges
        GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE) |>
        # coerse to VRanges
        VariantAnnotation::makeVRangesFromGRanges()

    # change seqlevel style
    GenomeInfoDb::seqlevelsStyle(vr) <- 'UCSC'

    return(vr)
}
