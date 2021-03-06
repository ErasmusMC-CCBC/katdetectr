# Internal - Import supplied genomic variants. ----

.importGenomicVariants <- function(x, aggregateRecords){

    # Input validation ----

    checkmate::assert(
        checkmate::checkClass(x, "VRanges"),
        checkmate::checkAccess(x, access = 'r')
    )

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
            base::stop('Provide a VRanges object or path to MAF or VCF file (extension: vcf, vcf.gz, maf, maf.gz)')
        )}

    # Check multiple samples. ----
    if(base::length(base::unique(Biobase::sampleNames(genomicVariants))) != 1){
        if(!aggregateRecords){
            base::stop('Your genomic variant data contains multiple sample names, select only a single sample to avoid overlapping mutations!
            I.e. pre-filter as a VRanges or set aggregateRecords = TRUE')
        }else{
            Biobase::sampleNames(genomicVariants) <- base::paste(base::unique(Biobase::sampleNames(genomicVariants)), collapse = ', ')
            genomicVariants <- IRanges::unique(genomicVariants)
        }
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

    # convert maf to VRanges
    vr <- maf@data |>
        tibble::as_tibble() |>
        # rename columns
        dplyr::rename(dplyr::any_of(columnNames)) |>
        # coerce to GRanges
        GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE) |>
        # coerse to VRanges
        VariantAnnotation::makeVRangesFromGRanges()

    # change seqlevel style
    GenomeInfoDb::seqlevelsStyle(vr) <- 'UCSC'

    return(vr)
}
