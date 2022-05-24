## ----knitr_init, echo = FALSE, cache = FALSE, results = 'hide', warning=FALSE, message=FALSE----
suppressPackageStartupMessages(library(knitr))
suppressPackageStartupMessages(library(katdetectr))

## ---- label = 'Importing genomic variants', eval = TRUE-----------------------
# Genomic variants stored within the VCF format.
pathToVCF <- system.file(package = "katdetectr", "extdata/CPTAC_Breast.vcf")

# Genomic variants stored within the MAF format.
pathToMAF <- system.file(package = "katdetectr", "extdata/APL_primary.maf")

# In addition, we can generate synthetic genomic variants with interjected kataegis regions
# using generateSyntheticData(). This will output a VRanges object.
syntheticData <- generateSyntheticData(nBackgroundVariants = 2500, nKataegisFoci = 1)

## ---- label = 'Detection of kataegis foci', eval = TRUE-----------------------
# Detect kataegis foci within the given VCF file.
kdVCF <- detectKataegis(genomicVariants = pathToVCF)

# # Detect kataegis foci within the given MAF file.
# As this file contains multiple samples, we set aggregateRecords = TRUE.
kdMAF <- detectKataegis(genomicVariants = pathToMAF, aggregateRecords = TRUE)

# Detect kataegis foci within our synthetic data.
kdSynthetic <- detectKataegis(genomicVariants = syntheticData)

## ---- label = 'Overview of KatDetect objects', eval = TRUE--------------------
summary(kdVCF)
print(kdVCF)
show(kdVCF)

# Or simply:
kdVCF

## ---- label = 'Retrieving information from KatDetect objects', eval = TRUE, results = 'hide'----
# Processed genomic variants used as input for changepoint detection.
getGenomicVariants(kdVCF)

# GRanges containing the segments as derived from changepoint detection.
getSegments(kdVCF)

# GRanges containing segments designated as putative kataegis foci.
getKataegisFoci(kdVCF)

# Supplementary information.
getInfo(kdVCF)

## ---- label = 'Visualisation of KatDetect objects', eval = TRUE---------------
rainfallPlot(kdVCF)

# With showSegmentation, the detected segments (changepoints) as visualized with their mean IMD.
rainfallPlot(kdMAF, showSegmentation = TRUE)

# With showSequence, we can display specific chromosomes or all chromosomes in which a putative kataegis foci has been detected.
rainfallPlot(kdSynthetic, showKataegis = TRUE, showSegmentation = TRUE, showSequence = "Kataegis")

## ---- label = 'Session Information', eval = TRUE------------------------------
utils::sessionInfo()

