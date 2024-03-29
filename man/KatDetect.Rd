% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllClasses.R
\docType{class}
\name{KatDetect-class}
\alias{KatDetect-class}
\title{KatDetect-class: KatDetect objects}
\value{
(KatDetect) KatDetect object.
}
\description{
The katdetectr package introduces a new S4 object which stores all relevant information regarding kataegis detection.
}
\section{Slots}{

\describe{
\item{\code{kataegisFoci}}{(\link[GenomicRanges]{GRanges}): Contains all
annotated putative kataegis foci.}

\item{\code{genomicVariants}}{(\link[VariantAnnotation]{VRanges}): Contains all
processed and annotated genomic variants.}

\item{\code{segments}}{(\link[GenomicRanges]{GRanges}): Contains all segments
detected after changepoint analysis.}

\item{\code{info}}{(list): Contains some general information and model parameters used for kataegis detection.}
}}

\examples{

syntheticData <- generateSyntheticData()
kd <- detectKataegis(syntheticData)

getKataegisFoci(kd)
getGenomicVariants(kd)
getSegments(kd)
getInfo(kd)

}
\author{
Daan Hazelaar

Job van Riet
}
