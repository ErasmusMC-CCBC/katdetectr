# katdetectr

![license](https://img.shields.io/badge/license-GPL--3-blue.svg) [![GitHub issues](https://img.shields.io/github/issues/ErasmusMC-CCBC/katdetectr.svg)](https://github.com/ErasmusMC-CCBC/katdetectr/issues) ![rversion](https://img.shields.io/badge/R%20version-%3E4.2.0-lightgrey.svg)

## Introduction

`katdetectr` is an *R* package for the detection, characterization and visualization of localized hypermutated regions, often referred to as *kataegis*.

Please see the [Application Note](https://www.biorxiv.org/content/10.1101/2022.07.11.499364v3) (under submission) for additional background, details and performance evaluations of `katdetectr`.

The general workflow of `katdetectr` can be summarized as follows:

1. Import of genomic variants; VCF, MAF or VRanges objects.
2. Detection of kataegis foci.
3. Visualization of segmentation and kataegis foci.

Please see the vignette for an overview of the workflow in a step-by-step manner on publicly-available datasets which are included within this package.

## Citation

Not yet published.

## Installation

### BioConductor

```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("katdetectr")
```

### Development version

```R
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

devtools::install_github("https://github.com/ErasmusMC-CCBC/katdetectr")
```
