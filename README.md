# MACanalyzeR
[![bioAxiv Article](https://img.shields.io/badge/Article-bioRxiv-red)](https://www.biorxiv.org/content/10.1101/2024.07.29.605588v1)
[![bioAxiv Article](https://img.shields.io/badge/GitHub-Matherials&Methods-blue)](https://github.com/andreeedna/MACanalyzeR_MaterialsMethods)

*Macanalyzer* is an R package designed and developed for characterizing macrophages in single-cell RNA-seq data (scRNAseq). It provides three modules that enable in-depth analysis of this myeloid component.

**FoamSpotteR**: This module investigates the foamingness of macrophages. It is based on a machine learning model trained to recognize foamy macrophages from atherosclerotic plaques.

**MacPolarizeR**: This module investigates macrophage polarization. It is based on a clustering method built on pro- and anti-inflammatory genes.

**PathAnalyzeR**: This module calculates a score for each condition or single cell to allow comparison of pathways.

The vignette is available on GitHub at [andreeedna.github.io/projects/](https://andreeedna.github.io/projects/): 

## Cite MACanalyzeR
If you use MACanalyzeR don't forget to cite the relative paper:
>**MACanalyzeR: scRNA-seq Analysis Tool Reveals PPARγHI Lipid-Associated Macrophages Facilitate Thermogenic Expansion in BAT** (2024) Ninni A., Zaccaria F., Verteramo L., Sciarretta F. [...] and Lettieri-Barbato D. 
>doi.org/10.1101/2024.07.29.605588 
## Installation
To ensure a correct functionality, it’s essential to set up the required dependencies. 
Below are the steps to install them:
``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("ggplot2", "ggpubr", "gridExtra", "pheatmap", "caret"))
```
The development version of `MACanalyzeR` can be installed using the `devtools` package:

``` r
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

devtools::install_github("andreeedna/MACanalyzeR")
```
The latest stable release can also be installed from the BioConductor repository \[NOT PUBLISHED YET\]:
``` r
BiocManager::install("MACanalyzeR")
```
## Contributions
*MACanalyzeR* is a product conceived, designed, and programmed by the tSNEland team. We would like to express our gratitude to the following individuals for their contributions:

-   **Fabio Zaccaria (PhD Student):** Conceptualization and Alpha Testing
-   **Luca Verteramo (Master's Student):** Development of the Human Version
-   **Francesca Sciarretta (PostDoc):** Experimental Validations
