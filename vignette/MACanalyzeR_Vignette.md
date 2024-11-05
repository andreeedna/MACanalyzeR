---
layout: page
title: MACanalyzeR
description: MACanalyzer vignette
img: #assets/img/12.jpg
importance: 1
category: Vignettes
---
-----------------------------------

# MACanalyzeR
*Macanalyzer* is an R package designed and developed for characterizing macrophages in single-cell RNA-seq data (scRNAseq). It provides three modules that enable in-depth analysis of this myeloid component.

**FoamSpotteR**: This module investigates the foamingness of macrophages. It is based on a machine learning model trained to recognize foamy macrophages from atherosclerotic plaques.

**MacPolarizeR**: This module investigates macrophage polarization. It is based on a clustering method built on pro- and anti-inflammatory genes.

**PathAnalyzeR**: This module calculates a score for each condition or single cell to allow comparison of pathways.

The package is available on GitHub at [andreeedna/MACanalyzeR](https://github.com/andreeedna/MACanalyzeR):

## Cite MACanalyzeR
If you use *MACanalyzeR* don't forget to cite the relative paper:
>**Characterization of Monocyto Macrophages in Obese Brown Adipose Tissue** Ninni A., Zaccaria F., Verteramo L., Sciarretta F. and Lettieri-Barbato D.

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

## Create MACanalyzeR object
*MACanalyzeR* employs a package-specific R object to store both the data and the results of the analysis. Given the widespread adoption of Seurat in single-cell analysis, there is a convenient function to instantiate an object directly from a Seurat one. It is *highly recommended* to utilize a preprocessed `Seurat` object containing only clustered monocyte/macrophages and desired metadata as input to ensure maximum compatibility.  
>MACanalyzeR is compatible with Seurat v4 and v5

The `CreateMacObj()` function creates a MACanalyzeR object from a Seurat object. It extracts the information from the Seurat metadata. It supports the following parameters:
- `seurat_obj`: the Seurat object;

-   `sample`: column name in the Seurat object’s metadata containing the experimental design for the differential analysis;

-   `cluster`: the column name in the Seurat object’s metadata storing cluster information or cell annotations;

-   `reduction`:  the name of the reduction utilized for cell visualization;
``` r
seurat_obj <- readRDS("vignette/mac.rds")
head(seurat_obj@meta.data)

mac_obj <- CreateMacObj(seurat_obj = seurat_obj,
						sample = "Sample",
						cluster = "Cells",
						reduction = "umap"
						)
```
    ##	Seurat version: 5.0.1
    ##	Warning in asMethod(object): sparse->dense coercion: allocating vector of size 1.3 GiB

The MACanalyzeR object will have four main metadata:
- Sample: The experimental design
- Cluster: The macrophage subtypes (or any other annotation)
- Foam: The FoamSpotteR() classification
- Mac: The MacPolarizeR() classification

> In the MACanalyzeR functions, the `plot.by` and `meta` flags take one of these four arguments as input (the name used in the Seurat object is not preserved!)  
### Visualize MACanalyzeR object
The MACanalyzeR object can be plotted with `MacPlot()` function. This function creates a 2D scatter plot for visualizing single-cell. Each point in the plot represents a cell, positioned based on its embedding coordinates obtained from a dimensionality reduction technique.
```
MacPlot(mac_obj,
		plot.by="Cluster",
		pt.size=2,
		txt.size=17)
```
![]()
The MacPlot can be also splitted by one of the metadata
```
MacPlot(mac_obj, split.by="Sample")
```
![]()

## *FoamSpotteR*: Foaming-like Macrophage Classification
The *FoamSpotteR* module enables the identification of cells exhibiting foam-like characteristics within scRNAseq data. This module employs a `randomforest` binary classifier on the data, assigning a label (**fMAC+** for foam-like macrophage and **fMAC-** for non foam-like macrophage) to each cell. Additionally, it computes a score, `FoamDEX`, for each cell, representing the probability of being a foam-like macrophage.  

> The use of this module is recommended on myeloid cells, as foam cells are a typical phenotype of this cell type under certain conditions.

``` r
mac_obj <- FoamSpotteR(mac_obj)
```
The details regarding the `FoamDEX` score and the prediction results are stored in the MACanalyzeR object, at `mac_obj@FoamSpotteR`.

### Visualize FoamSpotteR results
To visualize the *FoamSpotteR* prediction, *MACanalyzeR* integrates a series of plotting functionalities. Here’s a guide on how to plot the prediction results.
A basic approach for visualizing fMAC+/fMAC- classification is through the `MacPlot()` function:
``` r
MacPlot(mac_obj, plot.by = "Foam")
```
![]()

We can improve the visualization of the prediction by plotting the FoamDEX on the plot. This functionality is enabled via the `FoamPlot()` function:
``` r
FoamPlot(mac_obj)
```
![]()
> The FoamDEX is a probability score! It is a value between 0 and 1

The `FoamPlot()` function offers functionality similar to `MacPlot()`, enabling customization of point size, font size, and data splitting based on a chosen metadata.
``` r
FoamPlot(mac_obj,
		 split.by = "Sample",
		 shade="magma",
		 txt.size= 17,
		 pt.size=1.5
		 )
```
![]()

An additional feature to visualize the distribution of FoamDEX values. The `FoamLine()` function creates density plot for either individual samples or clusters. This graphical representation allows users to easily compare how FoamDEX levels are spread across different conditions, providing insights into potential variations.
``` r
FoamLine(mac_obj,
		 plot.by = "Sample",
		 fill=T,
		 a=0.05
		 )
```
![]()

## *MacPolarizeR*: Macrophages Polarization Analysis
The *MacPolarizeR* module classify monocytes and macrophages based on their polarization phenotype. This module operates by clustering cells according to the expression of genes associated with M1 and M2 polarization and cluster cells into 3 different class: inflammatory (*M1-like*), healing (*M2-like*), and transitional (*M0-like*).
In addition to classification, *MacPolarizeR* calculates an M1 and M2 polarization score for each cell based on the expression of the gene markers. This score system facilitates a more comprehensive and nuanced characterization of macrophages polarization state.
``` r
mac_obj <- MacPolarizeR(mac_obj)
```
The outupt of the cell polarization clustering are stored within the metadata of the MACanalyzer object accessible at `mac_obj@MetaData$Mac` and the polarization score is accessible at `mac_obj@MacPolarizeR`.

### Visualize MacPolarizeR result
A basic way for visualizing MacPolarizeR classification is through the `MacPlot()` function:
> MacPolarizeR implements preset colors for the three polarization classes
> Inflammatory - #C5283D, Transitional - #E9724C", Healing - #FFC857
``` r
MacPlot(mac_obj, plot.by = "Mac")
```
![]()

``` r
MacPlot(mac_obj, plot.by = "Mac", split.by = "Sample")
```
![]()

It’s also possible to represent the proportion of the 3 phenotypes for metadata through barplots:
``` r
MacBarplot(mac_obj, plot.by = "Sample")
```
![]()

or using a radarplot:
``` r
MacRadar(mac_obj, plot.by = "Sample")
```
![]()

*Will be implemented a function for polarization score plotting
## *PathAnalyzeR*: Pathway Enrichment Analysis
The *PathAnalyzeR* module enables the enrichment analysis of cells based on gene expression. The function `PathAnalyzeR()` performs an enrichment analysis for each cell individually and computes a score for each cluster or sample. The technique for calculating the enrichment score for each gene set is elaborated in the MACanalyzeR paper. This function allows the following arguments:

- `mac_obj`: the MACanalyzeR object

-   `meta`: the metadata used to perform the patway analysis

-   `pathway`:  a list of pathways of interest. Each element must be a vector containing the genes involved in a pathway. These vectors must be named with the corresponding pathway names. By default, it uses a subset of KEGG pathways of interest in biological processes or cellular components of macrophages.

> A collection of KEGG, GO and Reactome for mouse (*Mus musculus*) is already incorporated into the package. However, users have the flexibility to utilize any category extracted from other databases by providing a list of genes with a named pathway associated.
``` r
mac_obj <- PathAnalyzeR(mac_obj)
```
PathAnalyzeR allows you to define custom pathways for analysis. These pathways can be:
- Subsets of existing databases like KEGG, Gene Ontology (GO), and Reactome pathways (which are already implemented as lists in MACanalyzeR);
- User-created pathways: you can build your own pathway by specifying a list of genes;
``` r
# subset of KEGG, GO or REACTOME
ptw <- kegg[c(1:10)]
ptw <- go[c(1:10)]
ptw <- reactome[c(1:10)]

# user-created pathways
ptw <- list("path_name1" = c(...),
			"path_name2" = c(...),
			...
			)

# PathAnalyzeR performed on Cluster
mac_obj <- PathAnalyzeR(mac_obj, pathway=ptw, meta="Cluster")
```
### Visualize PathAnalyzeR results

To visualize the results of *PathAnalyzeR* enrichment, *MACanalyzeR* provides various functions for both per-metadata and single-cell analysis.

`PathHeat()` function plots the results of *PathAnalyzeR* metadata analysis as an heatmap that allows a general view of the enriched pathways between the specified groups:
``` r
PathHeat(mac_obj, plot.by = "Sample")
PathHeat(mac_obj, plot.by = "Cluster")
```
![]()
> Before plotting PathHeat for a specific metadata, you must first launch `PathAnalyzeR()` function with the appropriate `meta` flag corresponding to that metadata.

In *PathAnalyzeR* enrichment single-cell visualization, all functions use pathway identifiers in form of numbers (subsequent functions will also reference pathways using these numbers). To establish this association between pathways and their corresponding numbers, the `PathDisplay()` function is implemented.
``` r
PathDisplay(mac_obj, meta = "Sample")
```
> In the next functions, the `pathway` flag will always correspond to the pathways' number you choose to visualize using the `PathDisplay()` function. For educational purposes, we'll focus on pathways related to OXPHOS (#1) and Glycolysis (#2).

PathAnalyzeR implement a novel and unconventional visualization approach called Cartplot, which utilizes 2D density plots. The `PathCart()` function allows to simultaneously observe the expression distribution of two pathways at the single-cell level. This unique visualization provides valuable insights into the coordinated upregulation of pathways within the context of the analyzed metadata.
``` r
PathCart(mac_obj,
		 plot.by = "Sample",
		 pathway = c(1, 2),
		 fill = T
		 )
```
![]()
> Is available also single pathway density plot, specifying one number in the `pathway` flag

Through the `PathPlot()` function, we can visualize the enrichment score for each cell independently within a MacPlot. This enables the visualization of pathway enrichment patterns across the entire cell embedding.
``` r
## base PathPlot
PathPlot(mac_obj, pathway = c(1,2))

## splitted PathPlot
PathPlot(mac_obj,
				 pathway = c(1,2),
				 split.by="Sample"
				 )
```
![]()

Another interesting way to visualize the expression of single-cell pathways is with violin plots. The `PathViolin()` function allows you to visualize the distribution of expression for individual pathways and enables statistical comparisons between these expressions across different metadata.
``` r
PathViolin(mac_obj,
					 plot.by = "Sample",
					 pathway = c(1,2))
```

![]()

## *MitoScanneR*: Mitochondrial Characterization

The *MitoScanneR* module of MACanalyzeR is dedicated to the analysis of the oxidative and mitochondrial profile. This module allows the identification of cells with high oxidative and mitochondrial activity. It is based on PathAnalyzeR pathway score on MitoCarta.
> This version of MACanalyzeR uses pathway from MitoCarta 3.0
``` r
mac_obj <- MitoScanneR(mac_obj, meta = "Sample")
```
### Visualize MitoScanneR results
`HeatMito()` function allows to visualize the oxidative and mitochondrial activity of cells in a heatmap for a series the different MitoCarta pathways. The `plot.by` argument specifies the group to be analyzed (e.g., “Sample” or “Cluster”), passed to `MitoScanner`.
``` r
HeatMito(mac_obj, plot.by = "Sample")
```
![]()

The `MitoBalance()` function generates a violin plot to visualize potential shifts in the expression levels of mitochondrial and mitonuclear genes. This plot depicts the ratio between the expression of these two gene types.
``` r
MitoBalance(mac_obj, plot.by = "Sample")
```
![]()

The `GliOxBalance()` function determine whether there is a shift in the expression of genes related to glycolytic or oxidative phosphorylation.
``` r
GliOxBalance(mac_obj, plot.by = "Sample")
```
![]()
