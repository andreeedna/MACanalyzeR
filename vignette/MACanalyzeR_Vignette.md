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
*MACanalyzeR* is an R package for the in-depth characterization of monocytes and macrophages in single-cell RNA-seq (scRNA-seq) data. Starting from a clustered `Seurat` object (mouse or human), it provides four analysis modules, each with its own plotting functions:

| Module | Question it answers | Main function | Plots |
|---|---|---|---|
| **FoamSpotteR** | Which cells are foam-like (lipid-laden) macrophages? | `FoamSpotteR()` | `FoamPlot()`, `FoamLine()` |
| **MacPolarizeR** | Are cells inflammatory, transitional or healing? | `MacPolarizeR()` | `MacBarplot()`, `MacRadar()`, `PolCart()` |
| **PathAnalyzeR** | Which pathways are active in each group and in each cell? | `PathAnalyzeR()` | `PathHeat()`, `PathCart()`, `PathPlot()`, `PathViolin()` |
| **MitoScanneR** | What is the mitochondrial and metabolic profile? | `MitoScanneR()` | `HeatMito()`, `MitoBalance()`, `GliOxBalance()`, `GeneBalance()` |

In addition, `macSpectrum()` places cells on the MacSpectrum and AtheroSpectrum activation maps.

The package is available on GitHub at [andreeedna/MACanalyzeR](https://github.com/andreeedna/MACanalyzeR).

## Cite MACanalyzeR
If you use *MACanalyzeR* don't forget to cite the paper:
>**MACanalyzeR scRNAseq analysis tool reveals PPARγ<sup>HIGH</sup>/GDF15<sup>HIGH</sup> lipid-associated macrophages facilitate thermogenic expansion in BAT** Ninni A., Zaccaria F., Verteramo L., Sciarretta F. *et al.* & Lettieri-Barbato D. *Nature Communications* 16, 5063 (2025). [doi:10.1038/s41467-025-60295-2](https://doi.org/10.1038/s41467-025-60295-2)

## Installation
Install the dependencies from CRAN:
``` r
install.packages(c("ggplot2", "ggpubr", "gridExtra", "pheatmap",
                   "caret", "randomForest", "Matrix", "SeuratObject"))
```
Then install the development version of `MACanalyzeR` from GitHub with `devtools`:
``` r
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

devtools::install_github("andreeedna/MACanalyzeR")
library(MACanalyzeR)
```

-----------------------------------

## 1. Create the MACanalyzeR object
*MACanalyzeR* stores data and results in its own object. The `CreateMacObj()` function builds it from a `Seurat` object (v4 or v5). It is *highly recommended* to use a preprocessed (normalized, reduced and clustered) `Seurat` object containing **only monocytes/macrophages**, with all genes (not only the variable features).

Arguments:
- `seurat_obj`: the Seurat object;
- `id`: metadata column with clusters or cell annotations. It becomes the **default grouping** used by the `plot.by` and `meta` arguments of all the other functions (default `"seurat_clusters"`);
- `reduction`: reduction used for plotting, e.g. `"umap"` or `"tsne"` (default `"umap"`). Only the first two dimensions are kept;
- `org`: organism, `"mm"` (mouse, default) or `"hs"` (human). It selects the organism-specific models and gene sets;
- `assay`: assay to use (default: the default assay of the Seurat object, e.g. `"RNA"` or `"SCT"`);
- `sparse`: keep the expression matrix sparse to save memory (default `TRUE`).

``` r
seurat_obj <- readRDS("mac.rds")
head(seurat_obj@meta.data)

mac_obj <- CreateMacObj(seurat_obj,
                        id = "Cells",
                        reduction = "umap",
                        org = "mm")
mac_obj
```
    ##  An object of class MACanalyzeR (v1.1.0)
    ##  <n> features across <n> cells - Mus musculus
    ##  Default grouping: Cells
    ##  FoamSpotteR: not run | MacPolarizeR: not run

> With Seurat v5, if the `data` layer is split (e.g. after integration), join it first with `seurat_obj <- SeuratObject::JoinLayers(seurat_obj)`.

**Output.** The object contains the normalized expression (`mac_obj@CountData`), the embedding (`mac_obj@Reduction`) and the cell metadata (`mac_obj@MetaData`). All the metadata columns of the Seurat object are kept **with their original names** (in this example `Sample`, the experimental condition, and `Cells`, the macrophage subtypes): any of them can be passed to `plot.by`, `split.by` and `meta`. The modules add new columns (`Foam`, `fMAC+`, `Mac`, ...) and fill the slots `@MacPolarizeR`, `@PathAnalyzeR` and `@MitoScanneR`.

A human example dataset of monocytes/macrophages is shipped with the package:
``` r
human_obj <- readRDS(system.file("extdata", "mono_macs_human.rds", package = "MACanalyzeR"))
```

### Visualize the object: `MacPlot()`
`MacPlot()` draws the cells on the embedding, coloured by any metadata column.

- `plot.by`: column used for the colours (default: the `id` grouping);
- `split.by`: optional column used to split the plot into panels (with `ncol` columns);
- `pt.size`, `txt.size`: point and text size;
- `col`: optional vector with one colour per group.

``` r
MacPlot(mac_obj)                                   # coloured by the default grouping
MacPlot(mac_obj, plot.by = "Sample", pt.size = 1)  # coloured by condition
MacPlot(mac_obj, plot.by = "Cells", split.by = "Sample", ncol = 2)
```
All plotting functions return a `ggplot` object (except heatmaps and `PathViolin()`), so they can be further customized with `ggplot2` and saved with `ggsave()`.

-----------------------------------

## 2. *FoamSpotteR*: foam-like macrophages
*FoamSpotteR* identifies cells with a foam-like (lipid-laden) phenotype. A random forest classifier, trained on foamy and non-foamy macrophages from atherosclerotic plaques, assigns each cell to **fMAC+** (foam-like) or **fMAC-** and computes the **FoamDEX**, the probability (0-1) of being a foam-like macrophage.

``` r
mac_obj <- FoamSpotteR(mac_obj)
```
    ##  FoamSpotteR: <n> fMAC+ and <n> fMAC- cells

- `max.missing`: maximum fraction of model genes that may be missing from the data (default `0.1`); missing genes are set to 0 with a warning. Use `0` to require all of them.

**Output.** Three new columns in `mac_obj@MetaData`:
- `Foam`: predicted class (`fMAC+` / `fMAC-`);
- `fMAC+`: the FoamDEX;
- `fMAC-`: `1 - FoamDEX`.

``` r
table(mac_obj@MetaData$Sample, mac_obj@MetaData$Foam)
```

### Visualize FoamSpotteR results
The classification can be shown on the embedding with `MacPlot()`:
``` r
MacPlot(mac_obj, plot.by = "Foam")
```

`FoamPlot()` colours each cell by its FoamDEX (cells with the highest values on top; colour scale fixed between 0 and 1). `shade` sets the viridis palette (`"magma"`, `"inferno"`, `"plasma"`, `"viridis"`, `"cividis"`, `"rocket"`, `"mako"`, `"turbo"`).
``` r
FoamPlot(mac_obj)
FoamPlot(mac_obj, split.by = "Sample", shade = "magma", pt.size = 1.5)
```

`FoamLine()` compares the FoamDEX distribution between groups with one density curve per group (curves scaled to the same height). Use `fill = FALSE` for lines only, `a` for the transparency and `col` for custom colours.
``` r
FoamLine(mac_obj, plot.by = "Sample")
FoamLine(mac_obj, plot.by = "Cells", split.by = "Sample", fill = FALSE)
```

-----------------------------------

## 3. *MacPolarizeR*: macrophage polarization
*MacPolarizeR* classifies cells along the M1 (pro-inflammatory) / M2 (anti-inflammatory) axis:

1. each cell gets an **M1 and an M2 score** from curated marker genes (1 = average of the dataset);
2. the **polarization balance** `z(M1) - z(M2)` is computed, which removes the activation shared by both programs;
3. cells are split by k-means into three states: **Inflammatory** (M1-like), **Transitional** and **Healing** (M2-like);
4. marker genes are resampled to measure how **stable** each assignment is.

``` r
mac_obj <- MacPolarizeR(mac_obj)
```
    ##  MacPolarizeR: Inflammatory = <n>, Transitional = <n>, Healing = <n> | mean stability <x>

- `seed`: random seed, so the result is identical at every run (default `42`); the RNG state of your session is not changed;
- `nstart`: number of k-means starts (default `25`);
- `n.boot`: number of bootstrap replicates for the stability (default `20`; `0` to skip).

**Output.**
- `mac_obj@MetaData$Mac`: polarization state of each cell;
- `mac_obj@MetaData$Mac.stability`: fraction of replicates confirming the assignment (low values = cells between two states);
- `mac_obj@MacPolarizeR$SingleCell`: M1 and M2 scores of each cell;
- `mac_obj@MacPolarizeR$Total`: mean M1 and M2 scores of each state.

``` r
table(mac_obj@MetaData$Sample, mac_obj@MetaData$Mac)
summary(mac_obj@MetaData$Mac.stability)

# keep only confidently assigned cells
confident <- mac_obj@MetaData$Mac.stability >= 0.8
```

### Visualize MacPolarizeR results
`MacPlot()` uses a preset palette for the three states:
``` r
MacPlot(mac_obj, plot.by = "Mac")
MacPlot(mac_obj, plot.by = "Mac", split.by = "Sample")
```

`MacBarplot()` shows the percentage of each state in each group:
``` r
MacBarplot(mac_obj, plot.by = "Sample")
```

`MacRadar()` draws a circular barplot: one slice per group, with the number of cells coloured by state. `sort = TRUE` orders the groups by size.
``` r
MacRadar(mac_obj, plot.by = "Cells", sort = TRUE)
MacRadar(mac_obj, plot.by = "Cells", split.by = "Sample")
```
> In `MacBarplot()` and `MacRadar()`, `col` takes three colours for Inflammatory, Transitional and Healing.

`PolCart()` places each cell on the M1 score (x) vs M2 score (y) plane, as 2D density contours (`style = "density"`) or points (`style = "point"`). By default scores are shown on a **log2 scale** (`scale = TRUE`): reference lines at 0 split the plane into four quadrants (low/high M1, low/high M2), and +1/-1 mean twice/half the average. Use `scale = FALSE` for the original scores (reference lines at 1).
``` r
PolCart(mac_obj, plot.by = "Sample")
PolCart(mac_obj, plot.by = "Mac", style = "point")
PolCart(mac_obj, plot.by = "Sample", split.by = "Cells", fill = FALSE, scale = FALSE)
```

-----------------------------------

## 4. *PathAnalyzeR*: pathway scores
*PathAnalyzeR* scores gene sets both **per group** and **per cell**. For each pathway, the expression of each gene is divided by its average across the groups of `meta`, and the ratios are averaged weighting each gene by its variance. A score of **1 means average activity**, above 1 higher and below 1 lower activity (details in the paper).

Arguments:
- `pathway`: a named list of gene vectors (default: a curated selection of 26 macrophage-relevant KEGG pathways). A single gene vector is also accepted;
- `meta`: the metadata column defining the groups (default: the `id` grouping);
- `verbose`: show the progress bar.

``` r
mac_obj <- PathAnalyzeR(mac_obj, meta = "Sample")
```

Pathways with fewer than 5 detected genes (or fewer than 3 genes expressed in every group) are skipped and reported.

**Output.** `mac_obj@PathAnalyzeR[["Sample"]]`, a list of two data frames:
- `Total`: score of each pathway in each group, plus `pval` and `padj` (Benjamini-Hochberg). The p-value tests whether the pathway differs between groups (Wilcoxon signed-rank test with 2 groups, Friedman test with 3 or more);
- `SingleCell`: score of each pathway in each cell.

``` r
res <- mac_obj@PathAnalyzeR$Sample$Total
res[order(res$padj), ]
```

Results are stored by grouping, so you can run `PathAnalyzeR()` on several groupings and keep all of them. The plotting functions look them up by the same name: to plot by `"Cells"`, first run `PathAnalyzeR(mac_obj, meta = "Cells")`.

### Choosing the pathways
Any named list of gene vectors can be used. Several collections are included:
``` r
# subsets of the built-in KEGG, GO and Reactome collections (mouse)
ptw <- mouse_kegg[1:10]
ptw <- mouse_go[grep("MACROPHAGE", names(mouse_go))]
ptw <- mouse_reactome[1:10]

# MSigDB collections, available offline with PathCreatoR()
PathCreatoR("mm")                                   # list the available collections
ptw <- PathCreatoR("mm", "MH")                      # hallmark
ptw <- PathCreatoR("mm", "GO:BP", pattern = "lipid|foam", min.size = 10, max.size = 300)
ptw <- PathCreatoR("hs", "KEGG_LEGACY")             # human KEGG

# user-defined gene sets
ptw <- list("Lipid uptake"  = c("Cd36", "Lpl", "Msr1", "Scarb1", "Ldlr"),
            "LAM signature" = c("Trem2", "Gpnmb", "Lgals3", "Cd9", "Spp1"))

mac_obj <- PathAnalyzeR(mac_obj, pathway = ptw, meta = "Cells")
```
> `PathCreatoR()` matching is case-insensitive, accepts the last part of the name (`"GO:BP"`) and human codes for mouse (`"C5:GO:BP"`). Gene symbols must match the organism of the data.

### Visualize PathAnalyzeR results
`PathHeat()` shows the group-level scores as a heatmap (centred on 1), keeping only pathways with `pval` below the threshold (`adjust = TRUE` to filter on `padj`):
``` r
PathHeat(mac_obj, plot.by = "Sample")
PathHeat(mac_obj, plot.by = "Sample", pval = 0.01, adjust = TRUE)
```

The single-cell plots select pathways by **number** or by **name**. `PathDisplay()` lists the numbers:
``` r
PathDisplay(mac_obj, meta = "Sample")
```
    ##  Pathways used for the Sample pathway analysis:
    ##  1 - GLYCOLYSIS_GLUCONEOGENESIS
    ##  2 - CITRATE_CYCLE_TCA_CYCLE
    ##  ...
    ##  5 - OXIDATIVE_PHOSPHORYLATION

`PathCart()` draws the **Cartplot**: the single-cell scores of two pathways on the x and y axes, as 2D densities coloured by group, showing cells in which both pathways are co-upregulated. With one pathway it draws a 1D density. As in `PolCart()`, scores are on a **log2 scale** by default (`scale = TRUE`, reference lines at 0); use `scale = FALSE` for the original scores.
``` r
PathCart(mac_obj, pathway = c(5, 1), plot.by = "Sample")
PathCart(mac_obj, pathway = "OXIDATIVE_PHOSPHORYLATION", plot.by = "Sample", fill = FALSE)
PathCart(mac_obj, pathway = c(5, 1), plot.by = "Sample", split.by = "Cells", scale = FALSE)
```

`PathPlot()` colours the cells on the embedding by their score for one or more pathways (one panel each). `max.cutoff` caps extreme values so that a few cells do not compress the colour scale.
``` r
PathPlot(mac_obj, pathway = 5, meta = "Sample")
PathPlot(mac_obj, pathway = c(5, 1), meta = "Sample", max.cutoff = 3)
PathPlot(mac_obj, pathway = 5, meta = "Sample", split.by = "Sample")
```

`PathViolin()` compares the single-cell scores between groups with violin plots. `stat = TRUE` adds the p-values of all pairwise Wilcoxon tests, `intercept = TRUE` a reference line at 1.
``` r
PathViolin(mac_obj, pathway = c(5, 1), plot.by = "Sample", ncol = 2)
PathViolin(mac_obj, pathway = 5, plot.by = "Sample", stat = TRUE, intercept = TRUE)
```

-----------------------------------

## 5. *MitoScanneR*: mitochondrial profile
*MitoScanneR* applies the PathAnalyzeR scoring to the pathways of **MitoCarta 3.0**, organized in its main categories (OXPHOS, Metabolism, Mitochondrial central dogma, Protein import, Mitochondrial dynamics, ...).
``` r
mac_obj <- MitoScanneR(mac_obj, meta = "Sample")
```
**Output.** `mac_obj@MitoScanneR[["Sample"]]`: one data frame (pathways x groups) per MitoCarta category.
``` r
names(mac_obj@MitoScanneR$Sample)
head(mac_obj@MitoScanneR$Sample$OXPHOS)
```

### Visualize MitoScanneR results
`HeatMito()` draws one heatmap per MitoCarta category (`plot.by` must match the `meta` used in `MitoScanneR()`):
``` r
HeatMito(mac_obj, plot.by = "Sample")
```

### Metabolic balances
These functions compare two gene sets **in each cell** and plot `log10(score set1 / score set2)` as violin plots per group (0 = balanced; `intercept = TRUE` draws the line at 0). They do not need `MitoScanneR()` to be run first.

`MitoBalance()` compares the mitochondrially encoded OXPHOS subunits with the nuclear-encoded ones (positive = relative excess of mitochondrial subunits):
``` r
MitoBalance(mac_obj, plot.by = "Sample", intercept = TRUE)
```
> The 13 mitochondrial genes (`mt-` in mouse, `MT-` in human) are often removed during QC: keep them in the Seurat object to use `MitoBalance()`.

`GliOxBalance()` compares glycolysis with oxidative phosphorylation (positive = glycolytic, negative = oxidative):
``` r
GliOxBalance(mac_obj, plot.by = "Sample", intercept = TRUE)
GliOxBalance(mac_obj, plot.by = "Mac", split.by = "Sample")
```

`GeneBalance()` does the same for **any pair of gene sets**:
``` r
sel <- mouse_selected
GeneBalance(mac_obj,
            set1 = sel$GLYCOLYSIS_GLUCONEOGENESIS,
            set2 = sel$FATTY_ACID_METABOLISM,
            plot.by = "Sample", intercept = TRUE,
            ylab = "log10(glycolysis / FAO)")
```

-----------------------------------

## 6. MacSpectrum and AtheroSpectrum
`macSpectrum()` places each cell on two activation maps, based on its similarity to reference macrophage profiles:
- `mode = "macspectrum"`: **MPI** (polarization index; M1-like positive, M2-like negative) vs **AMDI** (activation-induced differentiation index);
- `mode = "atherospectrum"`: **MPI** vs **MDFI** (foaming index; foamy-like positive, control-like negative).

``` r
macSpectrum(mac_obj, plot.by = "Sample")
macSpectrum(mac_obj, mode = "atherospectrum", plot.by = "Cells", style = "point")
```

`MacSpectrumIndex()` stores the three indices in the metadata, to reuse or export them:
``` r
mac_obj <- MacSpectrumIndex(mac_obj)
summary(mac_obj@MetaData[, c("MPI", "AMDI", "MDFI")])
```

-----------------------------------

## Combining the modules
Every result is a metadata column, so modules can be combined freely, e.g. pathways in foam-like vs non-foam-like cells, or metabolism in each polarization state:
``` r
mac_obj <- PathAnalyzeR(mac_obj, meta = "Foam")
PathHeat(mac_obj, plot.by = "Foam")

mac_obj <- PathAnalyzeR(mac_obj, meta = "Mac")
PathViolin(mac_obj, pathway = "OXIDATIVE_PHOSPHORYLATION", plot.by = "Mac", stat = TRUE)

# export all per-cell results
write.csv(mac_obj@MetaData, "MACanalyzeR_cell_results.csv")
```

## Reference data
The gene sets and models used by the package can be inspected directly (see `?"MACanalyzeR-data"`):

| Object | Content |
|---|---|
| `mouse_selected`, `human_selected` | Default PathAnalyzeR pathways (macrophage-relevant KEGG) |
| `mouse_kegg`, `mouse_go`, `mouse_reactome` | Mouse KEGG, GO and Reactome pathways |
| `gsea_mouse`, `gsea_human` | MSigDB collections (via `PathCreatoR()`) |
| `mouse_MitoCarta`, `human_MitoCarta` | MitoCarta 3.0 pathways |
| `mouse_PolGenes`, `human_PolGenes` | M1 and M2 marker genes |
| `mouse_LAMprey`, `human_LAMprey` | FoamSpotteR random forest models |
