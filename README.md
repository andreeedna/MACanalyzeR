
# MACanalyzeR

### Single-cell characterization of macrophages in R

[![Nature Communications](https://img.shields.io/badge/Nature%20Communications-2025-1f6f8b?style=flat-square)](https://doi.org/10.1038/s41467-025-60295-2)
[![bioRxiv](https://img.shields.io/badge/bioRxiv-preprint-b31b1b?style=flat-square)](https://www.biorxiv.org/content/10.1101/2024.07.29.605588v1)
[![Materials & Methods](https://img.shields.io/badge/GitHub-Materials%20%26%20Methods-24292e?style=flat-square&logo=github)](https://github.com/andreeedna/MACanalyzeR_MaterialsMethods)
[![Vignette](https://img.shields.io/badge/docs-vignette-4c9a2a?style=flat-square)](https://andreeedna.github.io/projects/)
<br>
![R](https://img.shields.io/badge/R-%E2%89%A5%204.0-276DC3?style=flat-square&logo=r&logoColor=white)
![Seurat](https://img.shields.io/badge/Seurat-v4%20%7C%20v5-6a4c93?style=flat-square)
![Organisms](https://img.shields.io/badge/organisms-mouse%20%7C%20human-e07a5f?style=flat-square)
![Version](https://img.shields.io/badge/version-1.1.0-lightgrey?style=flat-square)

---

**MACanalyzeR** is an R package for the in-depth characterization of macrophages in single-cell RNA-seq (scRNA-seq) data. Starting from a clustered `Seurat` object, it spots foam-like macrophages with a machine-learning classifier, assigns pro- and anti-inflammatory polarization states, scores pathways at single-cell and group level, and profiles mitochondrial and oxidative metabolism — with publication-ready plots for every step.

## 🆕 What's new in v1.1.0

| | Change |
|---|---|
| 🧭 **MacPolarizeR** | New reproducible algorithm: k-means on the M1/M2 polarization balance with a fixed seed and multiple starts. Same result at every run, plus a per-cell `Mac.stability` score |
| 📈 **PathCart / PolCart** | Scores shown on a **log2 scale** by default (`scale = TRUE`), with `log2(...)` axis labels and reference lines at 0 |
| 🛤️ **PathCreatoR** | Fixed (it always failed); lists the available MSigDB collections and filters gene sets by name and size |
| ⚖️ **GeneBalance** | New: per-cell balance between any two gene sets. `MitoBalance()` and `GliOxBalance()` now work with any grouping |
| 📊 **PathAnalyzeR** | 5-15x faster with identical scores; correct p-value with 3+ groups (Friedman test) and BH-adjusted `padj` |
| 💾 **Memory** | Sparse matrices are kept sparse; Seurat v5 layers and any assay supported |
| 🐞 **Bug fixes** | `PathPlot()` with several pathways, PCA reductions, duplicated FoamSpotteR columns, clearer error messages |

Full details in [NEWS.md](NEWS.md).

## ✨ Modules

| Module | What it does | How |
|---|---|---|
| 🫧 **FoamSpotteR** | Identifies foam-like macrophages (**fMAC+** / **fMAC−**) and gives each cell a **FoamDEX** score (probability of being foamy, 0–1) | Random forest classifier trained on foamy macrophages from atherosclerotic plaques |
| 🧭 **MacPolarizeR** | Classifies cells as Inflammatory, Transitional or Healing, with M1/M2 scores and a per-cell stability | Reproducible clustering of the M1/M2 polarization balance of curated marker genes |
| 🛤️ **PathAnalyzeR** | Scores pathways per cell and per condition/cluster, for easy comparison across groups | Built-in KEGG, GO and Reactome collections, or any custom gene set |
| ⚡ **MitoScanneR** | Profiles mitochondrial and oxidative activity, mito-nuclear and glycolysis/OXPHOS balance | PathAnalyzeR scoring on MitoCarta 3.0 pathways |

MACanalyzeR works with mouse and human data and keeps sparse matrices sparse, so large datasets fit in memory.

## 📦 Installation

MACanalyzeR depends on a few CRAN packages. Install them first:

```r
install.packages(c("ggplot2", "ggpubr", "gridExtra", "pheatmap",
                   "caret", "randomForest", "Matrix", "SeuratObject"))
```

Then install the development version from GitHub:

```r
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

devtools::install_github("andreeedna/MACanalyzeR")
```

> [!NOTE]
> A Bioconductor release is planned but not yet available.

## 🚀 Quick start

MACanalyzeR works best on a **preprocessed `Seurat` object (v4 or v5) containing only clustered monocytes/macrophages** and the metadata you want to compare (e.g. condition, genotype, cluster annotation).

```r
library(MACanalyzeR)

# 1. Create the MACanalyzeR object from Seurat
mac_obj <- CreateMacObj(seurat_obj,
                        id        = "seurat_clusters",  # metadata column with clusters / cell annotations
                        reduction = "umap",             # embedding used for plots
                        org       = "mm")               # "mm" (mouse) or "hs" (human)

MacPlot(mac_obj, split.by = "Sample")

# 2. Foam-like macrophages
mac_obj <- FoamSpotteR(mac_obj)
FoamPlot(mac_obj)                       # FoamDEX on the embedding
FoamLine(mac_obj, plot.by = "Sample")   # FoamDEX distribution per group

# 3. Polarization (reproducible: fixed seed by default)
mac_obj <- MacPolarizeR(mac_obj)
summary(mac_obj@MetaData$Mac.stability)   # per-cell stability of the assignment
MacPlot(mac_obj, plot.by = "Mac")
MacBarplot(mac_obj, plot.by = "Sample")
MacRadar(mac_obj, plot.by = "Sample")

# 4. Pathway scoring (default: macrophage-relevant KEGG subset)
mac_obj <- PathAnalyzeR(mac_obj, meta = "Sample")
PathHeat(mac_obj, plot.by = "Sample")
PathDisplay(mac_obj, meta = "Sample")   # pathway numbers used by the plots below
PathCart(mac_obj, pathway = c(1, 2), plot.by = "Sample")
PathViolin(mac_obj, pathway = c(1, 2), plot.by = "Sample")

# 5. Mitochondrial profile
mac_obj <- MitoScanneR(mac_obj, meta = "Sample")
HeatMito(mac_obj, plot.by = "Sample")
MitoBalance(mac_obj, plot.by = "Sample")
GliOxBalance(mac_obj, plot.by = "Sample")
```

Any MSigDB collection is available offline through `PathCreatoR()`, and custom gene sets can be passed to `PathAnalyzeR()` as a named list:

```r
hallmark <- PathCreatoR("mm", "MH")
mac_obj  <- PathAnalyzeR(mac_obj, pathway = hallmark, meta = "Sample")

my_pathways <- list(
  "Lipid uptake"   = c("Cd36", "Lpl", "Msr1"),
  "LAM signature"  = c("Trem2", "Gpnmb", "Lgals3", "Cd9")
)
mac_obj <- PathAnalyzeR(mac_obj, pathway = my_pathways, meta = "Cluster")
```

📖 The full step-by-step tutorial is in the **[vignette](https://andreeedna.github.io/projects/)**.

## 🧰 Function overview

| Step | Analysis | Visualization |
|---|---|---|
| Setup | `CreateMacObj()` | `MacPlot()` |
| FoamSpotteR | `FoamSpotteR()` | `FoamPlot()`, `FoamLine()` |
| MacPolarizeR | `MacPolarizeR()` | `MacBarplot()`, `MacRadar()`, `PolCart()` |
| PathAnalyzeR | `PathAnalyzeR()`, `PathCreatoR()` | `PathHeat()`, `PathDisplay()`, `PathCart()`, `PathPlot()`, `PathViolin()` |
| MitoScanneR | `MitoScanneR()` | `HeatMito()`, `MitoBalance()`, `GliOxBalance()`, `GeneBalance()` |
| Spectrum | `MacSpectrumIndex()` | `macSpectrum()` (`mode = "macspectrum"` or `"atherospectrum"`) |

Built-in reference data include KEGG, GO and Reactome collections, the MSigDB collections for mouse and human (retrieved with `PathCreatoR()`), MitoCarta 3.0, polarization and LAM gene sets for both mouse and human, plus a human–mouse ortholog dictionary.

## 📝 Citation

If you use MACanalyzeR in your work, please cite:

> Ninni A., Zaccaria F., Verteramo L., Sciarretta F., *et al.* & Lettieri-Barbato D.
> **MACanalyzeR scRNAseq analysis tool reveals PPARγ<sup>HIGH</sup>/GDF15<sup>HIGH</sup> lipid-associated macrophages facilitate thermogenic expansion in BAT.**
> *Nature Communications* 16, 5063 (2025). [doi:10.1038/s41467-025-60295-2](https://doi.org/10.1038/s41467-025-60295-2)

<details>
<summary>BibTeX</summary>

```bibtex
@article{ninni2025macanalyzer,
  title   = {MACanalyzeR scRNAseq analysis tool reveals PPAR$\gamma^{HIGH}$/GDF15$^{HIGH}$ lipid-associated macrophages facilitate thermogenic expansion in BAT},
  author  = {Ninni, Andrea and Zaccaria, Fabio and Verteramo, Luca and Sciarretta, Francesca and others and Lettieri-Barbato, Daniele},
  journal = {Nature Communications},
  volume  = {16},
  pages   = {5063},
  year    = {2025},
  doi     = {10.1038/s41467-025-60295-2}
}
```

</details>

The code and analyses used in the paper are available at [MACanalyzeR_MaterialsMethods](https://github.com/andreeedna/MACanalyzeR_MaterialsMethods).

## 👥 Team

MACanalyzeR was conceived, designed and developed by **Andrea Ninni** and the **tSNEland** team, with contributions from:

- **Fabio Zaccaria** — conceptualization and alpha testing
- **Luca Verteramo** — development of the human version
- **Francesca Sciarretta** — experimental validation

## 💬 Contact & issues

Found a bug or have a feature request? Please open an [issue](https://github.com/andreeedna/MACanalyzeR/issues).
For other questions, contact the maintainer at andrea.ninni@uniroma2.it.
