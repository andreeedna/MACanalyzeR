# MACanalyzeR 1.1.0

## Breaking changes

* `MacPolarizeR()` has a new, reproducible algorithm. Each cell gets M1 and M2
  scores relative to the whole dataset; cells are classified by k-means
  (multiple starts, fixed seed) on the polarization balance `z(M1) - z(M2)`.
  The previous version clustered marker-gene expression with a single random
  start, so its results changed at every run and were often driven by overall
  activation rather than polarization. New arguments: `seed`, `nstart`,
  `n.boot`. A new `Mac.stability` metadata column gives the fraction of
  marker-gene bootstrap replicates agreeing with each assignment.
* `PathAnalyzeR()`: with three or more groups, the `pval` column is now a
  Friedman test (genes as blocks) across all groups. Before, the test only
  looked at the first group. Results for two groups are unchanged. A `padj`
  column (Benjamini-Hochberg) was added.
* `MitoBalance()`, `GliOxBalance()` and `macSpectrum()` now group cells by the
  default grouping of the object (`id` of `CreateMacObj()`) instead of the
  hard-coded `"Sample"`/`"Cluster"` columns, which often did not exist.
* `PathViolin(intercept = )` is now logical (`TRUE`/`FALSE`).

## New features

* `GeneBalance()`: per-cell balance between any two gene sets (the general
  version of `MitoBalance()` and `GliOxBalance()`).
* `MacSpectrumIndex()`: stores the MPI, AMDI and MDFI indices in the metadata.
* `PathCreatoR()` gained `pattern`, `min.size` and `max.size`; lists the
  available collections when called without `subcollection`; matching is
  case-insensitive, accepts the last part of the name (`"GO:BP"`) and the
  human/mouse codes interchangeably (`"C5:GO:BP"` for mouse).
* `PathCart()` and `PolCart()` show scores on a log2 scale by default
  (`scale = TRUE`), with `log2(...)` axis labels and reference lines at 0;
  `scale = FALSE` restores the original scale.
* `PathPlot()` accepts several pathways (one panel each) and a `meta` argument.
* `PathCart()`, `PathPlot()` and `PathViolin()` accept pathway names as well as
  indices.
* `PathHeat(adjust = TRUE)` filters on adjusted p-values.
* `CreateMacObj()` gained `assay` (default: the default assay of the Seurat
  object) and `sparse` (default `TRUE`): sparse matrices are no longer
  converted to dense, greatly reducing memory use.
* `FoamSpotteR(max.missing = )` tolerates a few missing model genes (set to 0
  with a warning) instead of failing.
* `MitoBalance()`, `GliOxBalance()` gained `split.by`; `ncol` is now used.
* A `show()` method summarises the object.
* `PathAnalyzeR()` and `MitoScanneR()` gained `verbose`.

## Bug fixes

* `PathCreatoR()` always failed: the bundled collections were saved under the
  wrong object names (`msgidb_mouse`, `msigdb_human`). They are now
  `gsea_mouse` and `gsea_human`, and the files are about 3 times smaller.
* `PathAnalyzeR()` and `MitoScanneR()` failed with a single gene set and when
  a gene set had only one expressed gene; they are now about 5-15 times faster
  (vectorised scoring) with identical scores.
* `PathPlot()` failed with more than one pathway.
* `MacPolarizeR()`: when the same cluster had the highest M1 and M2 score, the
  labels could be assigned to the wrong cluster.
* `MitoBalance()` silently returned only missing values when the
  mitochondrially encoded genes had been removed during QC; it now stops with
  an informative error.
* Cells with a score of 0 produced `NaN` balance values that were not removed.
* `CreateMacObj()`: reductions with more than two dimensions (e.g. PCA) were
  plotted on their last two dimensions; embeddings and metadata are now
  aligned to the cells of the expression matrix; Seurat v5 objects with split
  `data` layers give a clear error.
* `FoamSpotteR()` duplicated its columns when run twice.
* `PolCart()` ignored `fill`, `a` and `txt.size`.
* `MacBarplot()` failed with a single group; `PathHeat()`/`HeatMito()` failed
  when the colour range was zero and could duplicate rows.
* `macSpectrum()` reported a wrong error message for an invalid `mode`.
* `HeatMito()` assumed exactly 8 MitoCarta categories.
* Plotting functions no longer index columns by position.
* Clearer error messages when a module has not been run or a metadata column
  does not exist.

## Other

* `data/mono_macs_human.rds` moved to `inst/extdata/` (`.rds` files are not
  allowed in `data/`): use
  `system.file("extdata", "mono_macs_human.rds", package = "MACanalyzeR")`.
* Complete help pages for all functions and datasets.
