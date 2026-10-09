#' Macrophage Spectrum Analysis
#'
#' @param object A Seurat object
#' @param mode Mode: "macspectrum" or "atherospectrum"
#' @param group.by Metadata column for coloring (default: from CreateMacObj)
#' @param split.by Metadata column for facetting
#' @param style Plot style: "point" or "density"
#' @param base.size Font size
#' @param cols Custom colors
#' @param ncol Columns for facetting
#' @param assay Assay (default: DefaultAssay(object))
#' @param layer Layer
#' @return A Seurat object with spectrum coordinates
#' @import ggplot2
#' @importFrom Seurat GetAssayData DefaultAssay
#' @export
macSpectrum <- function(
  object,
  mode = "macspectrum",
  group.by = NULL,
  split.by = NULL,
  style = "density",
  base.size = 17,
  cols = NULL,
  ncol = 3,
  assay = NULL,
  layer = "data"
) {
  if (!inherits(object, "Seurat")) {
    stop("Input must be a Seurat object")
  }

  mac_info <- object@misc$MACanalyzeR
  if (is.null(mac_info)) {
    stop("MACanalyzeR object not initialized. Please run 'CreateMacObj' first")
  }

  if (is.null(group.by)) {
    group.by <- mac_info$ident
  }
  if (is.null(object@meta.data[[group.by]])) {
    stop("Metadata column '", group.by, "' not found")
  }

  organism <- mac_info$organism
  if (is.null(assay)) {
    assay <- Seurat::DefaultAssay(object)
  }
  mtx <- Seurat::GetAssayData(object, assay = assay, layer = layer)

  calc_ms <- function(mac_mtx, organism) {
    MFDI_score <- function(mac_mtx, organism) {
      con_mean <- MACanalyzeR::con_mean_new
      foam_mean <- MACanalyzeR::foam_mean_new

      if (organism != "mm") {
        dict <- MACanalyzeR::human_mouse_dict
        names(con_mean) <- dict[names(con_mean), paste0(organism, "_symbol")]
        con_mean <- con_mean[!is.na(names(con_mean))]
        names(foam_mean) <- dict[names(foam_mean), paste0(organism, "_symbol")]
        foam_mean <- foam_mean[!is.na(names(foam_mean))]
      }

      MFI_genes <- intersect(names(foam_mean), rownames(mac_mtx))
      con_mean <- con_mean[MFI_genes]
      foam_mean <- foam_mean[MFI_genes]
      MFI_mtx <- as.matrix(mac_mtx[MFI_genes, ])

      mac_sigma <- sqrt(colSums(MFI_mtx^2) / length(MFI_genes))
      con_sigma <- sqrt(sum(con_mean^2) / length(con_mean))
      foam_sigma <- sqrt(sum(foam_mean^2) / length(foam_mean))

      con_pearson <- colSums((MFI_mtx / mac_sigma) * (con_mean / con_sigma)) / length(MFI_genes)
      foam_pearson <- colSums((MFI_mtx / mac_sigma) * (foam_mean / foam_sigma)) / length(MFI_genes)

      a <- 1; b <- 1; c <- 0
      d_sqr <- (a * foam_pearson + b * con_pearson + c)^2 / (a^2 + b^2)
      x_start <- -1; y_start <- 1; x_end <- 1; y_end <- -1

      l <- sqrt(pmax(0, (foam_pearson - x_start)^2 + (con_pearson - y_start)^2 - d_sqr))
      l_max <- sqrt(pmax(0, (x_end - x_start)^2 + (y_end - y_start)^2 - d_sqr))
      return((l / l_max) * 100 - 50)
    }

    M1_mean <- MACanalyzeR::M1_mean_new
    M2_mean <- MACanalyzeR::M2_mean_new
    M0_mean <- MACanalyzeR::M0_mean_new

    if (organism != "mm") {
      dict <- MACanalyzeR::human_mouse_dict
      names(M1_mean) <- dict[names(M1_mean), paste0(organism, "_symbol")]
      M1_mean <- M1_mean[!is.na(names(M1_mean))]
      names(M2_mean) <- dict[names(M2_mean), paste0(organism, "_symbol")]
      M2_mean <- M2_mean[!is.na(names(M2_mean))]
      names(M0_mean) <- dict[names(M0_mean), paste0(organism, "_symbol")]
      M0_mean <- M0_mean[!is.na(names(M0_mean))]
    }

    # Spectrum genes subset
    spec_genes <- unique(c(names(M1_mean), names(M2_mean), names(M0_mean)))
    spec_genes <- intersect(spec_genes, rownames(mac_mtx))
    mtx_spec <- as.matrix(mac_mtx[spec_genes, ])
    mtx_spec <- mtx_spec - rowMeans(mtx_spec)

    MPI_genes <- intersect(names(M1_mean), rownames(mtx_spec))
    M1_mean <- M1_mean[MPI_genes]; M2_mean <- M2_mean[MPI_genes]
    mtx_mpi <- mtx_spec[MPI_genes, ]

    AMDI_genes <- intersect(names(M0_mean), rownames(mtx_spec))
    M0_mean <- M0_mean[AMDI_genes]
    mtx_amdi <- mtx_spec[AMDI_genes, ]

    mac_sigma <- sqrt(colSums(mtx_mpi^2) / length(MPI_genes))
    mac_sigma_m0 <- sqrt(colSums(mtx_amdi^2) / length(AMDI_genes))

    M0_sigma <- sqrt(sum(M0_mean^2) / length(M0_mean))
    M1_sigma <- sqrt(sum(M1_mean^2) / length(M1_mean))
    M2_sigma <- sqrt(sum(M2_mean^2) / length(M2_mean))

    M0_pearson <- colSums((mtx_amdi / mac_sigma_m0) * (M0_mean / M0_sigma)) / length(AMDI_genes)
    M1_pearson <- colSums((mtx_mpi / mac_sigma) * (M1_mean / M1_sigma)) / length(MPI_genes)
    M2_pearson <- colSums((mtx_mpi / mac_sigma) * (M2_mean / M2_sigma)) / length(MPI_genes)

    a <- 0.991414467; b <- 1; c <- -0.0185412856
    d_sqr <- (a * M1_pearson + b * M2_pearson + c)^2 / (a^2 + b^2)
    x_start <- -1; y_start <- 1; x_end <- 1; y_end <- -1
    l <- sqrt(pmax(0, (M1_pearson - x_start)^2 + (M2_pearson - y_start)^2 - d_sqr))
    l_max <- sqrt(pmax(0, (x_end - x_start)^2 + (y_end - y_start)^2 - d_sqr))

    return(data.frame(
      MPI = (l / l_max) * 100 - 50,
      AMDI = -M0_pearson * 50,
      MDFI = MFDI_score(mac_mtx, organism),
      row.names = colnames(mac_mtx)
    ))
  }

  spec_results <- calc_ms(mtx, organism)
  object@meta.data <- cbind(object@meta.data, spec_results)

  # Plotting
  dfplot <- object@meta.data
  if (mode == "macspectrum") {
    p <- ggplot(dfplot, aes(x = MPI, y = AMDI)) +
      geom_vline(xintercept = 0, colour = 'grey30', linetype = "dashed") +
      geom_hline(yintercept = 0, colour = 'grey30', linetype = "dashed")
  } else {
    p <- ggplot(dfplot, aes(x = MPI, y = MDFI)) +
      geom_vline(xintercept = 0, colour = 'grey30', linetype = "dashed") +
      geom_hline(yintercept = 0, colour = 'grey30', linetype = "dashed")
  }

  p <- p + theme_classic(base_size = base.size) + labs(color = group.by, fill = group.by)

  if (style == 'point') {
    p <- p + geom_point(aes(color = .data[[group.by]]), size = 1, alpha = 0.6)
  } else {
    p <- p + stat_density_2d(aes(color = .data[[group.by]], fill = .data[[group.by]]),
                             geom = "polygon", alpha = 0.2)
  }

  if (!is.null(split.by)) {
    p <- p + facet_wrap(vars(.data[[split.by]]), ncol = ncol)
  }
  if (!is.null(cols)) {
    p <- p + scale_color_manual(values = cols) + scale_fill_manual(values = cols)
  }

  print(p)
  return(object)
}
