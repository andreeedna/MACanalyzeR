PathCreatoR <- function(organism, subcollection = NULL, pattern = NULL, min.size = 1, max.size = Inf) {
  if (missing(organism) || !organism %in% c("mm", "hs")) {
    stop("The organism argument must be 'mm' (Mus musculus) or 'hs' (Homo sapiens)", call. = FALSE)
  }
  collections <- if (organism == "hs") MACanalyzeR::gsea_human else MACanalyzeR::gsea_mouse

  if (is.null(subcollection)) {
    message("Available ", if (organism == "hs") "human" else "mouse", " collections (number of gene sets):")
    n <- lengths(collections)
    message(paste0("  ", names(n), " (", n, ")", collapse = "\n"))
    return(invisible(names(collections)))
  }

  ## case-insensitive match; "GO:BP" also matches "M5:GO:BP" / "C5:GO:BP"
  key <- toupper(subcollection)
  ## accept the human codes for mouse (C5:GO:BP -> M5:GO:BP, H -> MH) and vice versa
  if (organism == "mm") key <- sub("^H$", "MH", sub("^C([0-9])", "M\\1", key))
  else key <- sub("^MH$", "H", sub("^M([0-9])", "C\\1", key))
  nm <- names(collections)
  hit <- nm[toupper(nm) == key]
  if (length(hit) == 0) hit <- nm[endsWith(toupper(nm), paste0(":", key))]
  if (length(hit) != 1) {
    stop("Subcollection '", subcollection, "' ", if (length(hit) > 1) "is ambiguous" else "not found",
         ". Available ", if (organism == "hs") "human" else "mouse", " collections: ",
         paste(nm, collapse = ", "), call. = FALSE)
  }

  gs <- collections[[hit]]
  if (!is.null(pattern)) gs <- gs[grepl(pattern, names(gs), ignore.case = TRUE)]
  size <- lengths(gs)
  gs <- gs[size >= min.size & size <= max.size]

  if (length(gs) == 0) warning("No gene set left after filtering", call. = FALSE)
  message("Collection ", hit, ": ", length(gs), " gene sets")
  gs
}
