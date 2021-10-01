#' Find differentially expressed genes by species for each clusters
#' @param object Seurat object
#' @param species.by String of species to divid cells into 
#'   two groups (Positive / Negative)
#' @param clusters Clusters specified for DGE, default clusters=NULL, means
#'   all clusters will be evaluated
#' @param min.cells Minimum number of cells for DGE in any of the two groups,
#'   default min.cells=20
#' @param min.pct Genes expressed in any of the two groups with a minium cell 
#'   fraction (min.pct) are reserved for DGE, default min.pct=0.1
#' @param logFC A positive value to set the cutoff of logFC, default logFC=0.25
#' @return DGE table
#'
#' @importFrom Seurat FindMarkers GetAssayData
#'
#' @references Inspired by Stuart and Butler et al, Cell (2019)
#'
#' @export
#'
#' @examples
#' data("H3N2_small")
#' x <- scPathogenDGE(object = H3N2_small, 
#'   species.by = "H3N2",
#'   min.cells = 5,
#'   min.pct = 0.1,
#'   logFC = 0.25,
#' )
#' head(x)
#'
scPathogenDGE <- function(object = NULL, species.by = NULL, 
                          clusters = NULL, min.cells = 20, 
                          min.pct = 0.1, logFC = 0.25) {
  if (is.null(object)) {
    stop("Parameter 'object' must be specified!\n")
  }
  if (is.null(clusters)) {
    warning("All clusters will be evaluated!\n")
    clusters = levels(object)
  }
  clusters = as.character(clusters)
  if (!species.by %in% names(object@meta.data)) {
    stop("The feature ", species.by, " does not exist in MetaData slot!\n")
  }
  seurat_clusters <- p_val_adj <- NULL
  results <- list()
  for (feature in species.by) {
    message("### ", "Analysis feature ", feature, " ...\n")
    if (is.null(results[[feature]])) {
      results[[feature]] = list()
    }
    if (!feature %in% colnames(object@meta.data)) {
      warning("all cells with 0 reads.\n")
      next
    }
    DE <- data.frame()
    for (cluster in clusters) {
      message("====== ", "Analysis cluster-", cluster, " ... ")
      Object <- subset(object, seurat_clusters == cluster)
      Object$group <- ifelse(Object@meta.data[, feature] > 0, "Pos", "Neg")
      Expr <- as.matrix(Seurat::GetAssayData(Object))
      if (length(Object$group[Object$group == "Pos"]) < min.cells) {
        warning("too few cells to process.\n")
        next
      }
      if (length(Object$group[Object$group == "Neg"]) < min.cells) {
        warning("too few cells to process.\n")
        next
      }
      tmpDE <- suppressMessages(
        suppressWarnings(
          Seurat::FindMarkers(object = Object, assay = "RNA", 
                              ident.1 = "Pos", ident.2 = "Neg", 
                              group.by = "group",
                              min.cells.group = min.cells, 
                              min.pct = min.pct, 
                              logfc.threshold = logFC, 
                              test.use = "MAST", 
                              only.pos = FALSE, 
                              verbose = FALSE
          )
        )
      )
      tmpDE$cluster <- cluster
      tmpDE$gene <- rownames(tmpDE)
      DE <- rbind(DE, tmpDE)
      message("done.\n")
    }
    results[[feature]] <- DE
  }
  # clean up results
  for (feature in species.by) {
    if (length(results[[feature]]) == 0) {
      results[[feature]] <- NULL
      next
    }
    for (cluster in clusters) {
      if (length(results[[feature]][[cluster]]) == 0) {
        next
      }
      results[[feature]][[cluster]] <- subset(results[[feature]][[cluster]],
        p_val_adj <= 1)
      if (nrow(results[[feature]][[cluster]]) <= 1) {
        results[[feature]][[cluster]] <- NULL
        next
      }
    }
    if (length(results[[feature]]) == 0) {
      results[[feature]] <- NULL
      next
    }
  }
  return(results)
}
