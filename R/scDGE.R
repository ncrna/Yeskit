#' Differential gene expression analysis for each cluster
#' Perform differential analysis (DE) between two groups for each clusters
#' @param object Seurat object.
#' @param comparison Vector of nominal variables from group.by. 
#'   Eg., comparison=c('Normal', 'Tumor')
#' @param group.by Regroup cells before performing DGE, 
#'   default group.by='group'
#' @param min.cells Minimum number of cells for DGE in any of the two groups, 
#'   default min.cells=20
#' @param min.pct Genes expressed in any of the two groups with a minium cell
#'   fraction (min.pct) are reserved for DGE, default min.pct=0.1
#' @param logFC A positive value to set the cutoff of logFC, default logFC=0.25
#' @param clusters Vector of clusters to perform DGE, default clusters=NULL,
#'   which means all clusters will be evaluated
#' @seealso \code{\link[Seurat]{FindMarkers}}
#'
#' @return DGE Table
#'
#' @import Seurat
#' @importFrom Seurat FindMarkers
#'
#' @export
#'
#' @examples
#' data("H3N2_small")
#' x <- scDGE(object = H3N2_small, 
#'   comparison = c("Infected", "Bystander"),
#'   group.by = "group", 
#'   min.cells = 10, 
#'   logFC = 0.25, 
#'   clusters = NULL
#' )
#' head(x)
#'
scDGE <- function(object = object, comparison = c("condA", "condB"), 
                  group.by = NULL, min.cells = 20, min.pct = 0.1,
                  logFC = 0.25, clusters = NULL) {
  results <- list()
  if (is.null(clusters)) clusters = levels(object)
  if (length(comparison) != 2) stop("Comparison must have 2 elements!")
  condA <- comparison[1]
  condB <- comparison[2]
  if (is.null(group.by)) group.by <- "group"
  clusters <- as.character(clusters)
  DE <- data.frame()
  for (cluster in clusters) {
    message("### ", "Comparing cluster-", cluster, " between ", 
            condA, " and ", condB, " ...\n")
    Object <- subset(object, clusters == cluster)
    if (all(comparison %in% names(table(Object[[group.by]])))) {
      Object <- Object[, Object@meta.data[[group.by]] %in% comparison]
    } else {
      warning("Cluster ", cluster, " has no cell in ", 
              condA, " or ", condB, ". Ignored.\n")
      next
    }
    if (!all(as.data.frame(table(Object[[group.by]]))$Freq >= min.cells)) {
      warning("Cell number in cluster ", cluster, 
              " is less than ", min.cells, " cells. Ignored.\n")
      next
    } else if (!condA %in% names(table(Object$group))) {
      warning("condA has no cell in cluster ", cluster, ". Ignored.\n")
      next
    } else if (!condB %in% names(table(Object$group))) {
      warning("condB has no cell in cluster ", cluster, ". Ignored.\n")
      next
    }
    tmpDE <- suppressWarnings(FindMarkers(object = Object, 
      assay = "RNA", ident.1 = condA, ident.2 = condB, 
      group.by = group.by, min.cells.group = min.cells, min.pct = min.pct,
      logfc.threshold = logFC, test.use = "MAST", only.pos = FALSE)
    )
    tmpDE$cluster <- cluster
    tmpDE$gene <- rownames(tmpDE)
    DE <- rbind(DE, tmpDE)
  }
  message("done.")
  return(DE)
}

