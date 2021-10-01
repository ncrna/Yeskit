#' GO annotation for single-cell
#' @title scGO
#' @param object Seurat object
#' @param key Keyword where DGEs are stored
#' @param clusters A vector of clusters, default clusters=NULL
#' @param logFC A positive value to set the cutoff of logFC
#' @param only.pos Only up-regulated genes are taken into account, 
#'   default only.pos=TRUE
#' @param reference Reference organism to set ['human' or 'mouse'], 
#'   default reference = 'human'
#' @param extra Extra slot used to fetch DGEs
#' @return A list object that stored GO results.
#'
#' @importFrom methods new
#' @importFrom org.Hs.eg.db org.Hs.egALIAS2EG 
#' @importFrom org.Mm.eg.db org.Mm.egALIAS2EG
#' @import topGO
#' 
#' @references Inspired by Yuan et al, Methods in Molecular Biology (2019)
#' @export
#'
#' @examples
#' data("H3N2_small")
#' require(topGO)
#' x <- scGO(object = H3N2_small,
#'   key = "Infected_vs_Bystander", 
#'   logFC = 0.25,
#'   only.pos = FALSE,
#'   reference = "human"
#' )
#' head(x)
#'
scGO <- function(object = NULL, key = NULL, clusters = NULL, logFC = 0.25, 
                 only.pos = TRUE, reference = "human", extra = NULL) {
  cluster <- avg_log2FC <- p_val_adj <- NULL
  if (is.null(object)) {
    stop("Parameter 'object' must be specified!\n")
  }
  if (is.null(key)) {
    stop("Parameter 'key' must be specified!\n")
  }
  if (!key %in% names(object@misc)) {
    stop("The 'key' ", key, " does not exist!\n")
  }
  directions <- c("up", "down")
  if (only.pos) {
    directions <- c("up")
  }
  if (is.null(extra)) {
    markers <- object@misc[[key]]
  } else {
    markers <- object@misc[[key]][[extra]]
  }
  if (!is.null(clusters)) {
    markers <- subset(markers, cluster %in% clusters)
  }
  # Fix the header of the DGE result with early version of 'MAST' package
  if ("avg_logFC" %in% colnames(markers)) {
    colnames(markers)[which(colnames(markers) == "avg_logFC")] <- "avg_log2FC"
  }
  x <- list()
  mappings <- ""
  if (reference == "human") {
    x <- as.list(org.Hs.eg.db::org.Hs.egALIAS2EG)
    mappings <- "org.Hs.eg.db"
  } else if (reference == "mouse") {
    x <- as.list(org.Mm.eg.db::org.Mm.egALIAS2EG)
    mappings <- "org.Mm.eg.db"
  } else {
    stop("Parameter 'reference' must be specified with 'human' or 'mouse'!\n")
  }
  geneList <- rep(0, length(rownames(object)))
  names(geneList) <- rownames(object)
  geneList <- geneList[intersect(names(geneList), names(x))]
  TotalGenes <- names(geneList)
  for (ii in seq_len(length(geneList))) {
    names(geneList)[ii] = x[[names(geneList)[ii]]][1]
  }
  go_enrichment_results = list()
  for (c in as.character(unique(markers$cluster))) {
    print(paste0("Running cluster ", c))
    go_enrichment_results[[c]] = list()
    for (direction in directions) {
      if (direction == "up") {
        queryMatrix = subset(markers[order(markers$pct.1/markers$pct.2, 
          decreasing = TRUE), ], cluster == eval(quote(c)) & 
          avg_log2FC > abs(logFC) & p_val_adj < 0.01)
      } else {
        queryMatrix = subset(markers[order(markers$pct.2/markers$pct.1, 
          decreasing = TRUE), ], cluster == eval(quote(c)) & 
          avg_log2FC < -abs(logFC) & p_val_adj < 0.01)
      }
      go_enrichment_results[[c]][[direction]] = list()
      queryGene = rownames(queryMatrix)
      if (length(queryGene) < 1) {
        warning("No ", direction, 
                "-regulated genes left to perform GO analysis!")
        next
      }
      # Run against topGO ####
      queryGeneList = geneList
      queryGeneList[which(TotalGenes %in% queryGene)] = 1
      # skip if no query gene left
      if (all(queryGeneList == 0)) {
        warning("No ", direction, 
                "-regulated genes left to perform GO analysis!")
        next
      }
      tab = list()
      for (ont in c("BP", "CC", "MF")) {
        GOdata <- suppressMessages(methods::new("topGOdata", ontology = ont, 
          allGenes = as.factor(queryGeneList), nodeSize = 10, 
          annot = topGO::annFUN.org, mapping = mappings, ID = "entrez"))
        resultTopGO.elim <- suppressMessages(topGO::runTest(GOdata, 
          algorithm = "elim", statistic = "Fisher"))
        tab[[ont]] <- rbind(tab[[ont]], topGO::GenTable(GOdata, 
          Fisher.elim = resultTopGO.elim, 
          orderBy = "Fisher.elim", topNodes = 200))
      }
      go_enrichment_results[[c]][[direction]] = tab
    }
  }
  results <- list()
  if (!is.null(extra)) {
    results[[extra]] <- go_enrichment_results
  } else {
    results <- go_enrichment_results
  }
  return(results)
}
