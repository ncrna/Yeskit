#' GO annotation for scPathogenDGE results
#' @title scPathogenGO annotation
#' @param object Seurat object
#' @param key The DGEs are stored in object@misc$key
#' @param clusters A vector of clusters, default all clusters will be used
#' @param species Species which is used to divid cells into two groups
#' @param logFC A positive value to set the cutoff of logFC
#' @param only.pos Only up-regulated genes are taken into account, 
#'   default only.pos=TRUE
#' @param reference Reference organism to set ['human' or 'mouse'], 
#'   default reference = 'human'
#' @return Table that stores scPathogenGO results
#'
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
#' x <- scPathogenGO(object = H3N2_small,
#'   key = "H3N2", 
#'   species = "H3N2",
#'   logFC = 0.25,
#'   only.pos = FALSE,
#'   reference = "human"
#' )
#' head(x)
#'
scPathogenGO <- function(object = NULL, key = NULL, species = NULL, 
                         clusters = NULL, logFC = 0.25, only.pos = TRUE, 
                         reference = "human") {
  avg_log2FC <- p_val_adj <- Data <- cluster <- NULL
  if (is.null(object)) {
    stop("Parameter 'object' must be specified!\n")
  }
  if (is.null(key)) {
    stop("Parameter 'key' must be specified!\n")
  }
  if (!key %in% names(object@misc)) {
    stop("The 'key' ", key, " does not exist!\n")
  }
  if (is.null(species)) {
    stop("Parameter 'species' must be specified!\n")
  }
  if (!species %in% names(object@misc[[key]])) {
    stop("The 'species' ", species, " does not exist!\n")
  }
  if (is.null(clusters)) {
    warning("All clusters will be evaluated!\n")
    clusters <- levels(object)
  }
  clusters <- as.character(clusters)
  direction <- c("up", "down")
  if (only.pos) {
    directions <- c("up")
  }
  x <- list()
  mappings <- NULL
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
  go_enrichment_results[[species]] = list()
  Data = object@misc[[key]][[species]]
  # Fix the header of the DGE result with early version of 'MAST' package
  if ("avg_logFC" %in% colnames(Data)) {
    colnames(Data)[which(colnames(Data) == "avg_logFC")] <- "avg_log2FC"
  }
  for (C in clusters) {
    DEdata <- subset(Data, cluster == C)
    if (nrow(DEdata) == 0) {
      warning("No DE genes between species ", species, 
              " pos/neg in cluster-", C)
      next
    }
    go_enrichment_results[[species]][[C]] = list()

    message("Running ", species, " in cluster-", C, " ... ")

    for (direction in c("up", "down")) {
      if (direction == "up") {
        queryMatrix = subset(DEdata[order(DEdata$pct.1/DEdata$pct.2, 
          decreasing = TRUE), ], avg_log2FC > abs(logFC) & p_val_adj < 0.01)
      } else {
        queryMatrix = subset(DEdata[order(DEdata$pct.2/DEdata$pct.1, 
          decreasing = TRUE), ], avg_log2FC < -abs(logFC) & p_val_adj < 0.01)
      }
      queryGene = rownames(queryMatrix)
      # next if queryGene with only one Gene
      if (length(queryGene) < 2) {
        next
      }
      # Run against topGO ####
      queryGeneList = geneList
      queryGeneList[which(TotalGenes %in% queryGene)] = 1
      if (length(table(queryGeneList)) < 2) {
        next
      }
      queryMatrix$gene = rownames(queryMatrix)
      queryGene = queryMatrix$gene
      go_enrichment_results[[species]][[C]][[direction]] = list()

      tab = list()
      for (ont in c("BP", "CC", "MF")) {
        GOdata <- suppressMessages(new("topGOdata", ontology = ont, 
          allGenes = as.factor(queryGeneList),
          nodeSize = 10, annot = topGO::annFUN.org, mapping = mappings, 
          ID = "entrez"))
        resultTopGO.elim <- suppressMessages(topGO::runTest(GOdata, 
          algorithm = "elim", statistic = "Fisher"))
        tab[[ont]] <- rbind(tab[[ont]], topGO::GenTable(GOdata, 
          Fisher.elim = resultTopGO.elim,
          orderBy = "Fisher.elim", topNodes = 200))
      }
      go_enrichment_results[[species]][[C]][[direction]] = tab
    }
    message("done.\n")
  }
  return(go_enrichment_results)
}
