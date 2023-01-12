#' MSigdb signature scoring for single-cells
#' @title scMsigdbScoring
#' @param object Seurat object
#' @param species Species name. Such as "human" or "mouse"
#' @param category String of MSigDB category to use, 
#'   default category='H' stands for 'hallmark gene sets'
#' @param genesets Vector of genesets under the specific category of 
#'   MSIGDB, such as "APOPTOSIS". Default genesets=NULL, which means 
#'   all genesets under the specified category will be used 
#' @return Seurat object
#'
#' @importFrom Seurat AddMetaData AddModuleScore Embeddings
#' @importFrom msigdbr msigdbr
#' @importFrom utils data globalVariables
#' @importFrom stats aggregate
#'
#' @export
#'
#' @examples
#' data("H3N2_small")
#' x <- scMsigdbScoring(object = H3N2_small, 
#'   species = "human", category = "H",
#'   genesets = c("HALLMARK_INFLAMMATORY_RESPONSE")
#' )
#' head(x)
#'
scMsigdbScoring <- function(object = NULL, species = "human", 
                            category = NULL, genesets = NULL) {
  if (is.null(species)){
    species <- "human"
  }
  if (is.null(category)){
    category <- "H"
  }
  MSIGDB <- NULL   
  if (category=="H" & (species == "human" | species == "Homo sapiens")){
    utils::data("HALLMARK")
    MSIGDB <- HALLMARK
  }else{
    all_gene_sets <- msigdbr::msigdbr(species=species, category=category)
    all_gene_sets <- as.data.frame(all_gene_sets[, 
                                   c("gs_cat", "gs_name", "gene_symbol")])
    all_gene_sets <- aggregate(. ~ gs_cat + gs_name, 
                               data = all_gene_sets, FUN = "paste")
    MSIGDB <- list()
    for (i in seq_len(nrow(all_gene_sets))){
      MSIGDB[[all_gene_sets$gs_name[i]]] <- 
        unlist(all_gene_sets$gene_symbol[i])
    }
  }

  if (length(names(MSIGDB)) == 0) {
    stop("MSigDB has no valid entries!")
  }

  if (is.null(genesets)) {
    genesets = names(MSIGDB)
  } else {
    genesets <- intersect(genesets, names(MSIGDB))
    if (length(genesets) == 0) {
      stop("genesets ", genesets, "does not exist in MSIGDB!")
    }
  }
  for (item in genesets) {
    if (item %in% names(object[[]])) {
      object[[item]] <- NULL
    }
    features <- MSIGDB[[item]]
    features <- list(Score = features)
    object.hallmark <- Seurat::AddModuleScore(object = object, 
      features = features, name = item, ctrl = 
        min(vapply(X = features, FUN = length, 
          FUN.VALUE = numeric(length = 1))))
    hallmark.columns <- grep(pattern = item, 
                             x = colnames(x = object.hallmark[[]]), 
                             value = TRUE)
    hallmark.scores <- object.hallmark[[hallmark.columns]]
    rm(object.hallmark)
    colnames(x = hallmark.scores) <- c(item)
    object[[colnames(x = hallmark.scores)]] <- hallmark.scores
  }
  reductions <- intersect(c("pca", "tsne", "umap"), names(object))
  for (reduction in reductions) {
    meta_ids <- gsub("coord", toupper(reduction), c("coord_1", "coord_2"))
    coord <- Seurat::Embeddings(object = object, reduction = reduction)[, 
                                c(1, 2)]
    object <- Seurat::AddMetaData(object = object, 
                                  metadata = coord, 
                                  col.name = meta_ids)
  }
  return(object)
}

utils::globalVariables(name = "HALLMARK")
