#' @title single-cell analysis of one sample
#' @description A wrapper of the standard workflow of Seurat
#' @param object Seurat object
#' @param nFeatures Number of variable features to use
#' @param nPC Number of principal components to use
#' @param resolution Value of the resolution parameter
#' @param perplexity Value of the perplexity parameter, set between 5 and 50 
#'
#' @return Seurat object
#'
#' @importFrom Seurat AddMetaData Embeddings FindClusters FindNeighbors 
#'   FindVariableFeatures Idents NormalizeData RunPCA RunTSNE RunUMAP 
#'   ScaleData VariableFeatures
#'
#' @references Inspired by Stuart and Butler et al, Cell (2019)
#' @export
#'
#' @examples
#' Bystander <- scRead(sample_name = "Bystander", 
#'                     data_dir = system.file(
#'                       "extdata/H3N2_10X_matrix/Bystander/", 
#'                     package="Yeskit"), 
#'                     gene_column = 2, project_name = "H3N2", 
#'                     group_name = "Bystander", 
#'                     meta_file = system.file(
#'                       "extdata/H3N2_10X_matrix/Bystander/microbes.tsv", 
#'                       package="Yeskit")
#'                     )
#' Infected <- scRead(sample_name = "Infected", 
#'                     data_dir = system.file(
#'                       "extdata/H3N2_10X_matrix/Infected/", 
#'                     package="Yeskit"), 
#'                     gene_column = 2, project_name = "H3N2", 
#'                     group_name = "Infected", 
#'                     meta_file = system.file(
#'                       "extdata/H3N2_10X_matrix/Infected/microbes.tsv", 
#'                       package="Yeskit")
#'                     )
#' Integrated <- merge(Bystander, Infected)
#' Integrated <- scOne(object = Integrated, 
#'   nFeatures = 2000,
#'   nPC = 30,
#'   resolution = 0.7
#' )
#'
scOne <- function(object = NULL, nFeatures = 2000, nPC = 30, resolution = 0.5, 
                  perplexity = NULL) {
  object <- Seurat::NormalizeData(object, 
                                  normalization.method = "LogNormalize",
                                  scale.factor = 10000)
  object <- Seurat::FindVariableFeatures(object, selection.method = "vst", 
    nfeatures = nFeatures)
  object <- Seurat::ScaleData(object, features = rownames(object))
  object <- Seurat::RunPCA(object, features = Seurat::VariableFeatures(object))
  object <- Seurat::FindNeighbors(object, dims = seq_len(nPC))
  object <- Seurat::FindClusters(object, resolution = resolution)
  object <- Seurat::RunUMAP(object, reduction = "pca", dims = seq_len(nPC))
  if (is.null(perplexity)) {
    object <- Seurat::RunTSNE(object, reduction = "pca", dims = seq_len(nPC))
  } else {
    object <- Seurat::RunTSNE(object, reduction = "pca", dims = seq_len(nPC),
      perplexity = perplexity)
  }
  cols <- NA
  if (length(levels(Seurat::Idents(object))) <= 36) {
    cols = c("#1660A7", "#FF6A00", "#219418", "#CD0C18", "#814BB2", "#794339",
             "#DC59B6", "#CC79A7", "#FF0000", "#11B3C6", "#AFB400", "#00FFFF", 
             "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", 
             "#D55E00", "#D55E00", "#CC79A7", "#00AFBB", "#E69F00", "#009E73", 
             "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#4477AA", 
             "#EE6677", "#228833", "#CCBB44", "#66CCEE", "#AA3377", "#BBBBBB")
    object@misc$cols = cols
  }
  reductions <- intersect(c("pca", "tsne", "umap"), names(object))
  for (reduction in reductions) {
    meta_ids <- gsub("coord", toupper(reduction), c("coord_1", "coord_2"))
    coord <- Seurat::Embeddings(object = object, 
                                reduction = reduction)[, c(1, 2)]
    object <- Seurat::AddMetaData(object = object, 
                                  metadata = coord, 
                                  col.name = meta_ids)
  }
  return(object)
}
