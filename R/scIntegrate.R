# Set a default value if the left value is null 
# @param lhs Default value to use
# @param rhs The value to use if lhs is null 
# @author Hadley Wickham 
# @references https://adv-r.hadley.nz/functions.html#missing-arguments
`%||%` <- NULL
`%||%` <- function(lhs, rhs) {
    if (!is.null(lhs)) {
        lhs
    } else {
        rhs
    }
}

#' Integrate seurat objects into one
#' @title scIntegrate
#' @param object.list A list of seurat objects
#' @param object.names An array of seurat object names
#' @param nVariable Number of features to select as top variable features
#' @param nPC Number of PCs to use
#' @param resolution Resolution parameter to set. Default resolution=0.5
#' @param batch.rm Remove batch effect with 'seurat' or 'harmony'. 
#'   Default batch.rm='harmony'
#' @return Seurat object.
#'
#' @importFrom Seurat AddMetaData DefaultAssay Embeddings FindClusters 
#'   FindNeighbors FindVariableFeatures Idents NormalizeData RunPCA 
#'   RunTSNE RunUMAP ScaleData VariableFeatures
#' @importFrom dplyr bind_rows
#' @importFrom harmony RunHarmony
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
#' Integrated <- scIntegrate(object.list = list(Bystander, Infected), 
#'   object.names = c("Bystander", "Infected"),
#'   batch.rm = "harmony",
#'   resolution = 0.7
#' )
#'
scIntegrate <- function(object.list = NULL, object.names = NULL, 
                        nVariable = 2000, nPC = 30, resolution = 0.5, 
                        batch.rm = "harmony") {
  # Merge each seurat object
  meta.list <- list()
  meta.table <- object.list[[1]]@meta.data
  for (i in 2:length(object.list)) {
    meta.table <- dplyr::bind_rows(meta.table, object.list[[i]]@meta.data)
  }
  meta.table[is.na(meta.table)] <- 0

  object <- merge(object.list[[1]], y = object.list[2:length(object.list)], 
    add.cell.ids = as.character(object.names), project = "seurat")
  rownames(meta.table) <- rownames(object@meta.data)
  object@meta.data <- meta.table

  object <- Seurat::NormalizeData(object, verbose = TRUE)
  object <- Seurat::FindVariableFeatures(object, selection.method = "vst", 
    nfeatures = nVariable)
  object <- Seurat::ScaleData(object, vars.to.regress = c("nCount_RNA", 
    "percent.mito"), verbose = TRUE)
  object <- Seurat::RunPCA(object, pc.genes = Seurat::VariableFeatures(object),
    npcs = nPC, verbose = TRUE)

  object <- harmony::RunHarmony(object = object, group.by.vars = "sample", 
    assay.use = Seurat::DefaultAssay(object), plot_convergence = TRUE)

  if (batch.rm == "harmony") {
    object <- Seurat::FindNeighbors(object, reduction = "harmony", 
      dims = seq_len(nPC))
    object <- Seurat::FindClusters(object, resolution = resolution)
    object <- Seurat::RunUMAP(object, reduction = "harmony", 
      dims = seq_len(nPC))
    #object <- Seurat::RunTSNE(object, reduction = "harmony", 
    #  dims = seq_len(nPC), do.fast = TRUE)
  } else if (batch.rm == "seurat") {
    object <- Seurat::FindNeighbors(object, reduction = "pca", 
      dims = seq_len(nPC))
    object <- Seurat::FindClusters(object, resolution = resolution)
    object <- Seurat::RunUMAP(object, reduction = "pca", dims = seq_len(nPC))
    #object <- Seurat::RunTSNE(object, reduction = "pca", dims = seq_len(nPC))
  } else {
    stop("batch.rm must be 'harmony' or 'seurat'!")
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
  # Add reduction coordinates into Meta.data slot
  reductions <- intersect(c("pca", "tsne", "umap"), names(object))
  for (reduction in reductions) {
    meta_ids <- gsub("coord", toupper(reduction), c("coord_1", "coord_2"))
    coord <- Seurat::Embeddings(object = object, 
                                reduction = reduction)[, c(1, 2)]
    object <- Seurat::AddMetaData(object = object, metadata = coord, 
      col.name = meta_ids)
  }
  # Add cluster information into Meta.data slot
  object[["clusters"]] <- Seurat::Idents(object)
  object <- JoinLayers(object)
  return(object)
}

