#' read scRNA matrix
#' @title scRead function
#' @description A function to read single-cell expression matrix
#' @details It supports read gene expression matrix and microbiome
#' quantification matrix simultaneously
#' @param sample_name Name of the sample
#' @param data_dir Directory containing the quantification matrix
#' (matrix.mtx, features.tsv, and barcodes.tsv)
#' @param gene_column Column of genes.tsv to use for gene names,
#'   default is 2
#' @param min_cells Minimum number of cells express this feature
#' @param min_features Minimum number of features expressed in this cell
#' @param min_rnas Minimum number of molecules detected within a cell.
#' @param percent_mito Maximum percent of mito in a cell
#' @param mito_pattern regex pattern for mitochondrial genes, 
#'   default is '^GRCh38_MT-|^mm10___mt-|^MT-|^mt-'
#' @param project_name Project name for the Seurat object
#' @param group_name Group name for the Seurat object
#' @param strip_suffix Remove trailing '-1' in cell barcodes
#' @param meta_file Meta file to use
#' @param human.prefix The prefix denoting rownames for human cells. 
#'   Default is 'GRCh38_'. Just for PDX samples
#' @param mouse.prefix The prefix denoting rownames for mouse cells. 
#'   Default is 'mm10___'. Just for PDX samples
#' @param organism.thres Threshould fraction of reads to separate organism 
#'   and background, Default is 0.9. Just for PDX samples
#' @param organism.use Which organism is used. Default is NULL, 
#'   means all organisms will be used. Just for PDX samples
#' @param meta_features Features used to store in the meta.data slot
#' @return A Seurat object.
#'
#' @importFrom Matrix colSums
#' @importFrom Seurat AddMetaData CreateSeuratObject PercentageFeatureSet 
#'   Read10X VlnPlot
#' @importFrom patchwork wrap_plots
#' @importFrom stats quantile
#' @importFrom utils read.table
#'
#' @references Inspired by Stuart and Butler et al, Cell (2019)
#' @export
#'
#' @examples
#' x <- scRead(sample_name = "Infected", 
#'   data_dir = system.file("extdata/H3N2_10X_matrix/Infected/",
#'                          package="Yeskit"), 
#'   gene_column = 2, 
#'   project_name = "H3N2", 
#'   group_name = "Infected",
#'   meta_file = system.file("extdata/H3N2_10X_matrix/Infected/microbes.tsv",
#'                           package="Yeskit"),
#' )    
#' x
#'
scRead <- function(sample_name = NULL, data_dir = NULL, gene_column = 2, 
                   min_cells = 5, min_features = 200, 
                   min_rnas = 1000, percent_mito = 20, 
                   mito_pattern = "^GRCh38_MT-|^mm10___mt-|^MT-|^mt-",
                   project_name = "Seurat", group_name = NULL, 
                   strip_suffix = TRUE, meta_file = NULL,
                   human.prefix = "GRCh38_", mouse.prefix = "mm10___", 
                   organism.thres = 0.9, organism.use = NULL,
                   meta_features = NULL) {
  object <- Seurat::Read10X(data.dir = data_dir, gene.column = gene_column, 
    strip.suffix = strip_suffix)
  if ('Gene Expression' %in% names(object)) {
    object <- object$`Gene Expression`
  }
  hsaFeatures <- grep(pattern = human.prefix, x = rownames(x = object), 
    value = TRUE)
  mmuFeatures <- grep(pattern = mouse.prefix, x = rownames(x = object), 
    value = TRUE)
  organism.info <- data.frame()
  if (length(hsaFeatures) > 0 | length(mmuFeatures) > 0) {
    hsa.object <- object[hsaFeatures, ]
    mmu.object <- object[mmuFeatures, ]
    nUMIs <- Matrix::colSums(object)
    hsaUMIs <- Matrix::colSums(hsa.object)
    mmuUMIs <- Matrix::colSums(mmu.object)
    hsa.object <- hsa.object[, hsaUMIs/nUMIs >= organism.thres]
    mmu.object <- mmu.object[, mmuUMIs/nUMIs >= organism.thres]
    organism.info <- data.frame(organism = c(rep("human", 
      length(colnames(hsa.object))),
      rep("mouse", length(colnames(mmu.object)))), 
      row.names = c(colnames(hsa.object),
      colnames(mmu.object)))
    if (is.null(organism.use)) {
      indcs <- which(colnames(object) %in% c(colnames(hsa.object), 
        colnames(mmu.object)))
      object <- object[, indcs, drop = FALSE]
    } else if (organism.use == "human") {
      rownames(hsa.object) <- gsub(human.prefix, "", rownames(hsa.object))
      object <- hsa.object
    } else if (organism.use == "mouse") {
      rownames(mmu.object) <- gsub(mouse.prefix, "", rownames(mmu.object))
      object <- mmu.object
    }
  }
  metadata <- data.frame()
  if (!is.null(meta_features)) {
    meta_features <- intersect(meta_features, rownames(object))
  }
  if (!is.null(meta_features)) {
    metadata <- data.frame(t(as.matrix(object[meta_features, ])))
    object <- Seurat::CreateSeuratObject(counts = object[!rownames(object) %in%
      meta_features, ], min.cells = min_cells, min.features = min_features,
      project = project_name, meta.data = metadata)
  } else {
    object <- Seurat::CreateSeuratObject(counts = object, 
      min.cells = min_cells, min.features = min_features, 
      project = project_name)
  }
  object$sample <- sample_name
  object$group <- group_name
  object$organism <- NA
  if (length(organism.info) > 0) {
    organism_df <- merge(data.frame(object$organism, 
      row.names = colnames(object)),
      organism.info, by = 0, order = FALSE)
    object$organism <- organism_df[, "organism"]
  }
  mito_pattern <- gsub("_", "-", mito_pattern)
  object[["percent.mito"]] <- Seurat::PercentageFeatureSet(object, 
    pattern = mito_pattern)
  p1 <- Seurat::VlnPlot(object = object, features = c("nFeature_RNA", 
    "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.01)

  if (!is.null(meta_file)) {
    if (file.exists(meta_file)) {
      metadata <- utils::read.table(meta_file, header = TRUE, row.names = 1, 
        stringsAsFactors = FALSE, sep = "\t", quote = "", comment.char = "")
      metadata[is.na(metadata)] <- 0
      rownames(metadata) <- gsub("'", "", rownames(metadata))
      rownames(metadata) <- gsub("\\[", "", rownames(metadata))
      rownames(metadata) <- gsub("\\]", "", rownames(metadata))
      rownames(metadata) <- gsub(" ", "_", rownames(metadata))
      metadata_ids <- rownames(metadata)
      if (nrow(metadata) > 0) {
        object <- Seurat::AddMetaData(object = object, metadata = t(metadata),
          col.name = metadata_ids)
      }
    } else {
      stop("Meta file does not exist! Please check the meta_file parameter")
    }
  }
  # QC and selecting cells for further analysis
  nFeature_RNA <- nCount_RNA <- percent.mito <- NULL
  object <- subset(object, subset = nFeature_RNA > max(min_features, 
    stats::quantile(object$nFeature_RNA, 0.01)) & 
    nFeature_RNA < stats::quantile(object$nFeature_RNA, 0.99) & percent.mito <=
    min(percent_mito, stats::quantile(object$percent.mito, 0.99)) & 
    nCount_RNA > max(min_rnas, stats::quantile(object$nCount_RNA, 0.01)) & 
    nCount_RNA < stats::quantile(object$nCount_RNA, 0.99))
  p2 <- Seurat::VlnPlot(object = object, features = c("nFeature_RNA", 
    "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.01)
  plot(patchwork::wrap_plots(list(p1, p2), ncol = 1))
  return(object)
}
