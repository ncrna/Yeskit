#' read scRNA matrix
#' @title scread function
#' @description A function to read single-cell expression matrix
#' @details It supports read gene expression matrix and antibody
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
#'   package="Yeskit"),
#' )    
#' x
#'
scread <- function(sample_name = NULL, data_dir = NULL, gene_column = 2, 
                   min_cells = 5, min_features = 200, 
                   min_rnas = 1000, percent_mito = 20, 
                   mito_pattern = "^GRCh38_MT-|^mm10___mt-|^MT-|^mt-",
                   project_name = "Seurat", group_name = NULL, 
                   strip_suffix = TRUE) {
  object <- Seurat::Read10X(data.dir = data_dir, gene.column = gene_column, 
			    strip.suffix = strip_suffix)
  object.rna <- object$`Gene Expression`
  object.adt <- object$`Antibody Capture`
  object <- Seurat::CreateSeuratObject(counts = object.rna, 
				       min.cells = min_cells, min.features = min_features, 
				       project = project_name)
  object$sample <- sample_name
  object$group <- group_name
  mito_pattern <- gsub("_", "-", mito_pattern)
  object[["percent.mito"]] <- Seurat::PercentageFeatureSet(object, 
							   pattern = mito_pattern)
  p1 <- Seurat::VlnPlot(object = object, features = c("nFeature_RNA", 
						      "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.01)
  # QC and selecting cells for further analysis
  nFeature_RNA <- nCount_RNA <- percent.mito <- NULL
  object <- subset(object, subset = nFeature_RNA > max(min_features, 
    stats::quantile(object$nFeature_RNA, 0.01)) & 
    nFeature_RNA < stats::quantile(object$nFeature_RNA, 0.99) & percent.mito <=
    min(percent_mito, stats::quantile(object$percent.mito, 0.99)) & 
    nCount_RNA > max(min_rnas, stats::quantile(object$nCount_RNA, 0.01)) & 
    nCount_RNA < stats::quantile(object$nCount_RNA, 0.99))
  object[['ADT']] <- Seurat::CreateAssayObject(counts = object.adt[, colnames(x = object)])
  p2 <- Seurat::VlnPlot(object = object, features = c("nFeature_RNA", 
						      "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.01)
  plot(patchwork::wrap_plots(list(p1, p2), ncol = 1))
  return(object)
}
