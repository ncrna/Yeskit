#' Volcano Plot for DGEs
#' @title scVolcanoPlot
#' @param object A seurat object.
#' @param key The DE matrix are stored in object@misc$key
#' @param cluster Cluster specified by the user
#' @param padj Cutoff of adjusted pvalue
#' @param logFC Log2FoldChange in DGEs
#' @param cols Colours to pain the dots, 
#'   default cols=c('#808080', '#810F7C', '#006D2C', '#000000')
#' @param pt.size Size of the dots, default pt.size=1
#' @param species.by Species which is used to divid cells into two groups
#' @param top_n Number of up- and down-regulated genes to show
#' @param labels An array of gene names to show
#' @param label.size Font size of the labels, default label.size=2.5
#' @param box Draws a rectangle behind the labels, default is TRUE
#' @param raster Convert points to raster format, default is TRUE
#' @return A ggplot object
#'
#' @importFrom ggplot2 aes geom_hline geom_point geom_vline ggplot 
#'    labs scale_color_manual scale_fill_manual unit xlim ylim
#' @importFrom cowplot theme_cowplot
#' @importFrom ggrastr rasterise
#' @importFrom ggrepel geom_label_repel geom_text_repel
#' @importFrom utils head
#'
#' @export
#'
#' @examples
#' data("H3N2_small")
#' scVolcanoPlot(object = H3N2_small,
#'   key = "Infected_vs_Bystander",
#'   cluster = "0",
#'   top_n = 5
#' )
#'
scVolcanoPlot <- function(object = NULL, key = NULL, cluster = NULL, 
                          padj = 0.05, logFC = 0.585, cols = NULL, 
                          pt.size = NULL, species.by = NULL, top_n = NULL,
                          labels = NULL, label.size = NULL, box = TRUE, 
                          raster = TRUE) {
  avg_log2FC <- color <- group <- label <- p_val_adj <- NULL
  if (is.null(object)) {
    stop("Parameter 'object' must be specified!\n")
  }
  if (is.null(key)) {
    stop("Parameter 'key' must be specified!\n")
  }
  if (is.null(cluster)) {
    stop("Parameter 'cluster' must be specified!\n")
  }
  cluster = as.character(cluster)
  if (!key %in% names(object@misc)) {
    stop("The 'key' ", key, " does not exist!\n")
  }
  if (is.null(top_n)) {
    top_n <- 5
  }
  if (is.null(cols)) {
    cols <- c("#808080", "#810F7C", "#006D2C", "#000000")
  }
  if (is.null(pt.size)) {
    pt.size <- 1
  }
  if (is.null(label.size)) {
    label.size <- 2.5
  }
  Data <- data.frame()
  if (is.null(species.by)) {
    Data <- object@misc[[key]]
  } else {
    Data <- object@misc[[key]][[species.by]]
  }
  # Fix the header of the DGE result with early version of 'MAST' package
  if ("avg_logFC" %in% colnames(Data)) {
    colnames(Data)[which(colnames(Data) == "avg_logFC")] <- "avg_log2FC"
  }
  cluster <- as.character(cluster)
  if (!cluster %in% unique(Data$cluster)) {
    stop("cluster ", cluster, " does not exist")
  }
  Data <- Data[Data$cluster == cluster, ]
  xrange <- range(Data$avg_log2FC)
  xrange <- xrange[is.finite(xrange)]
  xmin <- floor(min(xrange))
  xmax <- ceiling(max(xrange))
  xmin <- -max(abs(xmin), abs(xmax))
  xmax <- max(abs(xmin), abs(xmax))
  yrange <- -log10(Data$p_val_adj)
  yrange <- yrange[is.finite(yrange)]
  ymin <- 0
  ymax <- ceiling(max(yrange)) * 1.1
  Data$group = "no"
  if (!any(Data$p_val_adj < padj & Data$avg_log2FC > logFC) & 
    any(Data$p_val_adj < padj & Data$avg_log2FC < -logFC)) {
    stop("No up- or down-regulated genes exist in cluster-", cluster)
  }
  Data[Data$p_val_adj < padj & Data$avg_log2FC > logFC, ]$group <- "up"
  Data[Data$p_val_adj < padj & Data$avg_log2FC < -logFC, ]$group <- "down"
  index <- c()
  if (is.null(labels)) {
    labels <- c(utils::head(Data[Data$group == "up", ]$gene, top_n), 
      utils::head(Data[Data$group == "down", ]$gene, top_n))
  }
  index <- match(labels, Data$gene)
  index <- index[which(!is.na(index))]
  Data$color <- Data$group
  Data$color[index] <- "black"
  names(cols) = c("no", "up", "down", "black")
  Data$label <- ""
  Data$label[index] <- Data$gene[index]
  Data[Data$group == "no", ]$label <- ""
  p <- ggplot2::ggplot(Data, ggplot2::aes(x = avg_log2FC, 
                                          y = -log10(p_val_adj)))
  if (isTRUE(raster)) {
    p <- p + ggrastr::rasterise(dpi = 300, 
      ggplot2::geom_point(ggplot2::aes(fill = group), size = pt.size, 
        shape = 21, alpha = 0.6, show.legend = FALSE))
    p <- p + ggrastr::rasterise(dpi = 300, 
      ggplot2::geom_point(ggplot2::aes(colour = color), size = pt.size, 
        shape = 21, alpha = 0.6, show.legend = FALSE))
  } else {
    p <- p + ggplot2::geom_point(ggplot2::aes(fill = group), size = pt.size,
      shape = 21, alpha = 0.6, show.legend = FALSE)
    p <- p + ggplot2::geom_point(ggplot2::aes(colour = color), size = pt.size,
      shape = 21, alpha = 0.6, show.legend = FALSE)
  }
  p <- p + ggplot2::xlim(c(xmin, xmax)) + ggplot2::ylim(c(ymin, ymax))
  p <- p + ggplot2::scale_color_manual(values = cols)
  p <- p + ggplot2::scale_fill_manual(values = cols)
  p <- p + ggplot2::geom_hline(yintercept = -log10(padj), linetype = "dotted")
  p <- p + ggplot2::geom_vline(xintercept = c(-logFC, logFC),
    linetype = "dotted") + 
    if (isTRUE(box)) {
      p <- ggrepel::geom_label_repel(ggplot2::aes(label = label, color = group),
        show.legend = FALSE, fontface = "bold", size = label.size, 
        box.padding = ggplot2::unit(0.8, "lines"), 
        point.padding = ggplot2::unit(0.3, "lines"), segment.size = 0.3)
    } else {
      p <- ggrepel::geom_text_repel(ggplot2::aes(label = label, color = group),
        show.legend = FALSE, fontface = "bold", size = label.size, 
        box.padding = ggplot2::unit(0.8, "lines"), 
        point.padding = ggplot2::unit(0.3, "lines"), segment.size = 0.3)
    }
  p <- p + ggplot2::labs(x = "Log2 Fold change", y = "Log10 P-adjust", 
    title = paste0(key, " in cluster-", cluster)) + cowplot::theme_cowplot()
  return(p)
}
