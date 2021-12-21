#' DimPlot with rasterized point for single-cell visualization
#' @title scDimPlot
#' @param object Seurat object
#' @param cols Colors to plot, use ggplot2's default colors by default. 
#'   We include a pallete called 'sc' which consists of 36 colors
#' @param pt.size Adjust point size to plot, default pt.size=1
#' @param reduction Which dimensionality reduction to use
#' @param split.by Name of a metadata column to split plot by
#' @param label Whether to label the clusters
#' @param title Title of the plot
#' @param ncol Number of columns for display the plots
#' @param raster Convert points to raster format, default is TRUE
#' @param pt.shape Adjust point shape to plot, default pt.shape=21
#' @param pt.alpha Adjust point alpha to plot, default pt.alpha=0.7
#' @param pt.stroke Adjust point stroke size to plot, default pt.stroke=0.1
#' @seealso \code{\link[Seurat]{DimPlot}}
#' @return A ggplot object
#'
#' @importFrom RANN nn2
#' @importFrom Seurat AddMetaData Embeddings Idents
#' @importFrom cowplot theme_cowplot
#' @importFrom ggplot2 aes element_text geom_point geom_text ggplot 
#'   ggplotGrob ggtitle guide_legend guides labs scale_colour_manual 
#'   scale_x_continuous scale_y_continuous theme
#' @importFrom ggrastr rasterise
#' @importFrom gtable gtable_filter
#' @importFrom patchwork plot_layout wrap_plots
#' @importFrom scales hue_pal
#' @importFrom stats median
#'
#' @references Inspired by Stuart and Butler et al, Cell (2019)
#'
#' @export
#'
#' @examples
#' data("H3N2_small")
#' scDimPlot(object = H3N2_small,
#'   reduction = "umap",
#'   cols = NULL,
#'   split.by = "sample",
#'   ncol = 2,
#'   pt.size = 1
#' )
#'
scDimPlot <- function(object = NULL, cols = NULL, pt.size = 1, 
                      reduction = NULL, split.by = NULL, label = TRUE, 
                      title = NULL, ncol = NULL, raster = TRUE,
                      pt.shape = 21, pt.alpha = 0.7, pt.stroke = 0.1) {
  if (is.null(cols)) {
    cols <- (scales::hue_pal())(length(levels(Seurat::Idents(object))))
  } else if (cols == "sc") {
    if (length(levels(Seurat::Idents(object))) <= 36) {
      cols <- c("#1660A7", "#FF6A00", "#219418", "#CD0C18", "#814BB2", 
        "#794339", "#DC59B6", "#CC79A7", "#FF0000", "#11B3C6", "#AFB400", 
        "#00FFFF", "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", 
        "#0072B2", "#D55E00", "#D55E00", "#CC79A7", "#00AFBB", "#E69F00", 
        "#009E73", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", 
        "#4477AA", "#EE6677", "#228833", "#CCBB44", "#66CCEE", "#AA3377", 
        "#BBBBBB")
    } else {
      warning("Not enough colours provided for ", 
        length(levels(Seurat::Idents(object))),
        " clusters! Use ggplot2's default colors instead\n")
      cols <- (scales::hue_pal())(length(levels(Seurat::Idents(object))))
    }
  } else if (length(levels(Seurat::Idents(object))) > length(cols)) {
    stop("Not enough colours provided for ", 
      length(levels(Seurat::Idents(object))), " clusters!")
  }
  if (is.null(reduction)) {
    if ("umap" %in% names(object)) {
      reduction <- "umap"
    } else if ("tsne" %in% names(object)) {
      reduction <- "tsne"
    } else if ("pca" %in% names(object)) {
      reduction <- "pca"
    } else {
      stop("The reduction parameter does not support! 
        Please use 'umap', 'tsne', or 'pca' instead.\n")
    }
  }
  xmin <- xmax <- ymin <- ymax <- NULL
  ps <- function(data, title = NULL, subtitle = NULL) {
    p <- ggplot2::ggplot(data = data, ggplot2::aes(x = data[, 1], 
      y = data[, 2], z = data[, 3]))
    if (isTRUE(raster)) {
      p <- p + ggrastr::rasterise(dpi = 300, 
        ggplot2::geom_point(ggplot2::aes(fill = data[, 3]), 
          size = pt.size, shape = pt.shape, alpha = pt.alpha,
          stroke = pt.stroke, colour = "black"))
    } else {
      p <- p + ggplot2::geom_point(ggplot2::aes(fill = data[, 3]), 
        size = pt.size, shape = pt.shape, alpha = pt.alpha,
        stroke = pt.stroke, colour = "black")
    }
    if (isTRUE(label)) {
      p <- p + ggplot2::geom_text(ggplot2::aes(label = label), 
        na.rm = TRUE)
    }
    p <- p + ggplot2::scale_fill_manual(values = cols, na.value = "white")
    p <- p + ggplot2::guides(color = ggplot2::guide_legend(title = NULL, 
               override.aes = list(size = 3))) + 
             ggplot2::labs(x = colnames(data)[1], 
               y = colnames(data)[2], title = subtitle)
    p <- p + cowplot::theme_cowplot() + ggplot2::ggtitle(title) + 
      ggplot2::guides(fill = guide_legend(title = "cluster"))
    p <- p + ggplot2::scale_x_continuous(limits = c(xmin, xmax), 
      breaks = seq(floor(xmin/5) * 5, ceiling(xmax/5) * 5, by = 5)) + 
      ggplot2::scale_y_continuous(limits = c(ymin, ymax), 
        breaks = seq(floor(ymin/5) * 5, ceiling(ymax/5) * 5, by = 5))
    return(p)
  }
  pm <- function(data, title = NULL, subtitle = NULL) {
    p <- ggplot2::ggplot(data = data, ggplot2::aes(x = data[, 1], 
      y = data[, 2], z = data[, "color"]))
    if (isTRUE(raster)) {
      p <- p + ggrastr::rasterise(dpi = 300, ggplot2::geom_point(fill = 
        data[, "color"], size = pt.size, shape = pt.shape, alpha = pt.alpha,
        stroke = pt.stroke, colour = "black"))
    } else {
      p <- p + ggplot2::geom_point(fill = data[, "color"], size = pt.size,
        shape = pt.shape, alpha = pt.alpha,stroke = pt.stroke, 
        colour = "black")
    }
    if (isTRUE(label)) {
      p <- p + ggplot2::geom_text(ggplot2::aes(label = label), na.rm = TRUE)
    }
    p <- p + ggplot2::guides(color = ggplot2::guide_legend(title = NULL, 
      override.aes = list(size = 3))) +
      ggplot2::labs(x = colnames(data)[1], 
                    y = colnames(data)[2], 
                    title = title)
    p <- p + cowplot::theme_cowplot() + ggplot2::theme(plot.title = 
      ggplot2::element_text(size = 12))
    p <- p + ggplot2::scale_x_continuous(limits = c(xmin, xmax), 
      breaks = seq(floor(xmin/5) * 5, ceiling(xmax/5) * 5, by = 5)) + 
      ggplot2::scale_y_continuous(limits = c(ymin, ymax), 
        breaks = seq(floor(ymin/5) * 5, ceiling(ymax/5) * 5, by = 5))
    return(p)
  }
  GetXYCenter <- function(data) {
    cluster <- NULL
    result <- data.frame(x = NULL, y = NULL, ident = NULL)
    for (c in as.vector(unique(data[, "cluster"]))) {
      min_x <- max_x <- min_y <- max_y <- center_x <- center_y <- NULL
      sub_data <- subset(data, cluster == c)
      center_x <- median(sub_data[, 1])
      center_y <- median(sub_data[, 2])
      label <- c
      result <- rbind(result, data.frame(x = center_x, y = center_y, 
        label = label))
    }
    return(result)
  }
  reduction_ids <- gsub("coord", toupper(reduction), c("coord_1", "coord_2"))
  coord <- Seurat::Embeddings(object = object, 
                              reduction = reduction)[, c(1, 2)]
  object <- Seurat::AddMetaData(object = object, metadata = coord, 
    col.name = reduction_ids)
  object <- Seurat::AddMetaData(object = object, 
    metadata = as.vector(Seurat::Idents(object)), col.name = "cluster")
  if (is.null(x = split.by)) {
    Data <- object@meta.data[, c(reduction_ids, "cluster")]
    Data[, "cluster"] <- factor(Data[, "cluster"], 
      levels = levels(Seurat::Idents(object)))
    Data[, "color"] <- Data[, "cluster"]
    levels(Data[, "color"]) <- cols[seq_len(length(levels(
      Seurat::Idents(object))))]
    nearest.point <- RANN::nn2(data = Data[, c(1, 2)], 
      query = GetXYCenter(Data)[, c(1, 2)], k = 1)$nn.idx
    Data[, "label"] <- NA
    Data[nearest.point, "label"] <- as.vector(GetXYCenter(Data)[, "label"])
    xmin <- min(Data[, 1])
    xmax <- max(Data[, 1])
    ymin <- min(Data[, 2])
    ymax <- max(Data[, 2])
    return(ps(data = Data, title = title, subtitle = NULL))
  }
  plots <- list()
  if (!is.null(x = split.by)) {
    if (!split.by %in% colnames(object@meta.data)) {
      stop("The parameter 'split.by' ", split.by, 
        " does not exist in MetaData slot!\n")
    }
    Data <- object@meta.data[, c(reduction_ids, "cluster", split.by)]
    Data[, "cluster"] <- factor(Data[, "cluster"], levels = levels(
      Seurat::Idents(object)))
    Data[, "color"] <- Data[, "cluster"]
    levels(Data[, "color"]) <- cols[seq_len(length(levels(
      Seurat::Idents(object))))]
    xmin <- min(Data[, 1])
    xmax <- max(Data[, 1])
    ymin <- min(Data[, 2])
    ymax <- max(Data[, 2])
    if (is.null(ncol)) ncol = ceiling(sqrt(length(unique(Data[, split.by]))))
    legend <- ps(data = Data, title = NULL, subtitle = NULL)
    legend <- ggplot2::ggplotGrob(legend)
    legend <- gtable::gtable_filter(legend, "box", trim = FALSE)
    for (s in unique(Data[, split.by])) {
      data <- Data[Data[, split.by] == s, ]
      nearest.point <- RANN::nn2(data = data[, c(1, 2)], query = 
        GetXYCenter(data)[, c(1, 2)], k = 1)$nn.idx
      data[, "label"] <- NA
      data[nearest.point, "label"] <- as.vector(GetXYCenter(data)[, "label"])
      plots[[s]] <- pm(data, title = s, subtitle = NULL)
    }
  }
  return((patchwork::wrap_plots(plots, ncol = ncol) | 
          legend) + patchwork::plot_layout(widths = c(3, 1)))
}
