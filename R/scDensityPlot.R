#' @title scDensityPlot function
#' @description Visualize density plot in low dimensional embedding
#' @param object Seurat object
#' @param reduction Dimensional reduction to use, default reduction='umap'
#' @param title Figure title
#' @param split.by Name of a meta.data column to split plot by
#' @param ncol Number of columns for display the plots
#' @param colors Vector of colours to use for n-colour gradient
#' @param digits Integer representing the number of decimal places reserved for
#' the legend, default digits=3
#'
#' @return A ggplot object
#'
#' @importFrom ggplot2 aes ggplot ggplotGrob guide_legend guides labs 
#'   scale_fill_gradientn scale_x_continuous scale_y_continuous 
#'   stat_density2d theme
#' @importFrom patchwork wrap_plots plot_layout
#' @importFrom gtable gtable_filter
#' @importFrom MASS kde2d bandwidth.nrd
#' @importFrom Seurat Embeddings
#' @importFrom cowplot theme_cowplot
#'
#' @export
#'
#' @examples
#' scDensityPlot(object = H3N2_small,
#'   reduction = "umap",
#'   split.by = "sample",
#'   ncol = 2
#' )
#'
scDensityPlot <- function(object = NULL, reduction = NULL, title = NULL, 
                          split.by = NULL, ncol = NULL, colors = NULL, 
                          digits = NULL) {
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
  if (is.null(title)) {
    title <- "Density plot"
  }
  if (is.null(colors)) {
    colors <- c("white", "gray99", "#FFEDA0", "red", "darkred")
  }
  if (is.null(digits)) digits <- 3
  xmin <- xmax <- ymin <- ymax <- z.min <- z.max <- ..density.. <- NULL
  ps <- function(data, z.min = z.min, z.max = z.max, title = NULL) {
    p <- ggplot2::ggplot(data = data, 
      ggplot2::aes(x = data[, 1], y = data[, 2]))
    p <- p + ggplot2::stat_density2d(geom = "raster", interpolate = TRUE, 
      ggplot2::aes(fill = ..density..), contour = FALSE, na.rm = TRUE) +
      ggplot2::scale_fill_gradientn(colors = colors, 
        limits = c(z.min, z.max), 
        breaks = seq(z.min, z.max, length = 5), 
        labels = round(seq(z.min, z.max, length = 5), digits = digits)) + 
      cowplot::theme_cowplot() + ggplot2::guides(colour = 
        ggplot2::guide_legend(override.aes = list(size = 5))) +
      ggplot2::labs(x = colnames(data)[1], 
                    y = colnames(data)[2], 
                    title = title)
    p <- p + ggplot2::scale_x_continuous(limits = c(xmin, xmax), 
      breaks = seq(floor(xmin/5) * 5, ceiling(xmax/5) * 5, by = 5)) + 
      ggplot2::scale_y_continuous(limits = c(ymin, ymax), breaks = 
      seq(floor(ymin/5) * 5, ceiling(ymax/5) * 5, by = 5))
    return(p)
  }
  pm <- function(data, z.min = z.min, z.max = z.max, title = NULL) {
    p <- ggplot2::ggplot(data = data, 
      ggplot2::aes(x = data[, 1], y = data[, 2]))
    p <- p + ggplot2::stat_density2d(geom = "raster", interpolate = TRUE, 
      ggplot2::aes(fill = ..density..), contour = FALSE, na.rm = TRUE) + 
      ggplot2::scale_fill_gradientn(colors = colors, 
        limits = c(z.min, z.max), breaks = seq(z.min, z.max, length = 5), 
        labels = round(seq(z.min, z.max, length = 5), digits = digits)) + 
      cowplot::theme_cowplot() + ggplot2::guides(colour = 
        ggplot2::guide_legend(override.aes = list(size = 5))) +
      ggplot2::labs(x = colnames(data)[1], 
                    y = colnames(data)[2], 
                    title = title)
    p <- p + ggplot2::scale_x_continuous(limits = c(xmin, xmax), 
      breaks = seq(floor(xmin/5) * 5, ceiling(xmax/5) * 5, by = 5)) + 
      ggplot2::scale_y_continuous(limits = c(ymin, ymax), 
        breaks = seq(floor(ymin/5) * 5, ceiling(ymax/5) * 5, by = 5))
    p <- p + ggplot2::theme(legend.position = "none")
    return(p)
  }
  Data <- Seurat::Embeddings(object = object[[reduction]])[, c(1, 2)]
  Data <- as.data.frame(Data)
  m <- MASS::kde2d(x = Data[, 1], y = Data[, 2], 
    h = c(MASS::bandwidth.nrd(Data[, 1]), MASS::bandwidth.nrd(Data[, 2])),
    n = 100, lims = c(range(Data[, 1]), range(Data[, 2])))
    z.min <- min(m$z)
    z.max <- max(m$z)
    xmin <- min(Data[, 1])
    xmax <- max(Data[, 1])
    ymin <- min(Data[, 2])
    ymax <- max(Data[, 2])
  if (is.null(x = split.by)) return(ps(Data, z.min, z.max, title))
  plots <- list()
  z <- vector()
  if (!is.null(x = split.by)) {
    Data[, split.by] <- object[[split.by, drop = FALSE]]
    if (is.null(ncol)) {
      ncol = ceiling(sqrt(length(unique(Data[, split.by]))))
    }
    for (s in unique(Data[, split.by])) {
      data <- Data[Data[, 3] == s, ]
      m <- MASS::kde2d(x = data[, 1], y = data[, 2], 
        h = c(MASS::bandwidth.nrd(data[, 1]), 
        MASS::bandwidth.nrd(data[, 2])), n = 100, 
        lims = c(range(data[, 1]), range(data[, 2])))
      z <- c(z, range(m$z))
    }
    z <- sort(z)
    z.min <- z[1]
    z.max <- z[length(z)]
    xmin <- min(Data[, 1])
    xmax <- max(Data[, 1])
    ymin <- min(Data[, 2])
    ymax <- max(Data[, 2])
    legend <- ps(Data, z.min, z.max, title = NULL)
    legend <- ggplot2::ggplotGrob(legend)
    legend <- gtable::gtable_filter(legend, "box", trim = FALSE)
    for (s in unique(Data[, split.by])) {
      data <- Data[Data[, 3] == s, ]
      plots[[s]] <- pm(data, z.min, z.max, title = s)
    }
  }
  return((patchwork::wrap_plots(plots, ncol = ncol) | legend) + 
          patchwork::plot_layout(widths = c(3, 1)))
}
