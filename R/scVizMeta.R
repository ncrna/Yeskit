#' visualize data from MetaData slot
#' @title scVizMeta function
#' @description Visualize MetaData in low dimensional embedding
#' @param object Seurat object
#' @param signature Name of one gene set
#' @param pt.size Adjust point size to plot, default pt.size=0.5
#' @param reduction Which dimensionality reduction to use, default is 'umap'
#' @param title Text for the plot title, default title=NULL
#' @param split.by Name of a metadata column to split plot by
#' @param ncol Number of columns for display the plots
#' @param raster Convert points to raster format, default is TRUE
#' @param palette A palette name from RColorBrewer package
#' @param interval set MetaData intervals to plot
#' @return A ggplot object
#'
#' @importFrom RColorBrewer brewer.pal
#' @importFrom cowplot theme_cowplot
#' @importFrom ggplot2 aes geom_point ggplot ggplotGrob guide_legend guides 
#'   labs scale_color_manual scale_x_continuous scale_y_continuous theme
#' @importFrom ggrastr rasterise
#' @importFrom gtable gtable_filter
#' @importFrom patchwork plot_layout wrap_plots
#'
#' @references Inspired by Borcherding et al, F1000Research (2020)
#' @export
#'
#' @examples
#' data("H3N2_small")
#' scVizMeta(object = H3N2_small,
#'   reduction = "umap",
#'   signature = "H3N2",
#'   split.by = "sample",
#'   pt.size = 1,
#'   interval = c(
#'     Abundant = 1000, 
#'     Large = 500, 
#'     Medium = 100, 
#'     Small = 10, 
#'     Single = 1, 
#'     None = 0
#'   )
#' )
#'
scVizMeta <- function(object = NULL, signature = NULL, pt.size = NULL, 
                      reduction = NULL, title = NULL, split.by = NULL, 
                      ncol = NULL, raster = TRUE, palette = "YlOrRd",
                      interval = c(Max = 1e+05, Hyper = 10000, 
                                   Abundant = 1000, Large = 500, 
                                   Medium = 100, Small = 10, 
                                   Single = 1, None = 0)) {
  if (is.null(signature)) {
    stop("Parameter 'signature' must be specified!\n")
  } else if (length(signature) != 1) {
    stop("Parameter 'signature' must be one pathway!\n")
  }
  if (!signature %in% names(object@meta.data)) {
    stop("The signature parameter does not exist in MetaData slot!\n")
  }
  if (is.null(pt.size)) {
    pt.size <- 0.5
  }
  if (is.null(reduction)) {
    if ("umap" %in% names(object)) {
      reduction <- "umap"
    } else if ("tsne" %in% names(object)) {
      reduction <- "tsne"
    } else if ("pca" %in% names(object)) {
      reduction <- "pca"
    } else {
      stop("The reduction parameter does not support!", 
           "Please use 'umap', 'tsne', or 'pca' instead.\n")
    }
  }
  for (x in seq_along(interval)) {
    names(interval)[x] <- paste0(names(interval[x]), " (", interval[x + 1], 
      " < X <= ", interval[x], ")")
  }
  names(interval)[length(interval)] <- 
    gsub("\\(NA < X <", "\\(X", names(interval)[length(interval)])
  xmin <- xmax <- ymin <- ymax <- NULL
  ps <- function(data, cols, title = NULL, legend_title = NULL) {
    p <- ggplot2::ggplot(data = data, 
      ggplot2::aes(x = data[, 1], y = data[, 2]))
    if (isTRUE(raster)) {
      p <- p + ggrastr::rasterise(dpi = 300, 
        ggplot2::geom_point(ggplot2::aes(colour = class), 
          size = pt.size, na.rm = TRUE))
    } else {
    p <- p + ggplot2::geom_point(ggplot2::aes(colour = class), 
      size = pt.size, na.rm = TRUE)
    }
    p <- p + cowplot::theme_cowplot()
    p <- p + ggplot2::guides(colour = 
      ggplot2::guide_legend(override.aes = list(size = 5)))
    p <- p + ggplot2::scale_color_manual(values = cols)
    p <- p + ggplot2::labs(x = colnames(data)[1], 
      y = colnames(data)[2], title = title)
    p <- p + ggplot2::scale_x_continuous(limits = c(xmin, xmax), 
      breaks = seq(floor(xmin/5) * 5, ceiling(xmax/5) * 5, by = 5))
    p <- p + ggplot2::scale_y_continuous(limits = c(ymin, ymax), 
      breaks = seq(floor(ymin/5) * 5, ceiling(ymax/5) * 5, by = 5))
    return(p)
  }
  pm <- function(data, cols = cols, min.value = min.value, 
    max.value = max.value, title = NULL, legend_title = NULL) {
    p <- ggplot2::ggplot(data = data, 
      ggplot2::aes(x = data[, 1], y = data[, 2]))
    if (isTRUE(raster)) {
      p <- p + ggrastr::rasterise(dpi = 300, 
        ggplot2::geom_point(ggplot2::aes(colour = class), 
          size = pt.size, na.rm = TRUE))
    } else {
      p <- p + ggplot2::geom_point(ggplot2::aes(colour = class), 
        size = pt.size, na.rm = TRUE)
    }
    p <- p + cowplot::theme_cowplot() + 
      ggplot2::guides(colour = 
        ggplot2::guide_legend(override.aes = list(size = 5))
      ) + ggplot2::scale_color_manual(values = cols)
    p <- p + ggplot2::labs(x = colnames(data)[1], 
                           y = colnames(data)[2], title = title)
    p <- p + ggplot2::scale_x_continuous(limits = c(xmin, xmax), 
      breaks = seq(floor(xmin/5) * 5, ceiling(xmax/5) * 5, by = 5))
    p <- p + ggplot2::scale_y_continuous(limits = c(ymin, ymax), 
      breaks = seq(floor(ymin/5) * 5, ceiling(ymax/5) * 5, by = 5))
    p <- p + ggplot2::theme(legend.position = "none")
    return(p)
  }
  reduction_ids <- gsub("coord", toupper(reduction), c("coord_1", "coord_2"))
  if (is.null(x = split.by)) {
    Data <- object@meta.data[, c(reduction_ids, signature)]
    Data[, signature][is.na(Data[, signature])] <- 0
    Data[, "class"] <- NA
    for (i in seq_len(length(interval) - 1)) {
      Data$class <- ifelse(Data[, signature] > interval[i + 1] & 
        Data[, signature] <= interval[i], names(interval[i]), Data$class)
    }
    Data$class[is.na(Data$class)] <- names(interval)[length(interval)]
    Data$class <- factor(Data$class, 
      levels = intersect(names(interval), unique(Data$class)))
    Data <- Data[order(Data[, signature]), ]
    xmin <- min(Data[, 1])
    xmax <- max(Data[, 1])
    ymin <- min(Data[, 2])
    ymax <- max(Data[, 2])
    #cols <- c(rev(RColorBrewer::brewer.pal(length(levels(Data$class)) - 1, 
    #  palette))[1:length(levels(Data$class)) -1], "gray")
    cols <- c(rev(RColorBrewer::brewer.pal(length(levels(Data$class)) - 1, 
      palette))[seq_len(length(levels(Data$class)) -1)], "gray")
    # names(cols) <- names(interval)
    names(cols) <- names(interval)[(length(interval) - 
      length(levels(Data$class)) + 1):length(interval)]
    return(ps(data = Data, cols = cols, 
      title = title, legend_title = signature))
  }
  plots <- list()
  if (!is.null(x = split.by)) {
    if (!split.by %in% colnames(object@meta.data)) {
      stop("The parameter 'split.by' ", split.by, 
           " does not exist in MetaData slot!\n")
    }
    Data <- object@meta.data[, c(reduction_ids, signature, split.by)]
    Data[, signature][is.na(Data[, signature])] <- 0
    Data[, "class"] <- NA
    for (i in seq_len(length(interval) - 1)) {
      Data$class <- ifelse(Data[, signature] > interval[i + 1] & 
        Data[, signature] <= interval[i], names(interval[i]), Data$class)
    }
    Data$class[is.na(Data$class)] <- names(interval)[length(interval)]
    Data$class <- factor(Data$class, levels = intersect(names(interval), 
      unique(Data$class)))
    Data <- Data[order(Data[, signature]), ]
    xmin <- min(Data[, 1])
    xmax <- max(Data[, 1])
    ymin <- min(Data[, 2])
    ymax <- max(Data[, 2])
    #cols <- c(rev(RColorBrewer::brewer.pal(length(levels(Data$class)) - 1, 
    #  palette))[1:length(levels(Data$class)) -1], "gray")
    cols <- c(rev(RColorBrewer::brewer.pal(length(levels(Data$class)) - 1, 
      palette))[seq_len(length(levels(Data$class)) -1)], "gray")
    names(cols) <- names(interval)
    legend <- ps(data = Data, cols = cols, title = NULL, 
                 legend_title = signature)
    legend <- ggplot2::ggplotGrob(legend)
    legend <- gtable::gtable_filter(legend, "box", trim = FALSE)
    if (is.null(ncol)) {
      ncol = ceiling(sqrt(length(unique(Data[, split.by]))))
    }
    for (s in unique(Data[, split.by])) {
      data <- Data[Data[, split.by] == s, ]
      plots[[s]] <- pm(data = data, cols = cols, title = s, 
                       legend_title = signature)
    }
  }
  return((patchwork::wrap_plots(plots, ncol = ncol) | legend) + 
          patchwork::plot_layout(widths = c(3, 1)))
}
