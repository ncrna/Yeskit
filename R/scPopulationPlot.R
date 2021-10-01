#' Population plot by clusters or samples
#' @title scPopulationPlot
#' @param object Seurat object
#' @param by String used to separate cluster/sample, 
#'   only 'cluster' or 'sample' is accepted
#' @param order A vector to set the order of samples
#' @param cols Vector of colors, each color corresponds to an identity class. 
#'   Use ggplot2's default colors by default.
#' We include a pallete called 'sc' which consists of 36 colors
#' @return A dataframe that used to draw this plot
#'
#' @importFrom ggplot2 aes element_blank element_text geom_bar ggplot 
#'   guide_legend position_fill scale_fill_manual scale_y_continuous 
#'   theme xlab ylab
#' @importFrom Seurat Idents
#' @importFrom cowplot theme_cowplot
#' @importFrom scales hue_pal
#'
#' @export
#'
#' @examples
#' data("H3N2_small")
#' scPopulationPlot(object = H3N2_small,
#'   by = "sample",
#'   cols = NULL,
#'   order = c("Bystander", "Infected")
#' )
#'
scPopulationPlot <- function(object = NULL, by = c("sample", "cluster"), 
                             order = NULL, cols = NULL) {
  if (is.null(by)) by = "sample"
  if (is.null(order)) order = levels(factor(object$sample))
  if (is.null(cols)) {
    cols <- (scales::hue_pal())(length(levels(Seurat::Idents(object))))
  } else if (cols == "sc") {
    if (length(levels(Seurat::Idents(object))) <= 36) {
      cols <- c("#1660A7", "#FF6A00", "#219418", "#CD0C18", "#814BB2", 
                "#794339", "#DC59B6", "#CC79A7", "#FF0000", "#11B3C6", 
                "#AFB400", "#00FFFF", "#999999", "#E69F00", "#56B4E9", 
                "#009E73", "#F0E442", "#0072B2", "#D55E00", "#D55E00", 
                "#CC79A7", "#00AFBB", "#E69F00", "#009E73", "#56B4E9", 
                "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#4477AA",
                "#EE6677", "#228833", "#CCBB44", "#66CCEE", "#AA3377", 
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
  Sample <- Cells <- Cluster <- NULL
  data <- table(Seurat::Idents(object), object$sample)
  ClusterNames <- rownames(data)
  data <- as.data.frame(data)
  names(data) <- c("Cluster", "Sample", "Cells")
  data$Cluster <- as.character(data$Cluster)
  data$Cluster <- factor(data$Cluster, levels = ClusterNames)
  data$Sample <- factor(data$Sample, levels = order)
  if (by == "sample") {
    p <- ggplot2::ggplot(data = data, 
      ggplot2::aes(x = Sample, y = Cells, fill = Cluster))
    p <- p + ggplot2::geom_bar(stat = "identity", 
      position = ggplot2::position_fill(reverse = FALSE), 
        width = 0.6, size = 0.3, colour = NA)
    p <- p + ggplot2::xlab("Sample") + ggplot2::ylab("Fraction of cells (%)")
    p <- p + cowplot::theme_cowplot()
    p <- p + ggplot2::theme(axis.text.x = ggplot2::element_text(size = 10,
      angle = 45, vjust = 1, hjust = 1), 
      axis.title.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = 10), 
      axis.title.y = ggplot2::element_text(size = 12),
      panel.grid.minor = ggplot2::element_blank())
    p <- p + ggplot2::scale_fill_manual(values = cols,
      guide = ggplot2::guide_legend(keywidth = 0.5, keyheight = 0.5, ncol = 1))
    p <- p + ggplot2::scale_y_continuous(expand = c(0.01, 0.01), 
      labels = c(0, 25, 50, 75, 100))
    } else if (by == "cluster") {
      p <- ggplot2::ggplot(data = data, ggplot2::aes(x = Cluster, y = Cells, 
        fill = Sample))
      p <- p + ggplot2::geom_bar(stat = "identity", 
        position = ggplot2::position_fill(reverse = FALSE),
          width = 0.6, size = 0.3, colour = NA)
      p <- p + ggplot2::xlab("Cluster") + ggplot2::ylab("Fraction of cells (%)")
      p <- p + cowplot::theme_cowplot()
      p <- p + ggplot2::theme(axis.text.x = ggplot2::element_text(size = 10,
        angle = 45, vjust = 1, hjust = 1), 
        axis.title.x = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_text(size = 10), 
        axis.title.y = ggplot2::element_text(size = 12),
        panel.grid.minor = ggplot2::element_blank())
      p <- p + ggplot2::scale_fill_manual(values = cols, 
        guide = ggplot2::guide_legend(keywidth = 0.5, 
                                      keyheight = 0.5, 
                                      ncol = 1))
      p <- p + ggplot2::scale_y_continuous(expand = c(0.01, 0.01), 
        labels = c(0, 25, 50, 75, 100))
  } else {
    stop("Parameter 'by' must be 'cluster' or 'sample'!\n")
  }
  return(p)
}
