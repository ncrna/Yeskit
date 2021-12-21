#' pathogen population plot by clusters
#' @title scPathogenRatioPlot function
#' @param object Seurat object
#' @param species Name of pathogen species
#' @param cols Colors to plot
#' @param split.by Name of a metadata column to split plot by
#' @param ncol Number of columns for display the plots
#' @return A ggplot object
#'
#' @importFrom cowplot theme_cowplot
#' @importFrom ggplot2 aes geom_bar ggplot ggplotGrob labs 
#'   scale_fill_manual scale_x_discrete scale_y_continuous theme
#' @importFrom gtable gtable_filter
#' @importFrom patchwork plot_layout wrap_plots
#'
#' @export
#'
#' @examples
#' data("H3N2_small")
#' scPathogenRatioPlot(object = H3N2_small,
#'   species = "H3N2",
#'   split.by = "sample",
#'   ncol = 2
#' )
#'
scPathogenRatioPlot <- function(object = NULL, species = NULL, cols = NULL, 
                                split.by = NULL, ncol = NULL) {
  options(scipen = 5)
  if (is.null(cols)) cols <- c("gray80", "red")
  Cluster <- Frequency <- Status <- NULL
  pp <- function(data = NULL, title = NULL, Legend = FALSE) {
    data.Positive <- data.frame(colSums(table(data[data[, species] > 0, ])),
      Status = "Positive")
    if (nrow(data[data[, species] > 0, ]) == 0) {
      data.Positive <- data.frame(colSums(table(data[data[, species] > 0, ])),
        Counts = 0, Status = "Positive")
    }
    data.Positive[, "Cluster"] <- factor(rownames(data.Positive), 
      levels = as.character(rownames(data.Positive)))
    colnames(data.Positive) <- c("Frequency", "Status", "Cluster")
    data.Negative <- data.frame(colSums(table(data[data[, species] == 0, ])),
      Status = "Negative")
    if (nrow(data[data[, species] == 0, ]) == 0) {
      data.Negative <- data.frame(colSums(table(data[data[, species] == 0, ])),
        Counts = 0, Status = "Negative")
    }
    data.Negative[, "Cluster"] <- factor(rownames(data.Negative), 
      levels = as.character(rownames(data.Negative)))
    colnames(data.Negative) <- c("Frequency", "Status", "Cluster")
    data <- rbind(data.Positive, data.Negative)
    data[, "Status"] <- factor(data[, "Status"], 
      levels = c("Negative", "Positive"))
    p <- ggplot2::ggplot(data, 
      ggplot2::aes(x = Cluster, y = Frequency, fill = Status)) +
      ggplot2::geom_bar(position = "fill", stat = "identity") + 
      cowplot::theme_cowplot() + ggplot2::scale_fill_manual(values = cols) + 
      ggplot2::scale_x_discrete(expand = c(0.05, 0.05)) + 
      ggplot2::scale_y_continuous(expand = c(0.01, 0.01), 
        labels = c(0, 25, 50, 75, 100)) + 
      ggplot2::labs(x = "Clusters", y = "Fraction of cells (%)", title = title)
    if (isFALSE(Legend)) {
      p <- p + ggplot2::theme(legend.position = "none")
    }
    return(p)
  }
  if (is.null(x = split.by)) {
    Data <- object@meta.data[, c(species, "clusters")]
    return(pp(Data, title = NULL, Legend = TRUE))
  }
  plots <- list()
  if (!is.null(x = split.by)) {
    if (!split.by %in% colnames(object@meta.data)) {
      stop("The parameter 'split.by' ", split.by, 
           " does not exist in MetaData slot!\n")
    }
    Data <- object@meta.data[, c(species, "clusters", split.by)]
    if (is.null(ncol)) {
      ncol = ceiling(sqrt(length(unique(Data[, split.by]))))
    }
    legend <- pp(object@meta.data[, c(species, "clusters")], 
      title = NULL, Legend = TRUE)
    legend <- ggplot2::ggplotGrob(legend)
    legend <- gtable::gtable_filter(legend, "box", trim = FALSE)
    for (s in unique(Data[, split.by])) {
      data <- Data[Data[, split.by] == s, ]
      plots[[s]] <- pp(data, title = s, Legend = FALSE)
    }
  }
  return((patchwork::wrap_plots(plots, ncol = ncol) | legend) + 
          patchwork::plot_layout(widths = c(3, 1)))
}
