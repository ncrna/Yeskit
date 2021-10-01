#' GO BarPlot visualization for each cluster
#' @title scGOBarPlot visualization
#' @param object Seurat object
#' @param key The 'scGO' results that stored in object@misc$key
#' @param cluster Cluster to presentation
#' @param direction Regulation direction to use, only 'up' or 'down'
#' @param ont Ontology to use, only 'BP', 'CC', or 'MF'
#' @param top_n Only top_n entries to plot
#' @param extra Extra slot used to fetch GOs
#' @return A ggplot object
#'
#' @importFrom ggplot2 aes coord_flip element_blank element_line 
#'   element_text geom_bar geom_hline ggplot ggtitle labs 
#'   scale_y_continuous theme
#' @importFrom stats reorder
#'
#' @export
#'
#' @examples
#' scGOBarPlot(object = H3N2_small,
#'   key = "Infected_vs_Bystander.GO",
#'   ont = "BP",
#'   top_n = 6,
#'   direction = "up",
#'   cluster = "0"
#' )
#'

scGOBarPlot <- function(object = NULL, key = NULL, cluster = NULL, 
                        direction = "up", ont = "BP", top_n = 20, 
                        extra = NULL) {
  Term <- FDR <- NULL
  if (is.null(object)) stop("Parameter 'object' must be specified!\n")
  if (is.null(key)) stop("Parameter 'key' must be specified!\n")
  if (is.null(cluster)) stop("Parameter 'cluster' must be specified!\n")
  cluster = as.character(cluster)
  if (!key %in% names(object@misc)) {
    stop("The 'key' ", key, " does not exist!\n")
  }
  if (is.null(extra)) {
    GOdata = object@misc[[key]]
  } else {
    GOdata = object@misc[[key]][[extra]]
  }
  if (!cluster %in% names(GOdata)) {
    stop("The 'cluster' ", cluster, " does not exist!")
  }
  if (is.null(ont)) ont = "BP"
  if (!ont %in% c("BP", "CC", "MF")) {
    stop("ont must be 'BP', 'CC' or 'MF'.")
  }
  direction <- intersect(direction, c("up", "down"))
  if (length(direction) == 0) {
    warning("Only up-regulated GO terms are used for analysis!")
    direction <- "up"
  }
  colour <- "red"
  if (direction == "down") colour <- "blue"
  z <- GOdata[[cluster]][[direction]][[ont]]
  z <- z[!duplicated(z$GO.ID), ]
  z$Term <- paste(z$GO.ID, z$Term, sep = ": ")
  if (is.null(z)) stop("No data available!\n")
  z <- z[seq_len(top_n), ]
  z <- data.frame(Term = z$Term, 
                  FDR = as.numeric(gsub("< ", "", z$Fisher.elim)))
  p <- ggplot2::ggplot(z, 
    ggplot2::aes(x = stats::reorder(Term, -log10(FDR)), y = -log10(FDR)))
  p <- p + ggplot2::geom_bar(stat = "identity", fill = colour, 
                             color = "NA", width = 0.8)
  p <- p + ggplot2::geom_hline(yintercept = 2, linetype = "dashed", 
                               size = 0.5, color = "grey50")
  p <- p + ggplot2::coord_flip()
  p <- p + ggplot2::theme(panel.grid.major = ggplot2::element_blank(), 
    panel.grid.minor = ggplot2::element_blank(), 
    panel.background = ggplot2::element_blank(),
    axis.title = ggplot2::element_text(size = 12, colour = "black"),
    axis.text = ggplot2::element_text(colour = "black"),
    axis.line.x = ggplot2::element_line(colour = "black"), 
    axis.line.y = ggplot2::element_line(colour = "black"),
    axis.ticks.y = ggplot2::element_blank())
  p <- p + ggplot2::scale_y_continuous(expand = c(0, 0.015))
  p <- p + ggplot2::labs(x = "GO Terms", y = "-Log10(P.adjust)")
  if (is.null(extra)) {
    p <- p + ggplot2::ggtitle(paste0(ifelse(direction == "up", "Up", "Down"),
      "-regulated GO ", ont, " in cluster-", cluster))
  } else {
    p <- p + ggplot2::ggtitle(paste0(ifelse(direction == "up", "Up", "Down"),
      "-regulated GO ", ont, " in cluster-", cluster, " (", extra, ")"))
  }
  return(p)
}
