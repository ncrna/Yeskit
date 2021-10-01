#' scGODotPlot visualization for all clusters
#' @title scGODotPlot visualization
#' @param object Seurat object
#' @param key The keyword that is stored in object@misc$key
#' @param clusters Vectors of clusters to display, default clusters=NULL 
#'   means all clusters
#' @param direction Regulation direction to use, only 'up' or 'down' is allowed
#' @param ont Ontology to use, only 'BP', 'CC', or 'MF' is allowed
#' @param padj.cutoff The cutoff of adjusted Pvalue, default padj.cutoff=0.01
#' @param order.genes Order GO results by the number of significant genes, 
#'   default is FALSE
#' @param colors Colors to specify, default cols = c('white', 'red')
#' @param font.size Font size of axis
#' @param top_n Only top_n entries to plot
#' @param extra Extra slot used to fetch GOs
#' @return A ggplot object
#'
#' @importFrom ggplot2 aes element_blank element_line element_rect 
#'   element_text geom_point ggplot ggtitle labs scale_colour_gradient 
#'   scale_y_discrete theme
#'
#' @export
#'
#' @examples
#' scGODotPlot(object = H3N2_small,
#'   key = "Infected_vs_Bystander.GO",
#'   ont = "BP",
#'   top_n = 6,
#'   direction = "up",
#'   cluster = "0",
#'   font.size = 8
#' )
#'
scGODotPlot <- function(object = NULL, key = NULL, clusters = NULL, 
                        direction = "up", ont = "BP", padj.cutoff = 0.01, 
                        order.genes = FALSE, colors = NULL, font.size = 8,
                        top_n = 5, extra = NULL) {
  if (is.null(object)) stop("Parameter 'object' must be specified!\n")
  if (is.null(key)) stop("Parameter 'key' must be specified!\n")
  if (!key %in% names(object@misc)) 
    stop("The 'key' ", key, " does not exist!\n")
  direction <- intersect(direction, c("up", "down"))
  if (length(direction) == 0) {
    warning("Only up-regulated GO terms are used for analysis!")
    direction <- "up"
  }
  if (is.null(colors)) {
    if (direction == "up") {
      colors <- c("#FFC8C8", "#FF0000")
    } else {
      colors <- c("#C8C8FF", "#0000FF")
    }
  } else if (length(colors) < 2) {
    stop("Parameter 'cols' must be specified with two colors!\n")
  }
  if (is.null(font.size)) font.size = 8
  Cluster <- Term <- Value <- NULL
  if (is.null(extra)) {
    GOdata = object@misc[[key]]
  } else {
    GOdata = object@misc[[key]][[extra]]
  }
  if (is.null(top_n)) top_n = 5
  if (!is.null(clusters)) {
    clusters <- as.character(clusters)
    clusters <- intersect(clusters, names(GOdata))
  } else {
    clusters <- names(GOdata)
  }
  GO.terms = vector()
  for (i in clusters) {
    GO.tmp <- GOdata[[i]][[direction]][[ont]]
    if (is.null(GO.tmp)) next
    if (isTRUE(order.genes)) {
      GO.tmp <- GO.tmp[order(as.numeric(GO.tmp$Significant), 
        decreasing = TRUE), ]
      GO.tmp <- GO.tmp[as.numeric(gsub("< ", "", 
        GO.tmp$Fisher.elim)) <= padj.cutoff, ]
    }
    term <- GO.tmp$Term[seq_len(top_n)]
    term <- term[!is.na(term)]
    id <- GO.tmp$GO.ID[seq_len(top_n)]
    id <- id[!is.na(id)]
    term <- paste(id, term, sep = ": ")
    GO.terms <- unique(c(GO.terms, term))
  }
  GO.data <- data.frame(Cluster = character(), Term = character(), 
    GeneNum = character(), Value = numeric())
  for (i in clusters) {
    z = GOdata[[i]][[direction]][[ont]]
    if (is.null(z)) next
    z$Term <- paste(z$GO.ID, z$Term, sep = ": ")
    for (GO.term in GO.terms) {
      if (nrow(subset(z, Term == GO.term)) == 0) {
        GO.value = NA
        GeneNum = 0
      } else {
        GO.value = -log10(as.numeric(gsub("< ", "", 
          subset(z, Term == GO.term)[1, ]$Fisher.elim)))
        GeneNum = as.integer(subset(z, Term == GO.term)[1, ]$Significant)
      }
      GO.value <- ifelse(is.na(GO.value), 0, GO.value)
      GO.row <- data.frame(Cluster = i, Term = GO.term, GeneNum = GeneNum,
        Value = GO.value)
      GO.data <- rbind(GO.data, GO.row)
    }
  }
  if (nrow(GO.data) == 0) {
    stop("There is no infomation available for the clusters you specified!")
  }
  GO.data <- GO.data[GO.data$GeneNum > 0, ]
  # fix the order of go terms
  GO.data$Term <- factor(GO.data$Term, levels = GO.terms)
  # fix the order of clusters
  GO.data$Cluster <- factor(GO.data$Cluster, levels = clusters)
  breaks <- c(seq(-log10(padj.cutoff), ceiling(max(GO.data$Value)), 
    length.out = 5))
  p <- ggplot2::ggplot(GO.data, ggplot2::aes(x = Cluster, y = Term)) 
  p <- p + ggplot2::geom_point(ggplot2::aes(colour = Value, size = GeneNum), 
    stat = "identity")
  p <- p + ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", vjust = 1), 
    axis.title.x = ggplot2::element_text(face = "bold", vjust = -0.2),
    axis.title.y = ggplot2::element_text(face = "bold"), 
    axis.text.y = ggplot2::element_text(hjust = 1, colour = "black"), 
    axis.text.x = ggplot2::element_text(angle = -45, vjust = 1, 
      hjust = 0, colour = "black"), 
    panel.background = ggplot2::element_blank(), 
    axis.ticks = ggplot2::element_blank(), 
    axis.text = ggplot2::element_text(size = font.size),
    panel.border = ggplot2::element_rect(colour = "black", 
      fill = NA, size = 1),
    panel.grid.minor = ggplot2::element_line(colour = "grey", 
      linetype = "dotted", size = 0.1), 
    panel.grid.major = ggplot2::element_line(colour = "grey",
      linetype = "dotted", size = 0.1)) +
    ggplot2::scale_colour_gradient(low = colors[1], high = colors[2], 
      na.value = "white", breaks = breaks, labels = as.character(breaks),
      limits = c(-log10(padj.cutoff), max(breaks))) + 
    ggplot2::labs(x = "Cluster", y = "GO Terms", color = "-Log10(P.adjust)") + 
    ggplot2::scale_y_discrete(limits = rev(levels(GO.data$Term)))
    if (is.null(extra)) {
      p <- p + ggplot2::ggtitle(paste0(ifelse(direction == "up", "Up", "Down"),
        "-regulated GO ", ifelse(ont == "BP", "Biological Process", 
        ifelse(ont == "CC", "Cellular Component", "Molecular Function"))))
    } else {
      p <- p + ggplot2::ggtitle(paste0(ifelse(direction == "up", "Up", "Down"),
        "-regulated GO ", ifelse(ont == "BP", "Biological Process", 
        ifelse(ont == "CC", "Cellular Component", "Molecular Function")), 
        " (", extra, ")"))
    }
  return(p)
}
