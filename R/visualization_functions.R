#' Plot expression overlay on structure
#'
#' Creates a plot showing the reconstructed structure with cells colored by
#' gene expression levels (above/below a quantile threshold).
#'
#' @param spe SpatialExperiment object containing the spatial data
#' @param target_gene character(1) name of the gene to visualize
#' @param quantile_threshold numeric(1) quantile threshold for expression (default 0.75)
#' @param sost_output list containing xy coordinates and struct from runsost()
#' @param lwd numeric(1) line width for structure boundary (default 2)
#' @param assay_name character(1) name of assay to use (default "quant_norm")
#' @return Plot showing structure with expression overlay (side effect)
#' @examples
#' \dontrun{
#' # Assuming you have spe and sost_output from runsost()
#' plot_expression_overlay(spe, "INS", 0.75, sost_output, lwd = 2)
#' }
#' @export
plot_expression_overlay <- function(spe, target_gene, quantile_threshold = 0.75, 
                                   sost_output, lwd = 2, assay_name = "quant_norm") {
  # Validate inputs
  if (is.null(sost_output)) {
    stop("sost_output cannot be NULL")
  }
  if (!target_gene %in% rownames(spe)) {
    stop("target_gene '", target_gene, "' not found in spe rownames")
  }
  if (!assay_name %in% assayNames(spe)) {
    stop("assay_name '", assay_name, "' not found in spe")
  }
  
  # Extract expression values and calculate threshold
  assayv <- as.numeric(assay(spe[target_gene, ], assay_name))
  thr <- quantile(assayv, quantile_threshold)
  
  # Create the plot
  plot(sost_output$struct, reset = FALSE, lwd = lwd)
  points(sost_output$xy[, 1], sost_output$xy[, 2], 
         pch = 19, 
         col = ifelse(assayv > thr, "orange", "lightblue"))
}


#' Plot cell types by spatial coordinates
#'
#' Creates a ggplot showing the spatial distribution of cells colored by
#' their categorical annotation (e.g., cell type).
#'
#' @param spe SpatialExperiment object containing the spatial data
#' @param catname character(1) name of colData column containing cell categories
#' @return ggplot object showing spatial coordinates colored by cell type
#' @examples
#' \dontrun{
#' # Assuming you have filtered spe for a specific image
#' p <- plot_cell_types(spe, "cell_category")
#' print(p)
#' }
#' @import ggplot2
#' @export
plot_cell_types <- function(spe, catname) {
  # Validate inputs
  if (!catname %in% names(colData(spe))) {
    stop("catname '", catname, "' not found in spe colData")
  }
  
  # Get spatial coordinates
  sc <- spatialCoords(spe)
  
  # Create data frame for plotting
  ndf <- data.frame(
    x = sc[, 1], 
    y = sc[, 2], 
    type = colData(spe)[[catname]]
  )
  
  # Create ggplot
  ggplot(ndf, aes(x = x, y = y, colour = type)) + 
    geom_point()
}