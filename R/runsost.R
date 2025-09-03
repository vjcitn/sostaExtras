#' perform structure identification on a selected image
#' @param spe SpatialExperiment instance
#' @param imageId character(1) used for filtering
#' @param catname character(1) used to label cells
#' @param markSel character(1) a colData component identifying structure
#' @param imageColName character(1) name of colData component used for filtering
#' @param dim.recon numeric(1) passed to sosta::reconstructShapeDensityImage
#' @param lwd used for plotting boundary of structure
#' @export
runsost = function(spe, imageId = "E03", catname = "cell_category", markSel="islet", 
  imageColName = "image_name", dim.recon=500, lwd=2) {
 stopifnot("quant_norm" %in% assayNames(spe))
 stopifnot(catname %in% names(colData(spe)))

# need to propagate more arg vals up

spe = spe[, which(colData(spe)[[imageColName]] == imageId)]
if (ncol(spe)==0) return(NULL)
print(ncol(spe))

n <- estimateReconstructionParametersSPE(
    spe,
    marks = catname,
    imageCol = imageColName,
    markSelect = markSel,
    plotHist = FALSE
)

struct <- reconstructShapeDensityImage(
    spe,
    marks = catname,
    imageCol = imageColName,
    imageId = imageId,
    markSelect = markSel,
    bndw = n$bndwSPE,
    dim = dim.recon,
    thres = n$thresSPE
)

   xy = spatialCoords(spe)
#   assayv = as.numeric(assay(spe[input$targ,], "quant_norm"))
#   thr = quantile(assayv, input$qthresh)
#   table(assayv>thr)
#   plot(struct, reset=FALSE, lwd=lwd)
#   points(xy[,1], xy[,2], pch=19, col=ifelse(assayv>thr, "orange", "lightblue"))
#   lines(struct, lwd=lwd)
list(xy=xy, struct=struct)
}

