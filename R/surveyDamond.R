
#sc = spatialCoords(lit)
#plot(sc, cex=.2, col=factor(lit$cell_category))

#library(shiny)
#library(sosta)
#library(ggplot2)
#
#library(ExperimentHub)
#library(SpatialExperiment)
#library(SingleCellExperiment)
#library(ggplot2)
#library(sf)








#
#
#
# targs = rownames(spe)
# ui = fluidPage(
#  sidebarLayout(
#   sidebarPanel(
#    helpText("pick gene and display presence of cells with relatively 
#         high expression (normalized value exceeds selected quantile)"),
#    selectInput("targ", "target", sort(targs)),
#    numericInput("qthresh", "qthr", value = .75, min=.2, max=.9, step=.1),
#    width=2
#   ),
#   mainPanel(
#    tabsetPanel(
#     tabPanel("view", plotOutput("basic")),
#     tabPanel("by call", plotOutput("bycall")),
#     tabPanel("about", helpText("simple approach"))
#    )
#   )
#  )
# )
# server = function(input, output) {
#  output$basic = renderPlot({
#   })
#  output$bycall = renderPlot({
#   sc = spatialCoords(spe)
#   ndf = data.frame(x=sc[,1], y=sc[,2], type=colData(spe)[[catname]])
#   ggplot(ndf, aes(x=x, y=y, colour=type)) + geom_point()
#  })
#  }
# runApp(list(ui=ui, server=server))
#}



#' shiny app for running sosta on any image in the Damond SPE from sosta package
#' @import shiny
#' @import sosta
#' @import ggplot2
#' @import SpatialExperiment
#' @import SingleCellExperiment
#' @import ggplot2
#' @import sf
#' @param spe SpatialExperiment instance
#' @param qnasssayName character(1) defaults to "quant_norm"
#' @param catname character(1) used to label cells
#' @param markSel character(1) a colData component identifying structure
#' @param imageColName character(1) name of colData component used for filtering
#' @param dim.recon numeric(1) passed to sosta::reconstructShapeDensityImage
#' @param lwd used for plotting boundary of structure
#' @examples
#' if (interactive()) {
#'  requireNamespace("AnnotationHub")
#'  requireNamespace("ExperimentHub")
#'  eh <- try(ExperimentHub::ExperimentHub())
#'  if (inherits(eh, "try-error")) 
#'      eh = ExperimentHub::ExperimentHub(localHub=TRUE)
#'  qout = AnnotationHub::query(eh, "Damond_2019_Pancreas - sce - v1 - full")
#'  stopifnot(length(qout)==1)
#'  oid = names(qout)
#'  # Load single cell experiment object
#'  spe <- eh[[oid]]
#'  # Convert to spatial experiment object
#'  spe <- toSpatialExperiment(spe,
#'      sample_id = "image_name",
#'      spatialCoordsNames = c("cell_x", "cell_y"))
#'  surveyDamond(spe)
#' }
#' @export
surveyDamond = function(spe, qnassayName="quant_norm", catname = "cell_category", markSel="islet", 
   imageColName = "image_name", dim.recon=500, lwd=2) {
 stopifnot(qnassayName %in% assayNames(spe))
 stopifnot(catname %in% names(colData(spe)))

# need to propagate more arg vals up

#spe = spe[, which(spe$image_name == imageId)]
#stopifnot(ncol(spe)>0)
#print(ncol(spe))

if (FALSE) {
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
} # END FALSE BLOCK

 targs = rownames(spe)
 ui = fluidPage(
  sidebarLayout(
   sidebarPanel(
    uiOutput("getid"),
    uiOutput("getimgid"),
    helpText("pick gene and display presence of cells with relatively 
         high expression (normalized value exceeds selected quantile)"),
    selectInput("targ", "target", sort(targs)),
    numericInput("qthresh", "qthr", value = .75, min=.2, max=.9, step=.1),
    width=2
   ),
   mainPanel(
    tabsetPanel(
     tabPanel("view", plotOutput("basic")),
     tabPanel("by call", plotOutput("bycall")),
     tabPanel("about", helpText("simple approach"))
    )
   )
  )
 )
 server = function(input, output) {
  output$getid = renderUI({
   helpText("sample explorer")
   selectInput("id", "id", choices=unique(spe$patient_id))
  })
  output$getimgid = renderUI({
   spe = getdat()
   selectInput("imgid", "image id", choices=unique(colData(spe)[[imageColName]]))
  })
 getdat = reactive({
  validate(need(length(input$id)>0, "waiting for patient id"))
  spe[, which(spe$patient_id == input$id)]
  })


  output$basic = renderPlot({
   spe = getdat()
   validate(need(length(input$imgid)>0, "waiting for image id"))
   spe = spe[, which(colData(spe)[[imageColName]] == input$imgid)]
   xy = spatialCoords(spe)
   assayv = as.numeric(assay(spe[input$targ,], "quant_norm"))
   thr = quantile(assayv, input$qthresh)
   table(assayv>thr)
   validate(need(length(input$id)>0, "waiting for patient id"))
   validate(need(length(input$imgid)>0, "waiting for imgid"))
   sostOut = runsost(spe, imageId = input$imgid, catname = catname, markSel = markSel,
      imageColName = imageColName, dim.recon = dim.recon, lwd=lwd)
#runsost = function(spe, imageId = "E03", catname = "cell_category", markSel="islet",
#  imageColName = "image_name", dim.recon=500, lwd=2) {

   validate(need(!is.null(sostOut), "waiting for sosta"))
   plot(sostOut$struct, reset=FALSE, lwd=lwd)
   points(sostOut$xy[,1], sostOut$xy[,2], pch=19, col=ifelse(assayv>thr, "orange", "lightblue"))
#   lines(sostOut$struct, lwd=lwd)
   })
  output$bycall = renderPlot({
   spe = getdat()
   validate(need(length(input$imgid)>0, "waiting for image id"))
   spe = spe[, which(colData(spe)[[imageColName]] == input$imgid)]
   sc = spatialCoords(spe)
   ndf = data.frame(x=sc[,1], y=sc[,2], type=colData(spe)[[catname]])
   ggplot(ndf, aes(x=x, y=y, colour=type)) + geom_point()
  })
  }
 runApp(list(ui=ui, server=server))
}
