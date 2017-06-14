# 01-kmeans-app

palette(c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
  "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999"))

library(shiny)
library(geomorph)
library(ggplot2)
# rsconnect::deployApp('path/to/your/app')

# setwd("~/Documents/Website/")
source("readShinyData.R")

ui <- fluidPage(
  titlePanel("Choose your own analysis: Developmental System Drift edition"),
  column(3, wellPanel(
                radioButtons(
                  "cells", "Embryo Stage",
                  choices = c("2-cell" = "embryo2.3D",
                              "4-cell" = "embryo4.3D"),
                  selected = "embryo4.3D"),
                
                selectInput("analysis", "Analysis", 
                  c("landmarks","procrustes shape","PCA","cell sizes","trajectory analysis")), 
                
                uiOutput("ui")
                )),
           
  mainPanel(
    plotOutput("plot1")
  )
)

server <- function(input, output) {
  
  output$ui <- renderUI({
    if (is.null(input$analysis))
      return()
    
    # Depending on input$input_type, we'll generate a different
    # UI component and send it to the client.
    switch(input$analysis,
           "landmarks" = radioButtons("dynamic", "Options",
                                      choices = c("Specimens" = "specimens",
                                                  "Labels" = "labels"),
                                      selected = "labels"),
           
           "procrustes shape" = fluidRow(
                                  radioButtons("reference", "Reference Shape",
                                                     choices = c("elegans" = "elegans1-20",
                                                                 "briggsae" = "briggsae1-20",
                                                                 "elegans (stressed)" = "elegans1-30",
                                                                 "briggsae (stressed)" = "briggsae1-30"),
                                                      selected = "elegans1-20"),

                                  radioButtons("target", "Target Shape",
                                               choices = c("elegans" = "elegans2-20",
                                                           "briggsae" = "briggsae2-20",
                                                           "elegans (stressed)" = "elegans2-30",
                                                           "briggsae (stressed)" = "briggsae2-30"),
                                               selected = "briggsae2-20"),
                                  
                                  sliderInput("magnitude", "Exaggeration",
                                              min = 1, max = 5, value = 2, step = 1),
                                  
                                  selectInput("method", "Method",
                                              list("Points" = "points", 
                                                   "Vectors" = "vector", 
                                                   "Thin Plate Spline" = "TPS"),
                                              selected = "points")
                                  
                                  ),
           
           "PCA" =  fluidRow(
                      sliderInput("pcX", "X-Axis",
                                  min = 1, max = 10, value = 1, step = 1),

                      sliderInput("pcY", "Y-Axis",
                                  min = 1, max = 10, value = 2, step = 1)
                            ),
           
           "cell sizes" = fluidRow(
                      checkboxGroupInput("cellTypes", "Cells",
                                  choices = unique(cellNum()$type), 
                                  selected = unique(cellNum()$type)),
                                  
                      checkboxGroupInput("cellSpecies", "Comparisons",
                                  choices = unique(cellNum()$species),
                                  selected = unique(cellNum()$species))
                          ),
           
           "trajectory analysis" = checkboxInput("bp", "Add Color",
                                                  value = TRUE)
           )

  })
  
  
  embryo <- reactive({
    x <- get(input$cells)
  })
  
  outline <- reactive({
    switch(input$cells,
           embryo2.3D = cell.link.2,
           embryo4.3D = cell.link.4)
  })
  
  landmarks <- reactive({
    switch(input$cells,
           embryo2.3D = 12,
           embryo4.3D = 18)
  })
  
  dynamic <- reactive({
    switch(input$dynamic,
          "specimens" = c(FALSE,0.8),
          "labels" = c(TRUE,0) )
  })
  

  ref <- reactive({
    switch(input$reference,
           "elegans1-20" =  arrayspecs((rowsum(subset(two.d.array(embryo()), devTemp == "20C"), subset(ID, devTemp=="20C"))/as.vector(table(subset(ID, devTemp=="20C")))),landmarks(),2)[,,"N2"],
           "briggsae1-20" = arrayspecs((rowsum(subset(two.d.array(embryo()), devTemp == "20C"), subset(ID, devTemp=="20C"))/as.vector(table(subset(ID, devTemp=="20C")))),landmarks(),2)[,,"AF16"],
           "elegans1-30" =  arrayspecs((rowsum(subset(two.d.array(embryo()), devTemp == "30C"), subset(ID, devTemp=="30C"))/as.vector(table(subset(ID, devTemp=="30C")))),landmarks(),2)[,,"N2"],
           "briggsae1-30" = arrayspecs((rowsum(subset(two.d.array(embryo()), devTemp == "30C"), subset(ID, devTemp=="30C"))/as.vector(table(subset(ID, devTemp=="30C")))),landmarks(),2)[,,"AF16"])
  })
  
  tar <- reactive({
    switch(input$target,
           "elegans2-20" = arrayspecs((rowsum(subset(two.d.array(embryo()), devTemp == "20C"), subset(ID, devTemp=="20C"))/as.vector(table(subset(ID, devTemp=="20C")))),landmarks(),2)[,,"N2"],
           "briggsae2-20" = arrayspecs((rowsum(subset(two.d.array(embryo()), devTemp == "20C"), subset(ID, devTemp=="20C"))/as.vector(table(subset(ID, devTemp=="20C")))),landmarks(),2)[,,"AF16"],
           "elegans2-30" = arrayspecs((rowsum(subset(two.d.array(embryo()), devTemp == "30C"), subset(ID, devTemp=="30C"))/as.vector(table(subset(ID, devTemp=="30C")))),landmarks(),2)[,,"N2"],
           "briggsae2-30" = arrayspecs((rowsum(subset(two.d.array(embryo()), devTemp == "30C"), subset(ID, devTemp=="30C"))/as.vector(table(subset(ID, devTemp=="30C")))),landmarks(),2)[,,"AF16"])
  })
  

  
  TA <- reactive({
    switch(input$cells,
           embryo2.3D = TA2,
           embryo4.3D = TA4)
  })
  
  pts <- reactive({
    switch(input$bp,
           "FALSE" = c("darkgreen", "lightgray", "red"),
           "TRUE" = c("black", "gray","white") )
  })
  
  cellNum <- reactive({
    switch(input$cells,
           embryo2.3D = cells2,
           embryo4.3D = cells4)
  })


  C2 <- reactive(subset(cellNum(), cellNum()$type %in% input$cellTypes & cellNum()$species %in% input$cellSpecies))
  
  
  output$plot1 <- renderPlot({
    embryo <- embryo()
    outline <- outline()
    dynamic <- dynamic()
    plot1 <- switch(input$analysis,
                    "landmarks" = plotAllSpecimens(embryo, links = outline, label = dynamic[1], mean = TRUE, plot.param = list(pt.bg = as.factor(paste(ID, devTemp, sep="-")), pt.cex = dynamic[2], mean.bg = "gray", mean.cex = 1.2, link.col = "gray", link.lwd = 2, link.lty = 1)),
                    "procrustes shape" = plotRefToTarget(ref(), tar(), links = outline, mag = input$magnitude, method = input$method, gridPars = GP3)  
,
                    "PCA" = plotTangentSpace(embryo, groups = col.gp, axis1 = input$pcX, axis2 = input$pcY, legend = T)
,
                    "cell sizes" = #print(C2()$species), # print it out to the console for debugging.
                                      print(ggplot(C2(), aes(x=type, y=mean)) +
                                          geom_point(aes(color = factor(species))) +
                                          geom_errorbar(aes(ymin=low, ymax=high, color = factor(species)), width=0.1) +
                                          theme_classic()) ,

                    "trajectory analysis" = plot(TA(), pt.scale = 1.3, pt.seq.pattern = pts()) )
  })
  

  
 
}

shinyApp(ui = ui, server = server)
