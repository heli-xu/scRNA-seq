#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(dplyr)
library(Seurat)
library(tidyverse)

load("data/mouse MI/all_features_MI.rdata")

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("heart scRNA"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            selectInput("gene",
                        "Gene name:",
                        choices = all_features,
                        selected = "Pecam1"),
            selectInput("cell",
                        "Cell type:",
                        choices = c("EC", "MP"))
        ),

        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("DotPlot"),
           plotOutput("ViolinPlot")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

  #object name: TIP_MI_EC
  load("data/mouse MI/MI_EC_normalized.rdata")
  #TIP_MI_MP
  load("data/mouse MI/MI_MP_normalized.rdata")
  
      output$DotPlot <- renderPlot({
        # select object based on input
        object_to_plot = NULL
        if (input$cell == "EC"){object_to_plot = TIP_MI_EC} else 
          if  (input$cell == "MP"){
            object_to_plot = TIP_MI_MP
          }
        
        DotPlot(object_to_plot,
                features = input$gene,
                col.min = 0, col.max = 3,
                dot.min = 0, dot.scale = 6)+
          theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
        
      

    })
      
      output$ViolinPlot <- renderPlot({
        # select object based on input
        object_to_plot = NULL
        if (input$cell == "EC"){object_to_plot = TIP_MI_EC} else 
          if  (input$cell == "MP"){
            object_to_plot = TIP_MI_MP
          }
        
        VlnPlot(object_to_plot, features = input$gene)
        
      })
}

# Run the application 
shinyApp(ui = ui, server = server)
