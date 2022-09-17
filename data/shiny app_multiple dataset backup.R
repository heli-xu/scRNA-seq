{
  ## Dependencies
  library(shiny)
  library(dplyr)
  library(Seurat)
  library(ggplot2)
  library(waiter)
  library(tidyverse)
  
  
  ## Loader code: Global Scope
  loading_screen = div(
    tags$img(
      src = "logo.png",
      height = 60
    ),
    div(
      style = "padding-top: 50px;",
      spin_loader()) ,
    div("Thanks for your patience! scRNA datasets are large and take a few seconds to load!", id = "intro-note")
    
  )
}


load("R/ui_features.rdata")

# Define UI for application that draws a histogram
ui <- fluidPage(
  ## Loader code: UI (start)
  useWaiter(),
  waiterShowOnLoad(html = loading_screen,
                   color = 'white'),
  ## HTML Head
  tags$head(includeCSS("CSS/styles.css")),
  tags$head(includeCSS("CSS/Header.css")),
  
  ## Header
  includeHTML("HTML/Header.html"),
  
  # App 
  div(  titlePanel("Heart scRNA Exploratory Data Tool: mouse MI" ),id = 'navbar'),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      selectInput("dataset",
                  "Dataset:",
                  choices = c("2019 eLife mouse MI", 
                              "2020 Circulation mouse TAC"),
                  selected = "2019 eLife mouse MI"),
      selectInput("gene",
                  "Gene name:",
                  choices = " ",
                  selected = "Pecam1"),
      selectInput("cell",
                  "Cell type:",
                  choices = c("EC", "MP"),
                  selected = "EC")
    ),
    mainPanel(
      width = 9,
      column(6, plotOutput("DotPlot")),
      column(6, plotOutput("ViolinPlot")),
      a(href = "https://quarto.org/docs/authoring/article-layout.html#body-column", 
        "Publication Source")
    )
  ),
  
  ## Footer
  includeHTML("HTML/Footer.html")
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  
  observeEvent(input$dataset,
               updateSelectInput(session, "gene", "Gene name:",
                                 choices = ui_features$gene[ui_features$dataset == input$dataset],
                                 selected = "Pecam1")
  )
  
  
  ## Load data then remove loading screen
  ## Load data (takes long time)
  load("R/server_data.rdata")
  
  waiter_hide()
  
  
  output$DotPlot <- renderPlot({
    object_to_plot = server_data %>% 
      filter(dataset == input$dataset,
             cell == input$cell) %>% 
      pull(obj)
    
    
    DotPlot(object_to_plot[[1]],
            features = input$gene,
            col.min = 0, col.max = 3,
            dot.min = 0, dot.scale = 6)+
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
    
    
    
  })
  
  output$ViolinPlot <- renderPlot({
    object_to_plot = server_data %>% 
      filter(dataset == input$dataset,
             cell == input$cell) %>% 
      pull(obj)
    
    VlnPlot(object_to_plot[[1]], features = input$gene, pt.size = 0)
    
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
