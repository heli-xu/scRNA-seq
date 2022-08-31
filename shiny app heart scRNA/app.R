{
  ## Dependencies
  library(shiny)
  library(dplyr)
  library(Seurat)
  library(ggplot2)
  library(waiter)

  
  ## Loader code: Global Scope
  loading_screen = div(
    tags$img(
      src = "logo.png",
      height = 275
    ),
    div(
      # style = "padding-top: 50px;",
      spin_loader()) )
}


load("data/mouse MI/all_features_MI.rdata")

# Define UI for application that draws a histogram
ui <- fluidPage(
  ## Loader code: UI (start)
  useWaiter(),
  waiterShowOnLoad(html = loading_screen,
                   color = 'white'),
  ## HTML Head
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "style.css")
  ),
  
  
  
  
  # Application title
  # h1('test', id = 'navbar'),
  div(  titlePanel("heart scRNA" ),id = 'navbar'),

  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      width = 3,
      selectInput("dataset",
                  "Dataset:",
                  choices = "2019 eLife mouse MI"),
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
      width = 9,
      column(6, plotOutput("DotPlot")),
      column(6, plotOutput("ViolinPlot"))
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  ## Load data then remove loading screen
  ## Load data (takes long time)
  load("data/mouse MI/MI_EC_normalized.rdata") #object name: TIP_MI_EC
  load("data/mouse MI/MI_MP_normalized.rdata") #TIP_MI_MP
  waiter_hide()
  
  
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
    
    VlnPlot(object_to_plot, features = input$gene, pt.size = 0)
    
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
