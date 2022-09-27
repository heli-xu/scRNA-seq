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
      height = 60
    ),
    div(
      style = "padding-top: 50px;",
      spin_loader()) ,
    div("Thanks for your patience! scRNA datasets are large and take a few seconds to load!", id = "intro-note")
    
  )
}


load("R/all_features_humanHF.rdata")

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
  div(  titlePanel("Heart scRNA-Seq Dashboard: human heart failure" ),id = 'navbar'),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      selectInput("dataset",
                  "Dataset:",
                  choices = "2020 NatCellBio human HF",
                  selected = "2020 NatCellBio human HF"),
      selectInput("gene",
                  "Gene name:",
                  choices = all_features_humanHF,
                  selected = "MYH6"),
      selectInput("cell",
                  "Cell type:",
                  choices = c("CM", "EC", "MP"),
                  selected = "CM"),
      
      # Button for download data
      downloadButton("downloadData", "Export dotplot data")
      
    ),
    mainPanel(
      width = 9,
      column(6, plotOutput("DotPlot")),
      column(6, plotOutput("ViolinPlot")),
      a(href = "https://doi.org/10.1038/s41556-019-0446-7", 
        "Publication Source")
    )
  ),
  
  ## Footer
  includeHTML("HTML/Footer.html")
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  
  ## Load data then remove loading screen
  ## Load data (takes long time)
  load("R/server_data_humanHF.rdata")
  
  waiter_hide()
  
  
  output$DotPlot <- renderPlot({
    object_to_plot = server_data_humanHF %>%
      filter(cell == input$cell) %>%
      pull(obj)
    
    
    DotPlot(object_to_plot[[1]],
            features = input$gene,
            col.min = 0, col.max = 3,
            dot.min = 0, dot.scale = 6)+
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
    
    
    
  })
  
  output$ViolinPlot <- renderPlot({
    object_to_plot = server_data_humanHF %>%
      filter(cell == input$cell) %>%
      pull(obj)

    VlnPlot(object_to_plot[[1]], features = input$gene, pt.size = 0)
    
  })
  
  # Downloadable csv of dotplot data
  output$downloadData <- downloadHandler(
    filename= function() {
      paste('humanHF', ' ', input$cell, ' ', input$gene, '.csv', sep = '')
    },
    content = function(file) {
      object_to_plot = server_data_humanHF %>% 
        filter(cell == input$cell) %>% 
        pull(obj)
      
      plot <- DotPlot(object_to_plot[[1]],
                      features = input$gene)
      
      write.csv(plot$data, file)
    }
  )
}

# Run the application 
shinyApp(ui = ui, server = server)
