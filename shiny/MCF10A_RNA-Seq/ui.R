library(shiny)
library(markdown)
library(DT)

# Define UI for application that draws a histogram
shinyUI(
  fluidPage(
    
    # Application title
    titlePanel("MCF10A RNA-Seq analysis - Version 0.1"),
    
    mainPanel(
      tabsetPanel(type = "tabs",
                  tabPanel("Global DE - TGFb vs WT", DT::dataTableOutput("global_de_tgfb")),
                  tabPanel("Global DE - shZ KD vs WT", DT::dataTableOutput("global_de_shZ")),
                  tabPanel("Global Gene Expression  - TGFb vs WT", DT::dataTableOutput("global_expression_tgfb")),
                  tabPanel("Global Gene Expression  - shZ vs WT", DT::dataTableOutput("global_expression_shZ"))
      )
    )
  )
)