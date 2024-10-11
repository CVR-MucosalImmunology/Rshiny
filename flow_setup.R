ui <- fluidPage(
  column(3,
         wellPanel(
           titlePanel("Select Sample"),
           uiOutput("sample_selector"),
           actionButton("load_sample", "Update sample", width = "100%")
         ),
         wellPanel(
           titlePanel("Select Channel"),
           uiOutput("channel_selector")
         ),
         wellPanel(
           titlePanel("Select Y"),
           uiOutput("Y_selector"),
           actionButton("load_channel", "Update axes", width = "100%")       ,    
           actionButton("load_plot", "Update plot", width = "100%")

           )
  ),
  column(5,
         wellPanel(plotOutput("facs_plot", height = "600px"))
  ),
  column(4,
         fluidRow(
           column(6,
             tags$div(
               HTML("<strong>Phenotypic marker:</strong>"),
               checkboxInput("inc_pheno", "", value=T, width='100%'))),
           column(6,
             numericInput("cofactor", "Choose arcsinh cofactor:", value=50, width='100%'))
         ),
         fluidRow(
           column(6,
             tags$div(
               HTML("<strong>Cluster marker:</strong>"),
               checkboxInput("inc_clust", "", value=F, width='100%'))),
           column(6, 
                  numericInput("quant", "Filter Bottom quantile (%):", value=0.3, width='100%'))
                  
         ),         
         fluidRow(
           tags$div(
             HTML("<strong>Add data:</strong>"),
             actionButton("add", "Add Data", width='100%'))
           
         ),

         fluidRow(
           tags$div(
             HTML("<strong>Updated data table:</strong>")
         )), 
         fluidRow(
           DTOutput("csv_table")
         ),
         fluidRow(
           tags$div(
             HTML("<strong>Save data:</strong>"),
             actionButton("save", "Save", width='100%'))
           
         ),
  )  
);server <- function(input, output, session) {
  data <- reactiveVal(NULL)
  raw <- reactiveVal(NULL)
  trans <- reactiveVal(NULL)
  csv <- reactiveVal(NULL)
  temp <- reactiveVal(NULL)

  # Load sample -------------------------------------------------------------
  output$sample_selector <- renderUI({
    query_files <- list.files(path = "flow/", pattern = "*.csv", full.names = FALSE)
    selectInput("sample_data", "Choose flow sample data:", choices = query_files, selected = NULL)
  })
  
  observeEvent(input$load_sample, {
    req(input$sample_data)
    
    r <- read.csv(paste0("flow/", input$sample_data))
    coln <- colnames(r)
    df <- data.frame(pheno = rep(0, length(coln)), clust = rep(0, length(coln)), cofactor = NA, filt=NA,row.names = coln)
    
    
    csv(df)
    raw(r)
    
    updateSelectInput(session, "load_channel", choices = coln, selected = coln[7])
    updateSelectInput(session, "Y_selector", choices = coln, selected = coln[9])
  })
  
  output$channel_selector <- renderUI({
    req(raw())
    selectInput("channel_selector", "Select Channel to transform:", choices = rownames(csv()), selected = NULL)
  })
  
  output$Y_selector <- renderUI({
    req(raw())
    selectInput("Y_selector", "Select Y axis:", choices = rownames(csv()), selected = NULL)
  })
  output$csv_table <- renderDT({
    req(csv())
    datatable(csv(), options = list(pageLength = 10, autoWidth = TRUE))
  })
  # Load plot ---------------------------------------------------------------
  observeEvent(input$load_channel, {
    req(input$channel_selector, input$Y_selector, raw(), input$quant)
    # Apply the transformation
    r <- raw()
    if(input$quant>0){
      percent = input$quant/100
      lower_cut_off <- quantile(as.numeric(r[, input$channel_selector]), probs = percent, na.rm = TRUE)
      upper_cut_off <- quantile(as.numeric(r[, input$channel_selector]), probs = 0.999, na.rm = TRUE)
      filtered_r <- r[r[, input$channel_selector] > lower_cut_off & r[, input$channel_selector] < upper_cut_off, ]

    } else {
      filtered_r <- r
    }

    filtered_r <- as.data.table(filtered_r)
    
    transformed <- do.asinh(filtered_r, c(input$channel_selector, input$Y_selector), cofactor = input$cofactor)

    temp(transformed)
  })
  
observeEvent(input$load_plot, {
    output$facs_plot <- renderPlot({
      ggplot(temp(),aes_string(paste0(input$channel_selector, "_asinh"), paste0(input$Y_selector, "_asinh"))) +
        geom_bin_2d(bins = 200) +
        scale_fill_gradientn(colors = c("blue", "green", "yellow", "red")) +
        theme_bw()
    })
  })

observeEvent(input$add,{
  df <- csv()
  df[rownames(df)==input$channel_selector,1] = input$inc_pheno
  df[rownames(df)==input$channel_selector,2] = input$inc_clust
  df[rownames(df)==input$channel_selector,3] = input$cofactor
  df[rownames(df)==input$channel_selector,4] = input$quant
  csv(df)
})

observeEvent(input$save, {
  write.csv(csv(), paste0(wd,"/meta/", substr(input$sample_data,1,nchar(input$sample_data)-4), "_meta.csv"))
})  

}
shinyApp(ui,server);rm(ui,server)
