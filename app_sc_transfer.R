# Set working directory
wd <- getwd()

library(shiny)
library(Seurat)
library(ggplot2)
library(stringr)
library(ComplexHeatmap)

# UI for Single Cell App
ui <- fluidPage(
  fluidRow(
    column(3,
           wellPanel(
             titlePanel("Query Dataset"),
             uiOutput("dataset_selector")
           ),
           wellPanel(
             titlePanel("Reference Dataset"),
             uiOutput("dataset_selector2"),
             actionButton("load_datasets", "Update Dataset", width = "100%")
             
           ), 
           wellPanel(
             titlePanel("Group by (query):"),
             uiOutput("groupby")
             ),
           wellPanel(
             titlePanel("Group by (ref):"),
             uiOutput("groupby2"),
             actionButton("update_dataset", "Update Idents", width = "100%")
           )
    ),
    column(9,
           fluidRow(
             column(4, wellPanel(plotOutput("umap_plot", height = "400px"))),
             column(4, wellPanel(plotOutput("umap_plot2", height = "400px"))),
             column(4, wellPanel(plotOutput("umap_plot3", height = "400px")))
           ),
           fluidRow(             
             column(2, textInput("filename_prefix", "Filename prefix:", value = "user", width = "100%")),
             column(2, numericInput("umap_ptsize", "UMAP Pt Size", value = 1.2, min = 0.1, width='100%')),
             column(2, 
                    tags$div(
                      HTML("<strong>Include legend:</strong>"),
                      checkboxInput("legend", "", value=F, width='100%'))
                    ),
             column(2, 
                    tags$div(
                      HTML("<strong>Include labels:</strong>"),
                      checkboxInput("label", "", value=T, width='100%'))
             ),             
             column(4, 
                    tags$div(
                      HTML("<strong>Transfer Data:</strong>"),
                      actionButton("transfer", "Transfer", width = "100%"))
             ),
           ),
           fluidRow(
             column(6, wellPanel(plotOutput("heatmap", height = "400px"))),
             column(6,
                    tags$div(
                      HTML("<strong>Export data & plots:</strong>"),
                      actionButton("export_all", "Export UMAP & DotPlot", width = "100%"))
                    )
             )
           )
  )
)
server <- function(input, output, session) {
  
  # Reactive values for query and reference Seurat datasets
  query_data <- reactiveVal(NULL)
  reference_data <- reactiveVal(NULL)
  que_ref <- reactiveVal(NULL)
  heatmap <- reactiveVal(NULL)
  
  # Dataset selectors for query and reference datasets
  output$dataset_selector <- renderUI({
    query_files <- list.files(path = "sc/", pattern = "*.rds", full.names = FALSE)
    selectInput("selected_query_dataset", "Choose a Query dataset:", choices = query_files, selected = NULL)
  })
  
  output$dataset_selector2 <- renderUI({
    reference_files <- list.files(path = "sc/", pattern = "*.rds", full.names = FALSE)
    selectInput("selected_reference_dataset", "Choose a Reference dataset:", choices = reference_files, selected = NULL)
  })
  
  # Load datasets when the 'Load Dataset' button is clicked
  observeEvent(input$load_datasets, {
    req(input$selected_query_dataset, input$selected_reference_dataset)
    
    # Load query dataset
    query <- readRDS(file.path("sc/", input$selected_query_dataset))
    query$default <- Idents(query)
    
    query_data(query)
    
    # Load reference dataset
    reference <- readRDS(file.path("sc/", input$selected_reference_dataset))
    reference$default <- Idents(reference)
    reference_data(reference)

    # Output the metadata columns for the Group by dropdowns
    updateSelectInput(session, "groupby", choices = colnames(query@meta.data), selected = 'default')
    updateSelectInput(session, "groupby2", choices = colnames(reference@meta.data), selected = 'default')
  })
  
  # UI for selecting metadata columns for group.by in query and reference datasets
  output$groupby <- renderUI({
    req(query_data())
    selectInput("groupby", "Select Query Metadata Column:", choices = colnames(query_data()@meta.data), selected = NULL)
  })
  
  output$groupby2 <- renderUI({
    req(reference_data())
    selectInput("groupby2", "Select Reference Metadata Column:", choices = colnames(reference_data()@meta.data), selected = NULL)
  })
  
  # Render UMAP for query dataset
  output$umap_plot <- renderPlot({
    req(query_data())
    if(input$legend) {
      DimPlot(query_data(), reduction = "umap", group.by = input$groupby, pt.size = input$umap_ptsize, label=input$label)+NoAxes()+ggtitle("Query")
    } else {
      DimPlot(query_data(), reduction = "umap", group.by = input$groupby, pt.size = input$umap_ptsize, label=input$label)+NoAxes()+ggtitle("Query")+NoLegend()
    }  
    })
  
  # Render UMAP for reference dataset
  output$umap_plot2 <- renderPlot({
    req(reference_data())
    if(input$legend){
      DimPlot(reference_data(), reduction = "umap", group.by = input$groupby2, pt.size = input$umap_ptsize, label=input$label)+NoAxes()+ggtitle("Reference")
    }else {
      DimPlot(reference_data(), reduction = "umap", group.by = input$groupby2, pt.size = input$umap_ptsize, label=input$label)+NoAxes()+ggtitle("Reference")+NoLegend()
    }  
    })
  
  # Update UMAPs with selected metadata groupings when 'Update Idents' button is clicked
  observeEvent(input$update_dataset, {
    req(query_data(), reference_data(), input$groupby, input$groupby2)
    
    # Update query UMAP with the selected metadata column
    output$umap_plot <- renderPlot({
      if(input$legend) {
        DimPlot(query_data(), reduction = "umap", group.by = input$groupby, pt.size = input$umap_ptsize, label=input$label)+NoAxes()+ggtitle("Query")
      } else {
        DimPlot(query_data(), reduction = "umap", group.by = input$groupby, pt.size = input$umap_ptsize, label=input$label)+NoAxes()+ggtitle("Query")+NoLegend()
      }
    })
    
    # Update reference UMAP with the selected metadata column
    output$umap_plot2 <- renderPlot({
      if(input$legend){
        DimPlot(reference_data(), reduction = "umap", group.by = input$groupby2, pt.size = input$umap_ptsize, label=input$label)+NoAxes()+ggtitle("Reference")
      }else {
        DimPlot(reference_data(), reduction = "umap", group.by = input$groupby2, pt.size = input$umap_ptsize, label=input$label)+NoAxes()+ggtitle("Reference")+NoLegend()
      }
    })
  })
  observeEvent(input$transfer, {
    req(query_data(), reference_data(), input$groupby, input$groupby2)
    que <- query_data()
    ref <- reference_data()
    predictions <- TransferData(
      anchorset = FindTransferAnchors(reference = ref, query = que, features=VariableFeatures(ref)),
      refdata = FetchData(ref, input$groupby2)[,1]
    )
    que <- AddMetaData(que, metadata = predictions$predicted.id, col.name="predicted.id")
    
    predictions$orig =FetchData(que, input$groupby)[,1]
    df <- as.data.frame(matrix(data=NA,ncol=length(unique(predictions$orig)), nrow=length(unique(predictions$predicted.id))))
    colnames(df) = unique(predictions$orig); 
    rownames(df) = paste0("prediction.score.",gsub(" ", ".",unique(predictions$predicted.id)))
    df2 <- as.data.frame(matrix(data=NA,ncol=length(unique(predictions$orig)), nrow=length(unique(predictions$predicted.id))))
    colnames(df2) = unique(predictions$orig); 
    rownames(df2) = paste0("prediction.score.",gsub(" ", ".",unique(predictions$predicted.id)))
    
    for(col in 1:ncol(df)) {
      for(row in 1:nrow(df)){
        x = mean(predictions[predictions$orig==colnames(df)[col], rownames(df)[row]])
        x2 = mean(predictions[predictions$orig==colnames(df)[col], 
                              rownames(df)[row]]/
                    predictions[predictions$orig==colnames(df)[col], 
                                "prediction.score.max"])
        
        if(is.na(x) | is.infinite(x)){
          x=0
        }
        if(is.na(x2) | is.infinite(x2)){
          x2=0
        }
        df[row,col] <-x 
        df2[row,col] <-x2
      }
    }
    que_ref(que)
    
    # Update query UMAP with the selected metadata column
    output$umap_plot3 <- renderPlot({
      if(input$legend){
        DimPlot(que, reduction = "umap", group.by ="predicted.id", pt.size = input$umap_ptsize, label=input$label)+ggtitle("Transfered labels")+NoAxes()
      }else{
        DimPlot(que, reduction = "umap", group.by ="predicted.id", pt.size = input$umap_ptsize, label=input$label)+ggtitle("Transfered labels")+NoLegend()+NoAxes()
      }
    })
    
    heatmap(df)
    # Update reference UMAP with the selected metadata column
    output$heatmap <- renderPlot({
      ComplexHeatmap::Heatmap(na.omit(df))
      })
    
  })
  
  # Export Feature Plot
  observeEvent(input$export_all, {
    req(query_data(),reference_data(),que_ref())
    current_time <- format(Sys.time(), "%Y%m%d_%H%M")
    dataset_name <- gsub("\\.rds$", "", input$selected_query_dataset)
    dataset_name2 <- gsub("\\.rds$", "", input$selected_reference_dataset)
    
    # Create output directory if missing
    if (!dir.exists(file.path(wd, "output"))) {
      dir.create(file.path(wd, "output"))
    }
    
    prefix <- input$filename_prefix
    
    umap_filename <- file.path(wd, "output", paste0(current_time, "_", dataset_name2, "_to_",dataset_name,"_" , prefix, "_UMAPS.pdf"))
    heatmap_filename <- file.path(wd, "output", paste0(current_time, "_", dataset_name2, "_to_",dataset_name,"_" , prefix, "_Heatmap.pdf"))
    metadata_filename <- file.path(wd, "output", paste0(current_time, "_", dataset_name, "_from_",dataset_name2,"_" , prefix, "_metadata.csv"))
    heatmapcsv_filename <- file.path(wd, "output", paste0(current_time, "_", dataset_name, "_from_",dataset_name2,"_" , prefix, "_heatmapCSV.csv"))
    
    pdf(umap_filename, width = width, height = height)
    if(input$legend){
      print(
        DimPlot(query_data(), reduction = "umap", group.by =input$groupby, pt.size = input$umap_ptsize, label=input$label)+ggtitle("Query")+NoAxes(),
      )
      print(
          DimPlot(reference_data(), reduction = "umap", group.by =input$groupby2, pt.size = input$umap_ptsize, label=input$label)+ggtitle("Reference")+NoAxes(),
      )
      print(
        DimPlot(que_ref(), reduction = "umap", group.by ="predicted.id", pt.size = input$umap_ptsize, label=input$label)+ggtitle("Transfered labels")+NoAxes()
      )
    }else{      
      print(
        DimPlot(query_data(), reduction = "umap", group.by =input$groupby, pt.size = input$umap_ptsize, label=input$label)+ggtitle("Query")+NoAxes()+NoLegend(),
        )
      print(
        DimPlot(reference_data(), reduction = "umap", group.by =input$groupby2, pt.size = input$umap_ptsize, label=input$label)+ggtitle("Reference")+NoAxes()+NoLegend(),
        )
      print(
        DimPlot(que_ref(), reduction = "umap", group.by ="predicted.id", pt.size = input$umap_ptsize, label=input$label)+ggtitle("Transfered labels")+NoAxes()+NoLegend()      
        )
    }
    dev.off()
    #export heatmap
    pdf(heatmap_filename, width = widthheat, height = heightheat)
    print(ComplexHeatmap::Heatmap(na.omit(heatmap())))
    dev.off()
    
    write.csv(FetchData(que_ref(), "predicted.id"), metadata_filename)
    write.csv(heatmap(), heatmapcsv_filename)
  })
}


runApp(shinyApp(ui, server))
