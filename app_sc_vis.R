# Set working directory
wd <- dirname(rstudioapi::getActiveDocumentContext()$path) 

library(shiny)
library(Seurat)
library(ggplot2)
library(stringr)

default <- c("EPCAM,COL1A1,CD3E,MS4A1,CD79A,LYZ,MKI67")
height <- 5
width <- 5
dot_width <- 10
dot_height <- 5

# UI for Single Cell App
ui <- fluidPage(
  fluidRow(
    column(3,
           wellPanel(
             titlePanel("Datasets"),
             uiOutput("dataset_selector"),
             actionButton("update_dataset", "Update Dataset", width = "100%")
           ),
           wellPanel(
             titlePanel("Idents"),
             uiOutput("idents_panel")
           )
    ),
    column(9,
           fluidRow(
             column(6, wellPanel(plotOutput("umap_plot", height = "400px"))),
             column(6, wellPanel(plotOutput("feature_plot", height = "400px")))
           ),
           fluidRow(
             column(2,
                    tags$div(
                      HTML("<strong>Update UMAP:</strong>"),
                      actionButton("update_umap", "Update UMAP", width = "100%")
                    )
             ),
             column(2, textInput("filename_prefix", "Filename prefix:", value = "user", width = "100%")),
             column(2, numericInput("umap_ptsize", "UMAP Pt Size", value = 1.2, min = 0.1)),
             column(3,
                    tags$div(
                      HTML("<strong>Export Feature:</strong>"),
                      actionButton("export_feature", "Export Feature", width = "100%")
                    )
             ),
             column(3, textInput("feature_gene", "Gene to Plot:", value = "CD3E", width = "100%"))
           ),
           fluidRow(
             column(3,
                    tags$div(
                      HTML("<strong>Export UMAP & DotPlot:</strong>"),
                      actionButton("export_umap", "Export UMAP & DotPlot", width = "100%")
                    )
             ),
             column(9, textInput("dotplot_genes", "Genes for DotPlot (comma-separated):", value = "EPCAM,COL1A1,CD3E,MS4A1,CD79A,LYZ,MKI67", width = "100%"))
           ),
           fluidRow(column(12, wellPanel(plotOutput("dot_plot", height = "400px"))))
    )
  )
)

# Server for Single Cell App
server <- function(input, output, session) {
  
  # Reactive variable for Single Cell dataset
  seurat_data <- reactiveVal(NULL)
  
  # Dataset selector for Single Cell
  output$dataset_selector <- renderUI({
    rds_files <- list.files(path = "sc/", pattern = "*.rds", full.names = FALSE)
    selectInput("selected_dataset", "Choose a Single Cell dataset:", choices = rds_files, selected = NULL)
  })
  
  # Load dataset and update UI
  observeEvent(input$update_dataset, {
    req(input$selected_dataset)
    
    gc()  # Clean memory before loading new data
    
    data <- readRDS(file.path("sc/", input$selected_dataset))
    seurat_data(data)
    
    gc()  # Clean memory after loading new data
    
    # Unique Idents
    idents <- levels(Idents(seurat_data()))
    
    # Generate dynamic text boxes for Idents
    output$idents_panel <- renderUI({
      n <- length(idents)
      textboxes <- lapply(1:n, function(i) {
        fluidRow(
          column(2, strong(i)),
          column(10, textInput(paste0("ident_", i), label = NULL, value = idents[i]))
        )
      })
      do.call(tagList, textboxes)
    })
    
    # UMAP Plot
    output$umap_plot <- renderPlot({
      UMAPPlot(seurat_data(), label = TRUE, label.box = TRUE, pt.size = input$umap_ptsize) + NoLegend() + NoAxes()
    })
    
    # Feature Plot
    output$feature_plot <- renderPlot({
      feature_gene <- input$feature_gene
      if (!(feature_gene %in% rownames(seurat_data()))) {
        feature_gene <- "CD3E"
        ggtitle_text <- paste0(input$feature_gene, " is not found in the dataset, displaying default CD3E")
      } else {
        ggtitle_text <- feature_gene
      }
      FeaturePlot(seurat_data(), features = feature_gene, cols = c('grey', "red"), pt.size = input$umap_ptsize, order = TRUE, reduction='umap') + 
        ggtitle(ggtitle_text) + NoLegend() + NoAxes()
    })
    
    # DotPlot
    output$dot_plot <- renderPlot({
      genes <- str_split(input$dotplot_genes, ",")[[1]]
      genes <- str_trim(genes)
      genes <- genes[genes %in% rownames(seurat_data())]
      DotPlot(seurat_data(), features = genes, cols = c('grey', 'red')) + NoLegend() + RotatedAxis()
    })
  })
  
  # Update UMAP and DotPlot
  observeEvent(input$update_umap, {
    req(seurat_data())
    
    # Update Idents
    modified_idents <- sapply(1:length(unique(Idents(seurat_data()))), function(i) {
      input[[paste0("ident_", i)]]
    })
    
    names(modified_idents) <- levels(Idents(seurat_data()))
    seurat_data(RenameIdents(seurat_data(), modified_idents))
    
    # UMAP Plot
    output$umap_plot <- renderPlot({
      UMAPPlot(seurat_data(), label = TRUE, label.box = TRUE, pt.size = input$umap_ptsize) + NoLegend() + NoAxes()
    })
    
    # DotPlot
    output$dot_plot <- renderPlot({
      genes <- str_split(input$dotplot_genes, ",")[[1]]
      genes <- str_trim(genes)
      genes <- genes[genes %in% rownames(seurat_data())]
      DotPlot(seurat_data(), features = genes, cols = c('grey', 'red')) + NoLegend() + RotatedAxis()
    })
  })
  
  # Export UMAP & DotPlot
  observeEvent(input$export_umap, {
    req(seurat_data())
    current_time <- format(Sys.time(), "%Y%m%d_%H%M")
    dataset_name <- gsub("\\.rds$", "", input$selected_dataset)
    
    # Create output directory if missing
    if (!dir.exists(file.path(wd, "output"))) {
      dir.create(file.path(wd, "output"))
    }
    
    prefix <- input$filename_prefix
    idents_filename <- file.path(wd, "output", paste0(current_time, "_", dataset_name, "_", prefix, "_annotations.csv"))
    umap_filename <- file.path(wd, "output", paste0(current_time, "_", dataset_name, "_", prefix, "_UMAP.pdf"))
    dotplot_filename <- file.path(wd, "output", paste0(current_time, "_", dataset_name, "_", prefix, "_DotPlot.pdf"))
    
    modified_idents <- data.frame(cell = names(Idents(seurat_data())), annotation = Idents(seurat_data()))
    colnames(modified_idents)[2] <- paste0(prefix, "_annotations")
    write.csv(modified_idents, idents_filename, row.names = FALSE)
    
    # Export UMAP
    pdf(umap_filename, width = width, height = height)
    print(UMAPPlot(seurat_data(), label = TRUE, label.box = TRUE, pt.size = input$umap_ptsize) + NoLegend() + NoAxes())
    dev.off()
    
    # Export DotPlot
    pdf(dotplot_filename, width = dot_width, height = dot_height)
    genes <- str_split(input$dotplot_genes, ",")[[1]]
    genes <- str_trim(genes)
    genes <- genes[genes %in% rownames(seurat_data())]
    print(DotPlot(seurat_data(), features = genes, cols = c('grey', 'red')) + NoLegend() + RotatedAxis())
    dev.off()
  })
  
  # Export Feature Plot
  observeEvent(input$export_feature, {
    req(seurat_data())
    current_time <- format(Sys.time(), "%Y%m%d_%H%M")
    dataset_name <- gsub("\\.rds$", "", input$selected_dataset)
    
    # Create output directory if missing
    if (!dir.exists(file.path(wd, "output"))) {
      dir.create(file.path(wd, "output"))
    }
    
    prefix <- input$filename_prefix
    feature_filename <- file.path(wd, "output", paste0(current_time, "_", dataset_name, "_", prefix, "_", input$feature_gene, "_Feature.pdf"))
    
    pdf(feature_filename, width = width, height = height)
    print(
      FeaturePlot(seurat_data(), features = input$feature_gene, cols = c('grey', "red"), order = TRUE, pt.size = input$umap_ptsize, reduction='umap') + 
        ggtitle(ifelse(input$feature_gene %in% rownames(seurat_data()), input$feature_gene, paste0(input$feature_gene, " is not found in the dataset, displaying default CD3E"))) + 
        NoLegend() + NoAxes()
    )
    dev.off()
  })
}

shinyApp(ui, server)
