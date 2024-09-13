# Set working directory
wd <- dirname(rstudioapi::getActiveDocumentContext()$path) 

library(shiny)
library(Seurat)
library(ggplot2)
library(stringr)

default <- c("EPCAM,COL1A1,CD3E,MS4A1,CD79A,LYZ,MKI67")
height <- 5
width <- 5

# List of allowed metadata columns
allowed_metadata <- c("dapi_avg", "cd11c_avg", "lang_avg", "bcell_avg", "gland", "la_in", "la_out", "epi", "epi_dist", 
                      "sm_out", "wall", "wall_dist", "capsule", "capsule_dist", "percent.mt")

# UI for Spatial App
ui <- fluidPage(
  fluidRow(
    column(3,
           wellPanel(
             titlePanel("Datasets"),
             uiOutput("dataset_selector_spat"),
             actionButton("update_dataset_spat", "Update Dataset", width = "100%")
           ),
           wellPanel(
             titlePanel("Metadata Selector"),
             uiOutput("metadata_selector_spat")  # Dropdown for metadata selection
           ), 
           fluidRow(
             column(12, wellPanel(plotOutput("meta_plot", height = "600px")))
           ),
           fluidRow(
             column(6, textInput("pt_alpha3", "pt/alpha:", value = "1.6,1,1,1", width = "100%")),
             column(6, textInput("min_max3", "Min_Max:", value = "NA,NA", width = "100%"))
           )
             
    ),
    column(9,
           fluidRow(
             column(6, wellPanel(plotOutput("feature_plot_spat1", height = "600px"))),
             column(6, wellPanel(plotOutput("feature_plot_spat2", height = "600px")))
           ),
           fluidRow(
             column(3, textInput("filename_prefix_spat", "Filename prefix:", value = "user", width = "100%")),
             column(3, textInput("feature_gene_spat1", "Gene to Plot 1:", value = "CD3E", width = "100%")),
             column(3, textInput("feature_gene_spat2", "Gene to Plot 2:", value = "MS4A1", width = "100%")),
             column(3,
                    tags$div(
                      HTML("<strong>Export Feature Plots:</strong>"),
                      actionButton("export_feature_spat", "Export Feature Plots", width = "100%")
                    )
             )
           ),
           fluidRow(
             column(3, textInput("pt_alpha1", "pt/alpha:", value = "1.6,1,1,1", width = "100%")),
             column(3, textInput("min_max1", "Min_Max:", value = "NA,NA", width = "100%")),
             column(3, textInput("pt_alpha2", "pt/alpha:", value = "1.6,1,1,1", width = "100%")),
             column(3, textInput("min_max2", "Min_Max:", value = "NA,NA", width = "100%"))
           )
    )
  )
)

# Server for Spatial App
server <- function(input, output, session) {
  
  # Reactive variable for Spatial data
  seurat_data_spat <- reactiveVal(NULL)
  
  # Dataset selector for Spatial
  output$dataset_selector_spat <- renderUI({
    rds_files <- list.files(path = "vis/", pattern = "*.rds", full.names = FALSE)
    selectInput("selected_dataset_spat", "Choose a Spatial dataset:", choices = rds_files, selected = NULL)
  })
  
  # Load Spatial dataset and update UI
  observeEvent(input$update_dataset_spat, {
    req(input$selected_dataset_spat)
    
    gc()  # Clean memory before loading new data
    
    data <- readRDS(file.path("vis/", input$selected_dataset_spat))
    seurat_data_spat(data)
    
    gc()  # Clean memory after loading new data
    
    # Extract metadata columns that match the allowed list
    metadata_columns <- colnames(seurat_data_spat()@meta.data)
    available_metadata <- metadata_columns[metadata_columns %in% allowed_metadata]
    
    # If there are matching metadata columns, generate dropdown, else display a message
    if (length(available_metadata) > 0) {
      output$metadata_selector_spat <- renderUI({
        selectInput("selected_metadata_spat", "Choose a Metadata Column:", choices = available_metadata, selected = available_metadata[1])
      })
    } else {
      output$metadata_selector_spat <- renderUI({
        tags$p("No relevant metadata columns found in this dataset.")
      })
    }
    
    # First SpatialFeaturePlot
    output$feature_plot_spat1 <- renderPlot({
      feature_gene <- input$feature_gene_spat1
      if (!(feature_gene %in% rownames(seurat_data_spat()))) {
        feature_gene <- "CD3E"
        ggtitle_text <- paste0(input$feature_gene_spat1, " is not found in the dataset, displaying default CD3E")
      } else {
        ggtitle_text <- feature_gene
      }
      SpatialFeaturePlot(seurat_data_spat(), features = input$feature_gene_spat1, 
                         alpha = c(as.numeric(str_split(input$pt_alpha1, "\\,")[[1]][2]), as.numeric(str_split(input$pt_alpha1, "\\,")[[1]][3])), 
                         pt.size = as.numeric(str_split(input$pt_alpha1, "\\,")[[1]][1]), 
                         min.cutoff = if(is.na(str_split(input$min_max1, "\\,")[[1]][1])){NA}else{as.numeric(str_split(input$min_max1, "\\,")[[1]][1])}, 
                         max.cutoff = if(is.na(str_split(input$min_max1, "\\,")[[1]][2])){NA}else{as.numeric(str_split(input$min_max1, "\\,")[[1]][2])},
                         image.alpha = as.numeric(str_split(input$pt_alpha1, "\\,")[[1]][4])
      ) + 
        ggtitle(ggtitle_text) + NoLegend() + NoAxes()+RestoreLegend(position = "right")
    })
    
    # Second SpatialFeaturePlot
    output$feature_plot_spat2 <- renderPlot({
      feature_gene <- input$feature_gene_spat2
      if (!(feature_gene %in% rownames(seurat_data_spat()))) {
        feature_gene <- "MS4A1"
        ggtitle_text <- paste0(input$feature_gene_spat2, " is not found in the dataset, displaying default MS4A1")
      } else {
        ggtitle_text <- feature_gene
      }
      SpatialFeaturePlot(seurat_data_spat(), features = feature_gene, alpha = c(as.numeric(str_split(input$pt_alpha2, "\\,")[[1]][2]), as.numeric(str_split(input$pt_alpha2, "\\,")[[1]][3])), 
                         pt.size = as.numeric(str_split(input$pt_alpha2, "\\,")[[1]][1]),                         
                         min.cutoff = if(is.na(str_split(input$min_max2, "\\,")[[1]][1])){NA}else{as.numeric(str_split(input$min_max2, "\\,")[[1]][1])}, 
                         max.cutoff = if(is.na(str_split(input$min_max2, "\\,")[[1]][2])){NA}else{as.numeric(str_split(input$min_max2, "\\,")[[1]][2])},
                         image.alpha = as.numeric(str_split(input$pt_alpha2, "\\,")[[1]][4])
      ) + 
        ggtitle(ggtitle_text) + NoLegend() + NoAxes()+RestoreLegend(position = "right")
    })
  
  output$meta_plot <- renderPlot({
    feature_gene <- input$selected_metadata_spat
    ggtitle_text <- feature_gene
    
    SpatialFeaturePlot(seurat_data_spat(), features = feature_gene, alpha = c(as.numeric(str_split(input$pt_alpha3, "\\,")[[1]][2]), as.numeric(str_split(input$pt_alpha3, "\\,")[[1]][3])), 
                       pt.size = as.numeric(str_split(input$pt_alpha3, "\\,")[[1]][1]),                         
                       min.cutoff = if(is.na(str_split(input$min_max3, "\\,")[[1]][1])){NA}else{as.numeric(str_split(input$min_max3, "\\,")[[1]][1])}, 
                       max.cutoff = if(is.na(str_split(input$min_max3, "\\,")[[1]][2])){NA}else{as.numeric(str_split(input$min_max3, "\\,")[[1]][2])},
                       image.alpha = as.numeric(str_split(input$pt_alpha3, "\\,")[[1]][4])
    ) + 
      ggtitle(ggtitle_text) + NoLegend() + NoAxes()+RestoreLegend(position = "right")
  })
})
  
  # Export Feature Plots
  observeEvent(input$export_feature_spat, {
    req(seurat_data_spat())
    
    current_time <- format(Sys.time(), "%Y%m%d_%H%M")
    dataset_name_spat <- gsub("\\.rds$", "", input$selected_dataset_spat)
    
    # Create output directory if missing
    if (!dir.exists(file.path(wd, "output"))) {
      dir.create(file.path(wd, "output"))
    }
    
    prefix_spat <- input$filename_prefix_spat
    feature_filename_spat1 <- file.path(wd, "output", paste0(current_time, "_", dataset_name_spat, "_", prefix_spat, "_", input$feature_gene_spat1, "_Feature.pdf"))
    feature_filename_spat2 <- file.path(wd, "output", paste0(current_time, "_", dataset_name_spat, "_", prefix_spat, "_", input$feature_gene_spat2, "_Feature.pdf"))
    feature_filename_spat3 <- file.path(wd, "output", paste0(current_time, "_", dataset_name_spat, "_", prefix_spat, "_", input$selected_metadata_spat, "_Meta.pdf"))
    
    # Export first FeaturePlot
    pdf(feature_filename_spat1, width = width, height = height)
    print(
      SpatialFeaturePlot(seurat_data_spat(), features = input$feature_gene_spat1, 
                         alpha = c(as.numeric(str_split(input$pt_alpha1, "\\,")[[1]][2]), as.numeric(str_split(input$pt_alpha1, "\\,")[[1]][3])), 
                         pt.size = as.numeric(str_split(input$pt_alpha1, "\\,")[[1]][1]),                         
                         min.cutoff = if(is.na(str_split(input$min_max1, "\\,")[[1]][1])){NA}else{as.numeric(str_split(input$min_max1, "\\,")[[1]][1])}, 
                         max.cutoff = if(is.na(str_split(input$min_max1, "\\,")[[1]][2])){NA}else{as.numeric(str_split(input$min_max1, "\\,")[[1]][2])},
                         image.alpha = as.numeric(str_split(input$pt_alpha1, "\\,")[[1]][4]))+
        NoLegend() + NoAxes()+RestoreLegend(position = "right")
    )
    dev.off()
    
    # Export second FeaturePlot
    pdf(feature_filename_spat2, width = width, height = height)
    print(
      SpatialFeaturePlot(seurat_data_spat(), features = input$feature_gene_spat2, 
                         alpha = c(as.numeric(str_split(input$pt_alpha2, "\\,")[[1]][2]), as.numeric(str_split(input$pt_alpha2, "\\,")[[1]][3])), 
                         pt.size = as.numeric(str_split(input$pt_alpha2, "\\,")[[1]][1]),                         
                         min.cutoff = if(is.na(str_split(input$min_max2, "\\,")[[1]][1])){NA}else{as.numeric(str_split(input$min_max2, "\\,")[[1]][1])}, 
                         max.cutoff = if(is.na(str_split(input$min_max2, "\\,")[[1]][2])){NA}else{as.numeric(str_split(input$min_max2, "\\,")[[1]][2])},
                         image.alpha = as.numeric(str_split(input$pt_alpha2, "\\,")[[1]][4]))+
        NoLegend() + NoAxes()+RestoreLegend(position = "right")
    )
    dev.off()
    # Export Mera
    pdf(feature_filename_spat3, width = width, height = height)
    print(
      SpatialFeaturePlot(seurat_data_spat(), features = input$selected_metadata_spat, 
                         alpha = c(as.numeric(str_split(input$pt_alpha3, "\\,")[[1]][2]), 
                                   as.numeric(str_split(input$pt_alpha3, "\\,")[[1]][3])), 
                         pt.size = as.numeric(str_split(input$pt_alpha3, "\\,")[[1]][1]),                         
                         min.cutoff = if(is.na(str_split(input$min_max3, "\\,")[[1]][1])){NA}else{as.numeric(str_split(input$min_max3, "\\,")[[1]][1])}, 
                         max.cutoff = if(is.na(str_split(input$min_max3, "\\,")[[1]][2])){NA}else{as.numeric(str_split(input$min_max3, "\\,")[[1]][2])},
                         image.alpha = as.numeric(str_split(input$pt_alpha3, "\\,")[[1]][4])
      ) + NoLegend() + NoAxes()+RestoreLegend(position = "right")
    )
    dev.off()
  })
}

shinyApp(ui, server)

