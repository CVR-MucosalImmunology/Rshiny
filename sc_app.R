wd <- getwd()

library(shiny)
library(Seurat)
library(ggplot2)
library(stringr)
library(dplyr)

ui <- fluidPage(
  titlePanel("Single Cell RNA Analysis"),
  
  tabsetPanel(

# Single Cell Vis ---------------------------------------------------------
    tabPanel("scRNA Vis",
             fluidPage(
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
                          column(2, numericInput("ptsize", "UMAP Pt Size", value = 1.2, min = 0.1)),
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

# DE Analysis -------------------------------------------------------------

    ), tabPanel("DE Analysis",
                fluidPage(
                  fluidRow(
                    column(3,
                           wellPanel(
                             titlePanel("Ident 1:"),
                             uiOutput("ident1")
                           ),
                           wellPanel(
                             titlePanel("Ident 2:"),
                             uiOutput("ident2"),
                           )
                    ),
                    column(9,
                           fluidRow(
                             column(6, wellPanel(plotOutput("umap_plot2", height = "400px"))),
                             column(6, wellPanel(plotOutput("volcano", height = "400px")))
                           ),
                           fluidRow(     
                             column(2, 
                                    tags$div(
                                      HTML("<strong>Run DE:</strong>"),
                                      actionButton("DE", "DE", width = "100%"))
                             ),
                             column(2, numericInput("num", "# Genes", value = 30, min = 1, width='100%')),
                             column(2, numericInput("signif", "pval filter", value = 0.001, width='100%')),

                             column(1, numericInput("logfc", "logFC", value = 0.1, min = 0.01, width='100%')),
                             column(2, numericInput("min_pct", "Minimum % expression", value = 0.1, min = 0.01, width='100%')),
                             column(1, numericInput("min_diff_pct", "Minimum diff %", value = 0.1, min = 0.01, width='100%')),
                             column(1, numericInput("maxoverlap", "Label overlap", value = 10, min = 1, width='100%')),
                             
                             column(1,
                                    tags$div(
                                      HTML("<strong>Export:</strong>"),
                                      actionButton("export_allde", "Export", width = "100%"))
                             )
                           ),
                           
                           fluidRow(
                             column(6, wellPanel(plotOutput("heatmapde", height = "800px"))),
                             column(3, checkboxInput("include", "Include All Subsets in heatmap",value=T, width='100%')),
                             column(3, selectInput("scale", "Include Scaled data in heatmap", choices=c("data", 'scale.data'),selected='scale.data',width='100%'))
                             
                           )
                    )
                  )
                )

# Label Prediction --------------------------------------------------------

    )
,tabPanel("Label Prediction",
               fluidPage(
                 fluidRow(
                   column(3,
                          #wellPanel(
                            #titlePanel("Query Dataset"),
                            #uiOutput("dataset_selector")
                          #),
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
                            column(4, wellPanel(plotOutput("umap_plotq", height = "400px"))),
                            column(4, wellPanel(plotOutput("umap_plotr", height = "400px"))),
                            column(4, wellPanel(plotOutput("umap_plotqr", height = "400px")))
                          ),
                          fluidRow(
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
                            column(6, wellPanel(plotOutput("heatmaptrans", height = "400px"))),
                            column(6,
                                   tags$div(
                                     HTML("<strong>Export data & plots:</strong>"),
                                     actionButton("export_alltrans", "Export UMAP & DotPlot", width = "100%"))
                            )
                          )
                   )
                 )
               )

    )
  )
);server <- function(input, output, session) {
  
  # Reactive variable for Single Cell dataset
  seurat_data <- reactiveVal(NULL)
  hm <- reactiveVal(NULL)
  vol <- reactiveVal(NULL)
  mark <- reactiveVal(NULL)
  
  # Dataset selector for Single Cell
  output$dataset_selector <- renderUI({
    rds_files <- list.files(path = "sc/", pattern = "*.rds", full.names = FALSE)
    selectInput("selected_dataset", "Choose a Single Cell dataset:", choices = rds_files, selected = NULL)
  })
  
  # Load dataset and update UI
  observeEvent(input$update_dataset, {
    req(input$selected_dataset)
    data <- readRDS(file.path("sc/", input$selected_dataset))
    seurat_data(data)

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
      UMAPPlot(seurat_data(), label = TRUE, label.box = TRUE, pt.size = input$ptsize) + NoLegend() + NoAxes()
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
      FeaturePlot(seurat_data(), features = feature_gene, cols = c('grey', "red"), pt.size = input$ptsize, order = TRUE, reduction='umap') + 
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
    # After updating seurat_data(), update the ident1 and ident2 dropdowns
    updateSelectInput(session, "ident1", choices = unique(Idents(seurat_data())), selected = NULL)
    updateSelectInput(session, "ident2", choices = unique(Idents(seurat_data())), selected = NULL)
    
    # Update Idents
    modified_idents <- sapply(1:length(unique(Idents(seurat_data()))), function(i) {
      input[[paste0("ident_", i)]]
    })
    
    names(modified_idents) <- levels(Idents(seurat_data()))
    seurat_data(RenameIdents(seurat_data(), modified_idents))
    # After updating seurat_data(), update the ident1 and ident2 dropdowns
    updateSelectInput(session, "ident1", choices = unique(Idents(seurat_data())), selected = NULL)
    updateSelectInput(session, "ident2", choices = unique(Idents(seurat_data())), selected = NULL)
    
    # UMAP Plot
    output$umap_plot <- renderPlot({
      UMAPPlot(seurat_data(), label = TRUE, label.box = TRUE, pt.size = input$ptsize) + NoLegend() + NoAxes()
    })
    output$umap_plot2 <- renderPlot({
      UMAPPlot(seurat_data(), label = TRUE, label.box = TRUE, pt.size = input$ptsize) + NoLegend() + NoAxes()
    })
    output$umap_plotq <- renderPlot({
      UMAPPlot(seurat_data(), label = TRUE, label.box = TRUE, pt.size = input$ptsize) + NoLegend() + NoAxes()
    })
    
    # DotPlot
    output$dot_plot <- renderPlot({
      genes <- str_split(input$dotplot_genes, ",")[[1]]
      genes <- str_trim(genes)
      genes <- genes[genes %in% rownames(seurat_data())]
      DotPlot(seurat_data(), features = genes, cols = c('grey', 'red')) + NoLegend() + RotatedAxis()
    })
  })
  # Ident1 dropdown (to select a specific level from the groupby column)
  output$ident1 <- renderUI({
    req(seurat_data())
    selectInput("ident1", "Select Ident 1:", choices = NULL)
  })
  
  # Ident2 dropdown (to select a specific level from the groupby column)
  output$ident2 <- renderUI({
    req(seurat_data())
    selectInput("ident2", "Select Ident 2:", choices = NULL)
  })
  
  observeEvent(input$DE, {
    req(seurat_data(), input$ident1, input$ident2, input$include, input$signif, input$min_diff_pct, input$min_pct, input$logfc, input$scale, input$num)
    
    # Extract unique levels from the metadata
    data <- seurat_data()
    
    marks <- FindMarkers(
      data, 
      ident.1 = input$ident1, 
      ident.2 = input$ident2,
      min.pct = input$min_pct, 
      min.diff.pct = input$min_diff_pct, 
      logfc.threshold = input$logfc
    )
    if(input$num > nrow(marks)){n = nrow(marks)} else {n = input$num}
    # Add custom columns for visualization purposes
    marks <- marks %>%
      mutate(col = ifelse(p_val < input$signif, "signif", "ns")) %>%
      mutate(label = ifelse(p_val <= input$signif & abs(avg_log2FC) > 1, "label", "na")) %>%
      mutate(gene = rownames(marks)) %>%
      filter(!grepl("MT-", gene) & !grepl("^RPS", gene) & !grepl("^RPL", gene))
    
    # Update reactive value for the markers
    mark(marks)
    
    # Heatmap generation
    output$heatmapde <- renderPlot({
      if(input$include){
        hm(DoHeatmap(
          AggregateExpression(data, return.seurat = T),
          features=rownames(rbind(marks %>%top_n(n=n, wt=avg_log2FC), 
                                  marks %>%top_n(n=-n, wt=avg_log2FC))),
          draw.lines = F, slot=input$scale
        )+NoLegend())  
        hm()
      } else {
        hm(DoHeatmap(
          AggregateExpression(subset(data, idents=c(input$ident1, input$ident2)), return.seurat = T),
            features=rownames(rbind(marks %>%top_n(n=n, wt=avg_log2FC),  marks %>%top_n(n=-n, wt=avg_log2FC))), 
            draw.lines = F) +NoLegend())
        hm()
      }
    })
    
    # Volcano plot generation
    output$volcano <- renderPlot({
      volcano_plot <- ggplot(marks, aes(x = avg_log2FC, y = -log10(p_val), color = col)) +
        geom_point() +
        scale_color_manual(values = c("blue", "red")) +
        theme_bw() +
        ggtitle(paste(input$ident2, "vs", input$ident1))
      
      vol(LabelPoints(volcano_plot, rownames(marks)[marks$label == 'label'], repel = TRUE, max.overlaps = input$maxoverlap))
      vol()
    })
    

    
  })

# Export Vis --------------------------------------------------------------

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
    print(UMAPPlot(seurat_data(), label = TRUE, label.box = TRUE, pt.size = input$ptsize) + NoLegend() + NoAxes())
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
      FeaturePlot(seurat_data(), features = input$feature_gene, cols = c('grey', "red"), order = TRUE, pt.size = input$ptsize, reduction='umap') + 
        ggtitle(ifelse(input$feature_gene %in% rownames(seurat_data()), input$feature_gene, paste0(input$feature_gene, " is not found in the dataset, displaying default CD3E"))) + 
        NoLegend() + NoAxes()
    )
    dev.off()
  })
  ##### export DE ####
  # Export UMAP, Volcano, and Heatmap with the "Export All" button
  observeEvent(input$export_allde, {
    wd <- getwd()
    req(seurat_data(), vol(), hm(), mark())
    
    prefix <- input$filename_prefix
    current_time <- format(Sys.time(), "%Y%m%d_%H%M")
    dataset_name <- gsub("\\.rds$", "", input$selected_dataset)
    
    output_dir <- file.path(wd, "output")
    if (!dir.exists(output_dir)) {
      dir.create(output_dir)
    }
    
    # File paths
    umap_filename <- file.path(output_dir, paste0("/",current_time, "_", dataset_name, "_", prefix, "_UMAP.pdf"))
    vol_filename <- file.path(output_dir, paste0(current_time, "_", dataset_name, "_", input$ident2, "_vs_", input$ident1, "_vol.pdf"))
    mark_filename <- file.path(output_dir, paste0(current_time, "_", dataset_name, "_", input$ident2, "_vs_", input$ident1, "_mark.csv"))
    hm_filename <- file.path(output_dir, paste0(current_time, "_", dataset_name, "_", input$ident2, "_vs_", input$ident1, "_hm.png"))
    
    # Export UMAP plot
    pdf(umap_filename, width = 8, height = 8)
    print(DimPlot(seurat_data(), reduction = "umap", pt.size = input$ptsize, label = TRUE) + NoAxes() + NoLegend())
    dev.off()
    
    # Export Volcano plot
    pdf(vol_filename, width = 8, height = 8)
    print(vol())
    dev.off()
    
    # Export Heatmap
    png(hm_filename, width = 800, height = 800)
    print(hm())
    dev.off()
    
    # Export Markers (CSV)
    write.csv(mark(), mark_filename)
  })
};runApp(shinyApp(ui,server))



  
