wd <- getwd()

ui <- fluidPage(
  fluidRow(
    column(3,
           wellPanel(
             titlePanel("Load Dataset"),
             uiOutput("dataset_selector"),
             actionButton("update_dataset", "Load Dataset", width = "100%")
           ),
           wellPanel(
             titlePanel("Group by:"),
             uiOutput("groupby"),
             actionButton("update_idents", "Select Group", width = "100%")
             
           ),           
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
             column(6, wellPanel(plotOutput("umap_plot", height = "400px"))),
             column(6, wellPanel(plotOutput("volcano", height = "400px")))
           ),
           fluidRow(             
             column(1, textInput("filename_prefix", "Filename:", value = "user", width = "100%")),
             column(1, numericInput("umap_ptsize", "Pt Size", value = 1.2, min = 0.1, width='100%')),
             column(1, numericInput("num", "# Genes", value = 30, min = 1, width='100%')),
             column(2, numericInput("signif", "pval filter", value = 0.001, width='100%')),
             
             column(1, 
                    tags$div(
                      HTML("<strong>Run DE:</strong>"),
                      actionButton("DE", "DE", width = "100%"))
             ),
             column(1, numericInput("logfc", "logFC", value = 0.1, min = 0.01, width='100%')),
             column(2, numericInput("min_pct", "Minimum % expression", value = 0.1, min = 0.01, width='100%')),
             column(2, numericInput("min_diff_pct", "Minimum diff %", value = 0.1, min = 0.01, width='100%')),
             
             column(1,
                    tags$div(
                      HTML("<strong>Export:</strong>"),
                      actionButton("export_all", "Export", width = "100%"))
             )
           ),

           fluidRow(
             column(6, wellPanel(plotOutput("heatmap", height = "800px"))),
             column(3, checkboxInput("include", "Include All Subsets in heatmap",value=T, width='100%')),
             column(3, selectInput("scale", "Include Scaled data in heatmap", choices=c("data", 'scale.data'),selected='scale.data',width='100%'))
             
           )
    )
  )
)

server <- function(input, output, session) {
  
  # Reactive value for the Seurat dataset
  seurat_data <- reactiveVal(NULL)
  hm <- reactiveVal(NULL)
  vol <- reactiveVal(NULL)
  mark <- reactiveVal(NULL)
  
  # Dataset selector for loading the dataset
  output$dataset_selector <- renderUI({
    query_files <- list.files(path = "sc/", pattern = "*.rds", full.names = FALSE)
    selectInput("selected_dataset", "Choose a dataset:", choices = query_files, selected = NULL)
  })

# Read in Data ------------------------------------------------------------
  # Read in the dataset and populate the groupby dropdown
  observeEvent(input$update_dataset, {
    req(input$selected_dataset)
    print(paste0("Loading data:",input$selected_dataset))
    # Load the dataset
    data <- readRDS(file.path("sc/", input$selected_dataset))
    data$default <- Idents(data)
    seurat_data(data)
    updateSelectInput(session, "groupby", choices = colnames(data@meta.data), selected = 'default')
  })
  
  # Update Idents+UMAP ------------------------------------------------------
  observeEvent(input$update_idents, {
    req(seurat_data(), input$groupby)

    # Extract the unique levels from the selected metadata column
    selected_group <- seurat_data()@meta.data[[input$groupby]]
    idents <- unique(selected_group)
    
    # Populate the ident1 and ident2 dropdowns
    updateSelectInput(session, "ident1", choices = idents, selected = NULL)
    updateSelectInput(session, "ident2", choices = idents, selected = NULL)
    
    output$umap_plot <- renderPlot({
      req(seurat_data(), input$groupby)
        DimPlot(seurat_data(), reduction = "umap", group.by = input$groupby, pt.size = input$umap_ptsize, label=T)+NoAxes()+ggtitle("Dataset")+NoLegend()
    })
    
  })
  
  # Group by dropdown (to select metadata columns)
  output$groupby <- renderUI({
    req(seurat_data())
    selectInput("groupby", "Select Group By:", choices = colnames(seurat_data()@meta.data), selected = 'default')
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
  

# Output Diff analysis ----------------------------------------------------

  observeEvent(input$DE, {
    req(seurat_data(), input$groupby, input$ident1, input$ident2)
    
    # Extract the unique levels from the selected metadata column
    data <- seurat_data() 
    Idents(data) <- input$groupby

    marks <- FindMarkers(data,ident.1= input$ident1, ident.2=input$ident2,
                         min.pct = input$min_pct, min.pct.diff=input$min_diff_pct, logfc.threshold = input$logfc)
    marks <-marks %>%
      mutate(col =  ifelse(marks$p_val < input$signif, "signif", "ns")) %>% 
      mutate(label = ifelse(-log10(marks$p_val) > 1.3 & abs(marks$avg_log2FC) > 1, "label", "na")) %>% #here you can set your threshold. 
      mutate(gene = rownames(marks)) %>%
      filter(!grepl("MT-", gene) & !grepl("^RPS", gene) & !grepl("^RPL", gene))
    mark(marks)
    output$volcano <- renderPlot({
    vol(LabelPoints(ggplot(marks,aes(x=avg_log2FC, y = -log2(p_val), color = col))+geom_point()+scale_color_manual(values=c("blue", "red"))+theme_bw(), 
                rownames(marks)[marks$label == 'label'], repel=T, max.overlaps=15)+ggtitle(paste0( input$ident2, " vs ", input$ident1)))
    vol()
    })
    
    if(input$num > nrow(marks)){
      n = nrow(marks)
    } else {
      n = input$num
    }
      
    output$heatmap <- renderPlot({
      if(input$include){
        hm(DoHeatmap(
          AggregateExpression(data, return.seurat = T),
          features=rownames(rbind(marks %>%
                                    top_n(n=n, wt=avg_log2FC), 
                                  marks %>%
                                    top_n(n=-n, wt=avg_log2FC))
          ),draw.lines = F, slot=input$scale
        ) +NoLegend())  
        hm()
      } else {
       hm(DoHeatmap(
          AggregateExpression(subset(data, idents=c(input$ident1, input$ident2)), return.seurat = T),
          features=rownames(rbind(marks %>%
                                    top_n(n=n, wt=avg_log2FC), 
                                  marks %>%
                                    top_n(n=-n, wt=avg_log2FC))
          ),
          draw.lines = F, slot=input$scale
        ) +NoLegend())
        hm()
      }

    })
    
  })
  observeEvent(input$export_all, {
    req(seurat_data(),vol(),hm(),mark())
    prefix <- input$filename_prefix
    current_time <- format(Sys.time(), "%Y%m%d_%H%M")
    dataset_name <- gsub("\\.rds$", "", input$selected_query_dataset)
    
    # Create output directory if missing
    if (!dir.exists(file.path(wd, "output"))) {
      dir.create(file.path(wd, "output"))
    }
    umap_filename <- file.path(wd, "output", paste0(current_time, "_",dataset_name,"_" , prefix, "_UMAP.pdf"))
    vol_filename <- file.path(wd, "output", paste0(current_time, "_", dataset_name,"_" ,input$ident2, "_vs_",input$ident1,"_" , prefix, "_vol.pdf"))
    mark_filename <- file.path(wd, "output", paste0(current_time, "_",dataset_name,"_" , input$ident2, "_vs_",input$ident1,"_" , prefix, "_mark.csv"))
    hm_filename <- file.path(wd, "output", paste0(current_time, "_", dataset_name,"_" , input$ident2, "_vs_",input$ident1,"_" , prefix, "_hm.png"))
    
    pdf(umap_filename, width = width, height = height)
      print(DimPlot(seurat_data(), reduction = "umap", group.by = input$groupby, pt.size = input$umap_ptsize, label=T)+NoAxes()+ggtitle("Dataset")+NoLegend())
    dev.off()
    
    pdf(vol_filename, width = widthvol, height = heightvol)
      print(vol())
    dev.off()
    
    png(hm_filename, width=widthheat, height=heightheat)
      print(hm())
    dev.off()
    
    write.csv(mark(), mark_filename)
 })
}
runApp(shinyApp(ui,server))
