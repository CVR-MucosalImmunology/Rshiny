
ui <- fluidPage(
  column(2,
         wellPanel(
           titlePanel("Select Sample"),
           uiOutput("sample_selector"),
           actionButton("load_sample", "Update sample", width = "100%")
         ),
         wellPanel(
             titlePanel("Markers"),
             uiOutput("marker_selector"),
             actionButton("run_sample", "Run RPCA", width = "100%")
           )
  ),
  
  column(5,
         fluidRow(           
           wellPanel(titlePanel("Check integration (Pre)"),
                     plotOutput("umap_plot", height = "600px"))
         ),
         fluidRow(
           column(4,uiOutput("featselect")),
           column(4,actionButton("update_feat", "Update:", width='100%')),
           column(4,actionButton("printall", "PrintAll", width='100%'))
           
         ),
         fluidRow(
           wellPanel(plotOutput("feat_plot", height = "600px"))
         )
  ),
  
  column(5,
         fluidRow(           
           wellPanel(titlePanel("Check integration (Post)"),
                     plotOutput("umap_plot2", height = "600px"))
         ),
         fluidRow(
           column(3,numericInput("k", "Clusters (k):", value=10, width='100%')),
           column(3, textInput("prefix", "Save_prefix", value="user", width='100%')),
           column(6,actionButton("updateSoms", "Update K", value=F, width='100%'))
         ),
         fluidRow(
           column(12,wellPanel(plotOutput("som_plot", height = "600px"))
                  )         )
  )
);server <- function(input, output, session) {
  
  # Render checkboxes for available CSV files in the 'flow' directory
  output$sample_selector <- renderUI({
    query_files <- list.files(path = "flow/", pattern = "*.csv", full.names = FALSE)
    checkboxGroupInput("sample_data", "Choose flow sample data:", choices = query_files)
  })
  
  # Reactive to store the loaded data from selected files
  selected_data <- reactiveVal(list())
  selected_meta <- reactiveVal(list())
  working_data <- reactiveVal(list())
  image_data <- reactiveVal(NULL)
  clust_data <- reactiveVal(NULL)
  
  # When the load_sample button is clicked, load the selected files into a list
  observeEvent(input$load_sample, {
    if (!is.null(input$sample_data)) {
      files <- input$sample_data
      
      # Load data from the selected CSV files in 'flow/' directory
      data_list <- lapply(files, function(file) {
        filepath <- file.path("flow", file)
        read.csv(filepath)  # Load CSV file
      })
      
      # Load metadata files corresponding to the selected CSV files in 'meta/' directory
      meta_list <- lapply(files, function(file) {
        meta_filepath <- file.path("meta", paste0(substr(file, 1, nchar(file) - 4), "_meta.csv"))
        read.csv(meta_filepath)  # Load corresponding metadata CSV
      })
      
      # Store the loaded data and metadata in reactive values
      selected_data(data_list)
      selected_meta(meta_list)
      
      # Use the first column of the first data file to populate marker selection
      output$marker_selector <- renderUI({
        if (length(selected_meta()) > 0) {
          first_meta <- selected_meta()[[1]]  # Access the first metadata file
          markers <- first_meta[[1]]  # Assuming the first column contains marker names
          checkboxGroupInput("marker_data", "Choose markers to integrate:", choices = markers)
        } else {
          p("No data loaded yet.")
        }
      })
    }
  })
  output$featselect <- renderUI({
    req(selected_meta())
    selectInput("mark_selector2", "Select a marker:", choices = selected_meta()[[1]][[1]], selected = NULL)
  })
  # When the run_sample button is clicked, process the selected data
  observeEvent(input$run_sample, {
    req(input$marker_data)  # Ensure markers are selected
    
    markers <- input$marker_data  # Get selected markers
    data_list <- selected_data()  # Get selected data
    meta_list <- selected_meta()  # Get selected metadata
    file_names <- input$sample_data    
    
    processed_data_list <- lapply(1:length(data_list), function(i) {
      d <- data_list[[i]]    # Current data file
      d2 <- d[, 1, drop = FALSE]  # Keep the first column for merging
      d3 <- d[, 1, drop = FALSE]  # Initialize filter column for filtering later
      m <- meta_list[[i]]    # Current metadata file
      
      # Filter metadata by the selected markers
      m <- m[m$X %in% markers,]
      
      # Loop through each selected marker
      for (j in 1:nrow(m)) {
        marker_name <- m$X[j]
        val <- d[, colnames(d) == marker_name]
        
        # Apply the quantile filtering and transformation
        min_val <- quantile(val, probs = (m$filt[j] / 100))
        max_val <- quantile(val, probs = 0.999)
        filt <- ifelse(val > max_val | val < min_val, 1, 0)
        
        # Arcsinh transformation with cofactor
        val <- val / m$cofactor[j]
        val <- asinh(val)
        
        # Add the processed column to d2 and the filter column to d3
        val_df <- data.frame(x = val)
        colnames(val_df) <- marker_name
        filt_df <- data.frame(x = filt)
        colnames(filt_df) <- marker_name
        
        d2 <- cbind(d2, val_df)
        d3 <- cbind(d3, filt_df)
      }
      
      # Remove the first column from d2 and d3
      d2 <- d2[, -1]
      d3 <- d3[, -1]
      
      # Filter out rows where any marker fails the filter
      d4 <- d2[rowSums(d3) == 0, ]
      d4$sample <- gsub("\\.csv$", "", file_names[i])
      return(d4)
    })
    processed_data_list2 <- lapply(processed_data_list, function(df) {
      if (nrow(df) > 6000) {
        df[sample(1:nrow(df), 6000), ]
      } else {
        df  # Return the original df if it has fewer 
      }
    })
    # Store the processed data in working_data
    working_data(processed_data_list)
    # You can now access processed data in working_data for further processing
    combined_data <- do.call(rbind, processed_data_list2)
    d5 <- run.umap(as.data.table(combined_data),markers, umap.x.name = 'UMAP_X', umap.y.name = 'UMAP_Y')
    d6 <- run.rpca(dat = d5, use.cols = markers, batch.col = 'sample')
    d6 <- run.umap(d6,paste0(markers, "_rPCA_aligned"), umap.x.name = 'UMAP_X_rPCA', umap.y.name = 'UMAP_Y_rPCA')
    
    image_data(d6)
    
    output$umap_plot <- renderPlot({
      ggplot(d6, aes(UMAP_X,UMAP_Y, color=sample))+geom_point()+theme_bw();
    })
    
    output$umap_plot2 <- renderPlot({
      ggplot(d6, aes(UMAP_X_rPCA,UMAP_Y_rPCA, color=sample))+geom_point()+theme_bw();
    })

  })
  
  observeEvent(input$update_feat, {
    req(image_data(), input$mark_selector2)
    
    output$feat_plot <- renderPlot({
      ggplot(image_data(),
             aes_string('UMAP_X_rPCA','UMAP_Y_rPCA', 
                        
                        color=paste0(input$mark_selector2, "_rPCA_aligned")
                        #color=input$mark_selector2
                        
                        )) +
        geom_point() +
        scale_color_gradientn(colors = c("white","grey","blue","darkblue")) +
        theme_bw()
      })
  })
  
  observeEvent(input$updateSoms, {
    req(input$k)
    dt<- image_data()
    markers <- input$marker_data  # Get selected markers
    dt <- run.flowsom(dt, paste0(input$marker_data, "_rPCA_aligned"), meta.k = input$k)
    output$som_plot <- renderPlot({
      make.colour.plot(dt, "UMAP_X_rPCA", "UMAP_Y_rPCA", "FlowSOM_metacluster", col.type = 'factor', add.label = TRUE, blank.axis=T,save.to.disk = F)
    })

    clust_data(dt)
  })
  
  observeEvent(input$printall, {
    req(input$marker_data, clust_data())
    current_time <- format(Sys.time(), "%Y%m%d_%H%M")
    
    clustered_data_name = paste0("output/rds/",current_time,"_",input$prefix,"_clustered_data.rds")
    feature_name = paste0("output/features/",current_time,"_",input$prefix,"_featurePlot_")
    clustered_soms_name = paste0("output/",current_time,"_",input$prefix,"_soms_clusters.pdf")
    integrated_name = paste0("output/",current_time,"_",input$prefix,"_integrated_umap.pdf")
    info_name = paste0("output/rds/",current_time,"_",input$prefix,"_integration_info.rds")
    
    saveRDS(clust_data(), clustered_data_name)
    
    pdf(integrated_name, width=8,height=5)
      print(ggplot(clust_data(), aes(UMAP_X,UMAP_Y, color=sample))+geom_point()+theme_bw())
      print(ggplot(clust_data(), aes(UMAP_X_rPCA,UMAP_Y_rPCA, color=sample))+geom_point()+theme_bw());
    dev.off()
    dat <-   clust_data()
    
    make.colour.plot(dat, "UMAP_X_rPCA", "UMAP_Y_rPCA", "FlowSOM_metacluster", col.type = 'factor', add.label = TRUE, blank.axis=T, path='output')

    markers <- input$marker_data
    for(i in 1:length(markers)){
      pdf(paste0(feature_name, str_split(markers[i], "\\.")[[1]][1], ".pdf"), height=5, width=8)
      print(ggplot(image_data(),
                   aes_string('UMAP_X_rPCA','UMAP_Y_rPCA', 
                              color=paste0(markers[i], "_rPCA_aligned")
                   )) +
              geom_point() +
              scale_color_gradientn(colors = c("white","grey","blue","darkblue")) +
              theme_bw())
      dev.off()
    }
    info = list(samples = c(unique(image_data()$sample)), 
                markers = grep("aligned", colnames(image_data())))
    saveRDS(info, info_name)
    exp <- do.aggregate(dat, paste0(input$marker_data, "_rPCA_aligned"), by = "FlowSOM_metacluster")
    make.pheatmap(exp, "FlowSOM_metacluster", paste0(input$marker_data, "_rPCA_aligned"),path='output')
    
  })
}
