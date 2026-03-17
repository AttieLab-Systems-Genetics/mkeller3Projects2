launch_qtl_app <- function() {
  
  # ==========================================
  # 1. INSTALL & LOAD DEPENDENCIES
  # ==========================================
  required_packages <- c("shiny", "data.table", "ggplot2", "plotly", "GenomicRanges")
  bioc_packages <- c("GenomicScores", "phastCons35way.UCSC.mm39")
  
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  for (pkg in required_packages) { if(!(pkg %in% installed.packages()[,"Package"])) install.packages(pkg) }
  for (pkg in bioc_packages) { if(!(pkg %in% installed.packages()[,"Package"])) BiocManager::install(pkg) }
  
  library(shiny); library(data.table); library(ggplot2); library(plotly)
  library(GenomicScores); library(phastCons35way.UCSC.mm39); library(GenomicRanges)
  
  # ==========================================
  # 2. DATA PATHS
  # ==========================================
  base_dir <- "W:/General/main_directory/top_snps_high_moderate_impact/one_row_per_csq"
  gene_ref_path <- "W:/General/Projects2/Jeeyeon Cha/NEW integration MAFA SNPs and DEGs/mouse_genes_mm39_fixed.csv"
  
  # ==========================================
  # 3. LOAD GENOMIC REFERENCE (Fixed Logic)
  # ==========================================
  message("Loading genomic reference and phastCons scores...")
  if(!file.exists(gene_ref_path)) stop(paste("File not found:", gene_ref_path))
  
  gene_ref <- fread(gene_ref_path)
  
  # Standardize column names (handling the "fixed" CSV format)
  old_cols <- c("Gene stable ID", "Gene Symbol", "Gene start (bp)", "Gene end (bp)", "Gene type", "Chromosome")
  new_cols <- c("gene_id", "symbol", "start", "end", "type", "chr")
  
  # Check which columns actually exist to avoid setnames errors
  existing <- intersect(old_cols, colnames(gene_ref))
  if(length(existing) > 0) {
    setnames(gene_ref, old_cols, new_cols, skip_absent = TRUE)
  }
  
  # CLEAN CHROMOSOME COLUMN SAFELY
  # We use [[ ]] to ensure we are referencing the column, not an external object
  if("chr" %in% colnames(gene_ref)) {
    gene_ref[, chr := gsub("chr", "", as.character(get("chr")), ignore.case = TRUE)]
  } else {
    stop("Mapping Error: Could not find 'Chromosome' column. Check your CSV headers.")
  }
  
  mm39_lens <- data.table(
    chr_factor = factor(c(1:19, "X"), levels = c(1:19, "X")),
    length_mb = c(195.1, 181.7, 159.5, 155.9, 152.0, 149.3, 144.9, 130.1, 124.3, 130.5, 
                  121.9, 120.0, 120.9, 125.1, 104.0, 98.0, 95.3, 90.7, 61.3, 169.3)
  )
  mm39_lens[, offset := c(0, cumsum(length_mb)[-.N])]
  phast_obj <- getGScores("phastCons35way.UCSC.mm39")
  
  # ==========================================
  # 4. SCAN FOR QTL FILES
  # ==========================================
  message("Scanning for QTL data files...")
  all_files <- list.files(base_dir, pattern = "\\.csv$")
  trait_files <- all_files[!grepl("gene|isoform", all_files, ignore.case = TRUE)]
  
  get_trait_display_name <- function(fname) {
    name <- gsub("DO1200_", "", fname)
    name <- gsub("_all_mice_(additive|diet_interactive|sex_interactive).*", "", name)
    tools::toTitleCase(gsub("_", " ", name))
  }
  
  trait_map <- list()
  unique_traits <- unique(sapply(trait_files, get_trait_display_name))
  for (tr in unique_traits) {
    file_pattern <- gsub(" ", "_", tolower(tr))
    f_add  <- trait_files[grepl(paste0(file_pattern, ".*_additive"), trait_files)]
    f_diet <- trait_files[grepl(paste0(file_pattern, ".*_diet_interactive"), trait_files)]
    f_sex  <- trait_files[grepl(paste0(file_pattern, ".*_sex_interactive"), trait_files)]
    
    if (length(f_add) > 0 || length(f_diet) > 0 || length(f_sex) > 0) {
      models <- list()
      if(length(f_add) > 0)  models$additive <- f_add[1]
      if(length(f_diet) > 0) models$diet_interactive <- f_diet[1]
      if(length(f_sex) > 0)  models$sex_interactive <- f_sex[1]
      trait_map[[tr]] <- models
    }
  }
  
  # ==========================================
  # 5. UI & SERVER
  # ==========================================
  ui <- fluidPage(
    titlePanel("Interactive QTL Fine-Mapper v4.4"),
    sidebarLayout(
      sidebarPanel(
        width = 2,
        selectInput("trait_class", "Trait Class:", choices = names(trait_map)),
        uiOutput("model_selector"),
        hr(),
        sliderInput("lod_thresh", "Min Variant LOD:", min = 0, max = 25, value = 3, step = 0.5),
        sliderInput("prio_thresh", "Min Priority:", min = -1.1, max = 2.1, value = -1.1, step = 0.1),
        hr(),
        sliderInput("schematic_window", "Window (kb):", min = 20, max = 1000, value = 100, step = 20),
        actionButton("reset_plot", "Reset Selection", class = "btn-block btn-info"),
        uiOutput("metadataSidebar")
      ),
      mainPanel(
        width = 10,
        plotlyOutput("manhattanPlot", height = "520px"),
        hr(),
        plotOutput("geneSchematic", height = "480px")
      )
    )
  )
  
  server <- function(input, output, session) {
    curr_sel <- reactiveVal(NULL)
    observeEvent(input$reset_plot, { curr_sel(NULL) })
    
    output$model_selector <- renderUI({
      req(input$trait_class)
      avail_models <- names(trait_map[[input$trait_class]])
      choice_vec <- c("Additive" = "additive", "Diet Interactive" = "diet_interactive", "Sex Interactive" = "sex_interactive")
      radioButtons("model_type", "Scan Type:", choices = choice_vec[choice_vec %in% avail_models])
    })
    
    qtl_data <- reactive({
      req(input$trait_class, input$model_type)
      fname <- trait_map[[input$trait_class]][[input$model_type]]
      dt <- fread(file.path(base_dir, fname))[variant_type == "snp"]
      dt[, pos_bp := as.numeric(variant_pos) * 1e6]
      dt[, chr_factor := factor(gsub("chr", "", as.character(variant_chr)), levels = levels(mm39_lens$chr_factor))]
      gr <- GRanges(seqnames = paste0("chr", dt$chr_factor), ranges = IRanges(start = dt$pos_bp, width = 1))
      dt$phast_score <- gscores(phast_obj, gr)$default; dt[is.na(phast_score), phast_score := 0]
      sift_numeric <- suppressWarnings(as.numeric(dt$sift_score))
      dt[, priority_score := phast_score - ifelse(is.na(sift_numeric), 1.0, sift_numeric)]
      dt[impact == "HIGH", priority_score := phast_score + 1] 
      plot_dt <- dt[, .SD[which.max(priority_score)], by = .(variant_id)]
      plot_dt <- merge(plot_dt, mm39_lens[, .(chr_factor, offset)], by = "chr_factor")
      plot_dt[, global_pos_mb := variant_pos + offset]
      return(plot_dt)
    })
    
    output$manhattanPlot <- renderPlotly({
      dt <- qtl_data()[variant_lod >= input$lod_thresh & priority_score >= input$prio_thresh]
      req(nrow(dt) > 0)
      p <- ggplot(dt, aes(x = global_pos_mb, y = priority_score, key = variant_id)) +
        geom_rect(data = mm39_lens[as.numeric(chr_factor) %% 2 == 1], 
                  aes(xmin = offset, xmax = offset + length_mb, ymin = min(dt$priority_score)-0.2, ymax = max(dt$priority_score)+0.2), 
                  fill = "grey96", inherit.aes = FALSE) +
        geom_point(aes(fill = phast_score, size = variant_lod, alpha = phast_score,
                       text = paste0("Gene: ", gene_symbol, "<br>LOD: ", round(variant_lod, 2))), 
                   shape = 21, color = "black", stroke = 0.15) +
        scale_fill_gradientn(colors = c("#0571b0", "#f4a582", "#ca0020"), limits = c(0, 1)) +
        scale_size_continuous(range = c(2.5, 8.5), limits = c(0, 25)) + 
        scale_alpha_continuous(range = c(0.15, 1.0), limits = c(0, 1)) +
        scale_x_continuous(breaks = mm39_lens$offset + (mm39_lens$length_mb / 2), labels = mm39_lens$chr_factor) +
        theme_minimal() + labs(x = "Chromosome", y = "Priority Score") +
        theme(legend.position = "none", axis.text.x = element_text(size = 14, face = "bold"))
      
      if (!is.null(curr_sel())) {
        sel_dt <- dt[variant_id == curr_sel()]
        if (nrow(sel_dt) > 0) {
          p <- p + geom_point(data = sel_dt, aes(x = global_pos_mb, y = priority_score, size = variant_lod),
                              fill = "orange", color = "black", stroke = 0.8, alpha = 1, inherit.aes = FALSE, shape = 21)
        }
      }
      ggplotly(p, tooltip = "text", source = "manhattan") %>% layout(clickmode = "event") %>% style(hoverinfo = "text", traces = 1)
    })
    
    observe({
      ev <- event_data("plotly_click", source = "manhattan")
      if(!is.null(ev)) curr_sel(ev$key)
    })
    
    selected_snp <- reactive({ req(curr_sel()); qtl_data()[variant_id == curr_sel()][1] })
    
    output$metadataSidebar <- renderUI({
      req(selected_snp()); snp <- selected_snp()
      div(class="meta-card", div(class="gene-title", snp$gene_symbol), div(class="meta-section", "Variant Info"),
          div(class="meta-row", span(class="meta-label", "ID: "), span(class="meta-val", snp$rs_number)),
          div(class="meta-row", span(class="meta-label", "Pos: "), span(class="meta-val", sprintf("%.1f Mb", snp$variant_pos))),
          div(class="meta-row", span(class="meta-label", "Csq: "), span(class="meta-val", snp$csq)),
          div(class="meta-row", span(class="meta-label", "Var LOD: "), span(class="meta-val", round(snp$variant_lod, 1))),
          br(), div(class="meta-section", "Phenotype Info"),
          div(class="meta-row", span(class="meta-label", "Trait: "), span(class="meta-val", snp$phenotype)))
    })
    
    output$geneSchematic <- renderPlot({
      snp <- selected_snp(); req(snp); win_mb <- input$schematic_window / 1000
      local_genes <- gene_ref[chr == gsub("chr", "", snp$variant_chr) & start < (snp$variant_pos + win_mb)*1e6 & end > (snp$variant_pos - win_mb)*1e6]
      local_genes <- local_genes[order(start)]; local_genes[, y_level := 1]
      if(nrow(local_genes) > 1){
        for(i in 2:nrow(local_genes)){
          overlaps <- local_genes[1:(i-1)][end > (local_genes$start[i] - (win_mb * 0.15 * 1e6))]
          if(nrow(overlaps) > 0) local_genes$y_level[i] <- (max(overlaps$y_level) %% 5) + 1
        }
      }
      ggplot(local_genes) +
        annotate("text", x = snp$variant_pos, y = 7.5, label = paste(snp$phenotype, "QTL"), color = "darkred", fontface = "bold", size = 8.5) +
        geom_segment(aes(x = snp$variant_pos, xend = snp$variant_pos, y = 0.5, yend = 7.1), color = "darkred", linetype = "dashed", linewidth = 0.9) +
        geom_rect(aes(xmin = start/1e6, xmax = end/1e6, ymin = y_level*1.1, ymax = y_level*1.1 + 0.45, fill = (symbol == snp$gene_symbol)), color = "black") +
        geom_text(aes(x = (start+end)/2e6, y = y_level*1.1 + 0.75, label = symbol, size = (symbol == snp$gene_symbol), fontface = "bold")) +
        scale_size_manual(values = c("FALSE" = 5.5, "TRUE" = 9)) + scale_fill_manual(values = c("TRUE" = "orange", "FALSE" = "grey92")) +
        scale_x_continuous(labels = scales::comma, name = "Position (Mbp)") + scale_y_continuous(limits = c(0.5, 8.5)) +
        theme_minimal() + theme(axis.text.y = element_blank(), axis.title.y = element_blank(), panel.grid = element_blank(), legend.position = "none")
    })
  }
  
  shinyApp(ui, server)
}