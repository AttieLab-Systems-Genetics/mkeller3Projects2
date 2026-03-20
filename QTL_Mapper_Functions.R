launch_qtl_app <- function() {
  
  # ==========================================
  # 1. DEPENDENCIES
  # ==========================================
  library(shiny); library(data.table); library(ggplot2); library(plotly)
  library(GenomicScores); library(phastCons35way.UCSC.mm39); library(GenomicRanges)
  
  base_dir <- "W:/General/main_directory/top_snps_high_moderate_impact/one_row_per_csq"
  gene_ref_path <- "W:/General/Projects2/Jeeyeon Cha/NEW integration MAFA SNPs and DEGs/mouse_genes_mm39_fixed.csv"
  
  # ==========================================
  # 2. HARMONIZER: GENE REFERENCE
  # ==========================================
  message("Harmonizing Gene Reference...")
  gene_ref <- fread(gene_ref_path)
  
  # Greedy search: Find columns by keyword to handle variations in CSV headers
  curr_cols <- colnames(gene_ref)
  find_col <- function(pattern) curr_cols[grepl(pattern, curr_cols, ignore.case = TRUE)][1]
  
  # Map standard internal names
  target_map <- c(
    gene_id = find_col("stable|ID"),
    symbol  = find_col("Symbol"),
    chr     = find_col("Chrom|Chr"),
    start   = find_col("start"),
    end     = find_col("end"),
    strand  = find_col("Strand"),
    gene_type = find_col("Gene type|Biotype")
  )
  
  # Rename columns based on matches found
  for (internal_name in names(target_map)) {
    actual_col <- target_map[internal_name]
    if (!is.na(actual_col)) {
      setnames(gene_ref, actual_col, internal_name)
    }
  }
  
  # Final safety check for the 'chr' column
  if (!"chr" %in% names(gene_ref)) {
    stop("Could not find a Chromosome column. Headers found: ", paste(curr_cols, collapse=", "))
  }
  
  # Clean 'chr' safely to avoid scope bugs
  gene_ref[, chr := gsub("chr", "", as.character(.SD[["chr"]]), ignore.case = TRUE)]
  
  # ==========================================
  # 3. HARMONIZER: GENOMIC SCALES
  # ==========================================
  mm39_lens <- data.table(
    chr_std = factor(c(1:19, "X"), levels = c(1:19, "X")),
    length_mb = c(195.1, 181.7, 159.5, 155.9, 152.0, 149.3, 144.9, 130.1, 124.3, 130.5, 
                  121.9, 120.0, 120.9, 125.1, 104.0, 98.0, 95.3, 90.7, 61.3, 169.3)
  )
  mm39_lens[, offset := c(0, cumsum(length_mb)[-.N])]
  phast_obj <- getGScores("phastCons35way.UCSC.mm39")
  
  # ==========================================
  # 4. DYNAMIC SCANNING
  # ==========================================
  message("Scanning for QTL data files...")
  all_files <- list.files(base_dir, pattern = "\\.csv$")
  trait_files <- all_files[!grepl("gene|isoform", all_files, ignore.case = TRUE)]
  
  get_trait_name <- function(f) tools::toTitleCase(gsub("_", " ", gsub("DO1200_|_all_mice.*", "", f)))
  
  trait_map <- list()
  for (f in trait_files) {
    tr <- get_trait_name(f)
    model <- ifelse(grepl("additive", f), "additive", 
                    ifelse(grepl("diet", f), "diet_interactive", "sex_interactive"))
    if (is.null(trait_map[[tr]])) trait_map[[tr]] <- list()
    trait_map[[tr]][[model]] <- f
  }
  
  # ==========================================
  # 5. UI & SERVER
  # ==========================================
  ui <- fluidPage(
    titlePanel("Harmonized QTL Fine-Mapper v5.0"),
    sidebarLayout(
      sidebarPanel(
        width = 2,
        selectInput("trait", "Phenotype:", choices = names(trait_map)),
        uiOutput("model_ui"),
        hr(),
        sliderInput("lod", "Min LOD:", 0, 25, 5, 0.5),
        sliderInput("prio", "Min Priority:", -1.1, 2.1, -1.1, 0.1),
        sliderInput("window", "Schematic Zoom (± Mbp):", 0.01, 1.0, 0.1, 0.01),
        actionButton("reset", "Reset Selection", class = "btn-info btn-block"),
        br(), uiOutput("meta")
      ),
      mainPanel(
        width = 10,
        plotlyOutput("manhattan", height = "520px"),
        hr(),
        plotOutput("schematic", height = "400px")
      )
    )
  )
  
  server <- function(input, output, session) {
    curr_snp <- reactiveVal(NULL)
    observeEvent(input$reset, { curr_snp(NULL) })
    
    output$model_ui <- renderUI({
      req(input$trait)
      radioButtons("model", "Model:", choices = names(trait_map[[input$trait]]))
    })
    
    qtl_data <- reactive({
      req(input$trait, input$model)
      dt <- fread(file.path(base_dir, trait_map[[input$trait]][[input$model]]))[variant_type == "snp"]
      
      # Harmonize QTL Data
      curr_cols <- colnames(dt)
      
      # Robustly find 'chr' column
      chr_col <- curr_cols[grepl("^(variant_)?(chr|chrom|chromosome)$", curr_cols, ignore.case = TRUE)][1]
      if (!is.na(chr_col) && chr_col != "chr") {
        setnames(dt, chr_col, "chr")
      }
      
      # Robustly find 'symbol' column
      sym_col <- curr_cols[grepl("^(gene_)?symbol$", curr_cols, ignore.case = TRUE)][1]
      if (!is.na(sym_col) && sym_col != "symbol") {
        setnames(dt, sym_col, "symbol")
      }
      
      # Robustly find 'variant_pos' column
      pos_col <- curr_cols[grepl("^(variant_)?pos(ition)?$", curr_cols, ignore.case = TRUE)][1]
      if (!is.na(pos_col) && pos_col != "variant_pos") {
        setnames(dt, pos_col, "variant_pos")
      }
      
      # Validation check
      if (!"chr" %in% names(dt)) {
        stop("Could not find a chromosome column. Available headers: ", paste(curr_cols, collapse=", "))
      }
      
      # Clean 'chr' safely to avoid scope bugs
      dt[, chr := gsub("chr", "", as.character(.SD[["chr"]]), ignore.case = TRUE)]
      dt[, pos_bp := as.numeric(variant_pos) * 1e6]
      
      # PhastCons & Priority Scoring
      gr <- GRanges(seqnames = paste0("chr", dt$chr), ranges = IRanges(start = dt$pos_bp, width = 1))
      dt$phast <- gscores(phast_obj, gr)$default
      dt[is.na(phast), phast := 0]
      
      sift <- suppressWarnings(as.numeric(dt$sift_score))
      dt[, priority := phast - ifelse(is.na(sift), 1.0, sift)]
      dt[impact == "HIGH", priority := phast + 1]
      
      # Positioning
      dt[, chr_std := factor(chr, levels = levels(mm39_lens$chr_std))]
      dt <- merge(dt, mm39_lens, by = "chr_std")
      dt[, global_pos := variant_pos + offset]
      dt
    })
    
    output$manhattan <- renderPlotly({
      dt <- qtl_data()[variant_lod >= input$lod & priority >= input$prio]; req(nrow(dt) > 0)
      p <- ggplot(dt, aes(x = global_pos, y = priority, key = variant_id)) +
        geom_point(aes(fill = phast, size = variant_lod, 
                       text = paste("Trait:", phenotype, "\nChr:", chr, "\nPos (Mbp):", round(pos_bp/1e6, 3), 
                                    "\nPriority:", round(priority, 3), "\nConsequence:", csq),
                       alpha = ifelse(phast <= 0, 0.2, 1.0)), 
                   shape = 21, color = "black", stroke = 0.2) +
        scale_fill_gradientn(colors = c("#0571b0", "#f4a582", "#ca0020")) +
        scale_alpha_identity() +
        scale_x_continuous(breaks = mm39_lens$offset + (mm39_lens$length_mb/2), labels = mm39_lens$chr_std) +
        theme_minimal() + 
        theme(axis.title = element_text(size = 14), axis.text = element_text(size = 12))
        
      if (!is.null(curr_snp())) {
        p <- p + geom_point(data = dt[variant_id == curr_snp()], fill = "orange", size = 5, shape = 21)
      }
      
      suppressWarnings(
        ggplotly(p, source = "manhattan", tooltip = "text") %>% 
          layout(clickmode = "event", uirevision = "constant",
                 xaxis = list(titlefont = list(size = 14), tickfont = list(size = 12)),
                 yaxis = list(titlefont = list(size = 14), tickfont = list(size = 12))) %>%
          event_register("plotly_click")
      )
    })
    
    observe({ ev <- event_data("plotly_click", source = "manhattan"); if(!is.null(ev)) curr_snp(ev$key) })
    
    output$meta <- renderUI({
      req(curr_snp()); s <- qtl_data()[variant_id == curr_snp()][1]
      wellPanel(
        h4(s$symbol), 
        p(HTML(paste0("<b>Trait:</b> ", s$phenotype))),
        p(HTML(paste0("<b>QTL LOD:</b> ", round(s$qtl_lod, 2)))),
        p(HTML(paste0("<b>QTL Pos:</b> ", round(s$qtl_pos, 2), " Mbp"))),
        hr(style = "margin-top: 5px; margin-bottom: 5px;"),
        p("RSID: ", s$rs_number), 
        p("Variant LOD: ", round(s$variant_lod, 2)), 
        p("Csq: ", s$csq),
        p(HTML(paste0("<b>AA Change:</b> ", s$aa_change, "<br><b>AA Pos:</b> ", s$aa_pos)))
      )
    })
    
    output$schematic <- renderPlot({
      req(curr_snp(), input$window); s <- qtl_data()[variant_id == curr_snp()][1]
      
      # Use dynamic window size from UI slider
      local_genes <- gene_ref[chr == s$chr & start < (s$variant_pos + input$window)*1e6 & end > (s$variant_pos - input$window)*1e6]
      if ("gene_type" %in% names(local_genes)) {
        local_genes <- local_genes[gene_type == "protein_coding" | symbol == s$symbol]
      }
      
      if(nrow(local_genes) > 0) {
        # Stagger genes vertically
        local_genes <- local_genes[order(start)]
        local_genes[, stagger := (seq_len(.N) %% 4) * 1.5]
        local_genes[, ymin := stagger]
        local_genes[, ymax := stagger + 1]
        
        # Determine TSS and arrow direction
        # 0.005 Mbp length for arrow
        local_genes[, tss_x := ifelse(strand == -1, end/1e6, start/1e6)]
        local_genes[, arrow_xend := ifelse(strand == -1, tss_x - 0.005, tss_x + 0.005)]
      } else {
        local_genes[, `:=`(stagger=numeric(), ymin=numeric(), ymax=numeric(), tss_x=numeric(), arrow_xend=numeric())]
      }

      p_schem <- ggplot(local_genes) +
        geom_rect(aes(xmin=start/1e6, xmax=end/1e6, ymin=ymin, ymax=ymax, fill=(symbol == s$symbol)), color="black", linewidth=0.2) +
        geom_text(aes(x=(start+end)/2e6, y=ymax+0.4, label=symbol, size=ifelse(symbol == s$symbol, 7, 5.25))) +
        scale_size_identity()
        
      if(nrow(local_genes[symbol == s$symbol]) > 0) {
         p_schem <- p_schem + geom_segment(data=local_genes[symbol == s$symbol], aes(x=tss_x, xend=arrow_xend, y=ymin+0.5, yend=ymin+0.5), arrow=arrow(length=unit(0.25, "cm"), type="closed"))
      }
      
      # Add dashed line and AA annotation for the SNP
      p_schem <- p_schem + geom_vline(xintercept = s$variant_pos, linetype = "dashed", color = "red", alpha = 0.5)
      
      # Format AA change if available (e.g., C/* to C48*)
      if (!is.na(s$aa_change) && s$aa_change != "" && !is.na(s$aa_pos)) {
        parts <- strsplit(s$aa_change, "/")[[1]]
        if (length(parts) == 2) {
          aa_str <- paste0(parts[1], s$aa_pos, parts[2])
          max_y <- ifelse(nrow(local_genes) > 0, max(local_genes$ymax), 0)
          p_schem <- p_schem + annotate("text", x = s$variant_pos, y = max_y + 1.2, label = aa_str, color = "red", size = 7, fontface = "bold")
        }
      }
      
      p_schem + scale_fill_manual(values=c("TRUE"="orange", "FALSE"="grey82"), guide="none") + 
        theme_minimal() +
        theme(panel.grid.minor = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.title.y = element_blank(),
              axis.title.x = element_text(size = 14),
              axis.text.x = element_text(size = 12)) +
        labs(x=paste0("Chr ", s$chr, " (Mbp)"))
    })
  }
  shinyApp(ui, server)
}
