library(shiny)
library(data.table)
library(ggplot2)
library(GenomicRanges)

# ── 1. DATA ENGINE ────────────────────────────────────────────────────────
prepare_integrated_data <- function(base_dir) {
  chr_map <- data.table(
    Chr = paste0("Chr", c(1:19, "X", "Y")),
    Length = as.numeric(c(195154279, 181755017, 159745316, 156860686, 151754605, 149582041, 
                          144995196, 130127694, 124359700, 130530862, 121973369, 120092757, 
                          120885674, 125139656, 104073947, 98008968, 95294699, 90720763, 
                          61420004, 169476575, 91455967))
  )
  chr_map[, OffsetVal := cumsum(data.table::shift(Length, fill = 0, type = "lag"))]
  chr_map[, MidpointGlobal := OffsetVal + (Length / 2)]
  
  pc_genes <- fread(file.path(base_dir, "mouse_genes_mm39.csv"))
  master   <- fread(file.path(base_dir, "Master_DEG_Strain_Comparison.csv"))
  mafa     <- fread(file.path(base_dir, "MafA_Peaks_with_SNPs_V2.csv"))
  
  mafa[chr_map, on = "Chr", Offset_Int := i.OffsetVal]
  mafa[, GlobalPos := Mid + Offset_Int]
  
  return(list(degs = master, pc_genes = pc_genes[`Gene type`=="protein_coding"], mafa = mafa, chr_map = chr_map))
}

# ── 2. UI ──────────────────────────────────────────────────────────────────
ui <- fluidPage(
  titlePanel(span("MafA Discovery: Integrated Genomic Explorer", style = "font-weight: bold;")),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Search & Navigation"),
      selectizeInput("search_gene", "Search Gene Symbol:", choices = NULL),
      selectInput("zoom_mode", "View Mode:", choices = c("Genome-Wide", "Chromosome", "Locus Zoom")),
      uiOutput("chr_select_ui"),
      hr(),
      h4("Filters"),
      checkboxGroupInput("filter_cat", "Visible Categories:", 
                         choices = c("Shared", "C57_Specific", "SJL_Specific", "Discordant", "Non-DE"),
                         selected = c("Shared", "C57_Specific", "SJL_Specific", "Discordant")),
      numericInput("win_kb", "Search Window (Kbp):", value = 250), # Updated default
      sliderInput("min_score", "Min Peak Score:", min = 0, max = 10000, value = 0), # Updated default
      hr(),
      uiOutput("metadata_panel")
    ),
    mainPanel(
      width = 9,
      plotOutput("manhattan", height = "500px", click = "plot_click"),
      hr(),
      plotOutput("schematic", height = "400px")
    )
  )
)

# ── 3. SERVER ──────────────────────────────────────────────────────────────
server <- function(input, output, session) {
  
  data_core <- reactive({ 
    prepare_integrated_data("W:/General/Projects2/Jeeyeon Cha/NEW integration MAFA SNPs and DEGs") 
  })
  
  cat_colors <- c("Shared"="#228B22", "C57_Specific"="#0000CC", "SJL_Specific"="#CC0000", "Discordant"="#FF8C00", "Non-DE"="#D3D3D3")
  
  output$chr_select_ui <- renderUI({
    req(data_core()); selectInput("zoom_chr", "Select Chromosome:", choices = data_core()$chr_map$Chr)
  })
  
  observeEvent(data_core(), {
    updateSelectizeInput(session, "search_gene", choices = sort(unique(data_core()$pc_genes$Symbol)), server = TRUE)
  })
  
  v <- reactiveValues(active_peak = NULL)
  
  observeEvent(input$plot_click, {
    req(peaks_in_view())
    sel <- nearPoints(peaks_in_view(), input$plot_click, xvar = "GlobalPos", yvar = "Peak Score", threshold = 10, maxpoints = 1)
    if(nrow(sel) > 0) v$active_peak <- sel
  })
  
  peaks_in_view <- reactive({
    req(data_core(), input$filter_cat)
    d <- data_core()
    mafa_f <- d$mafa[`Peak Score` >= input$min_score]
    if(nrow(mafa_f) == 0) return(NULL)
    
    mafa_gr <- GRanges(mafa_f$Chr, IRanges(as.integer(mafa_f$Mid), width = 1))
    deg_gr  <- GRanges(d$degs$Chr, IRanges(as.integer(d$degs$Start), width = 1))
    
    hits <- findOverlaps(mafa_gr, deg_gr, maxgap = input$win_kb * 1000)
    if(length(hits) == 0) return(NULL)
    
    out <- copy(mafa_f[unique(queryHits(hits))])
    out_gr <- GRanges(out$Chr, IRanges(as.integer(out$Mid), width = 1))
    near_idx <- nearest(out_gr, deg_gr)
    
    out[, Category := d$degs$Category[near_idx]]
    out[, Direction := d$degs$Direction[near_idx]]
    return(out[Category %in% input$filter_cat])
  })
  
  # MANHATTAN PLOT
  output$manhattan <- renderPlot({
    req(data_core())
    chr_info <- data_core()$chr_map
    df <- peaks_in_view()
    
    x_range <- if(input$zoom_mode == "Genome-Wide") {
      c(0, max(chr_info$OffsetVal + chr_info$Length))
    } else if(input$zoom_mode == "Chromosome") {
      req(input$zoom_chr); r <- chr_info[Chr == input$zoom_chr]; c(r$OffsetVal, r$OffsetVal + r$Length)
    } else {
      req(v$active_peak); c(v$active_peak$GlobalPos - 4e5, v$active_peak$GlobalPos + 4e5)
    }
    
    y_max <- if(!is.null(df) && nrow(df[GlobalPos >= x_range[1] & GlobalPos <= x_range[2]]) > 0) {
      max(df[GlobalPos >= x_range[1] & GlobalPos <= x_range[2]]$`Peak Score`, na.rm = TRUE) * 1.1
    } else { 1000 }
    
    p <- ggplot() +
      geom_rect(data = chr_info, aes(xmin = OffsetVal, xmax = OffsetVal + Length, ymin = -Inf, ymax = Inf), 
                fill = rep(c("#FFFFFF", "#F9F9F9"), length.out = nrow(chr_info)), show.legend = FALSE) +
      scale_x_continuous(limits = x_range, breaks = chr_info$MidpointGlobal, labels = gsub("Chr", "", chr_info$Chr), expand = c(0,0)) +
      scale_y_continuous(limits = c(0, y_max)) +
      theme_minimal() +
      theme(panel.grid = element_blank(), 
            axis.line = element_line(color = "black"),
            axis.title = element_text(size = 14),
            axis.text = element_text(size = 14)) +
      labs(x = "Chromosome", y = "MafA Peak Score")
    
    if(!is.null(df) && nrow(df) > 0) {
      p <- p + geom_point(data = df, aes(x = GlobalPos, y = `Peak Score`, fill = Category, shape = Direction, size = SNP_Count), 
                          color = "black", stroke = 0.3, alpha = 0.8) +
        scale_fill_manual(values = cat_colors, drop = FALSE) +
        scale_shape_manual(values = c("UP" = 24, "DOWN" = 25))
    }
    
    if(!is.null(v$active_peak)) {
      p <- p + geom_point(data = v$active_peak, aes(x = GlobalPos, y = `Peak Score`), 
                          fill = "gold", shape = 23, size = 8, color = "black", stroke = 1.2)
    }
    p
  })
  
  # GENOMIC SCHEMATIC
  output$schematic <- renderPlot({
    req(v$active_peak, data_core())
    pk <- v$active_peak
    win <- input$win_kb * 1000
    
    genes_win <- data_core()$pc_genes[Chr == pk$Chr & abs(Start - pk$Mid) <= win]
    if(nrow(genes_win) == 0) return(NULL)
    
    plot_genes <- merge(genes_win, data_core()$degs[, .(GeneId, Category)], by = "GeneId", all.x = TRUE)
    plot_genes[is.na(Category), Category := "Non-DE"]
    plot_genes <- plot_genes[order(Start)][, y_lev := 1:.N * -1]
    
    ggplot() +
      geom_rect(aes(xmin = pk$Start, xmax = pk$End, ymin = -Inf, ymax = Inf), fill = "#FFF9C4", alpha = 0.5) +
      geom_vline(xintercept = pk$Mid, linetype = "dotted", color = "gold4", linewidth = 1) +
      geom_rect(data = plot_genes, aes(xmin = Start, xmax = End, ymin = y_lev - 0.3, ymax = y_lev + 0.3, fill = Category), 
                color = "black", linewidth = 0.1) +
      # RESTORED COLOR FOR TEXT LABELS
      geom_text(data = plot_genes, aes(x = (Start+End)/2, y = y_lev + 0.6, label = Symbol, color = Category), 
                fontface = "bold", size = 5, show.legend = FALSE) +
      scale_fill_manual(values = cat_colors, drop = FALSE) +
      scale_color_manual(values = cat_colors, drop = FALSE) +
      scale_x_continuous(labels = function(x) x/1e6) + 
      labs(x = paste(pk$Chr, "(Mbp)"), y = "", title = "Genomic Locus Detail") +
      theme_minimal() + 
      theme(axis.text.y = element_blank(), 
            axis.text.x = element_text(size = 14), # Standardized
            axis.title.x = element_text(size = 14), # Standardized
            panel.grid = element_blank(), 
            axis.line.x = element_line(color = "black")) 
  })
  
  # METADATA PANEL
  output$metadata_panel <- renderUI({
    req(v$active_peak)
    pk <- v$active_peak
    degs <- data_core()$degs[Chr == pk$Chr & abs(Start - pk$Mid) <= input$win_kb * 1000]
    
    wellPanel(
      style = "background: #fffdf5; border: 1.2px solid gold; padding: 10px;",
      h5("SELECTED PEAK", style="font-weight:bold; color:#856404;"),
      p(strong("ID: "), pk$PeakID),
      p(strong("Range: "), pk$Chr, ":", format(pk$Start, big.mark=","), "-", format(pk$End, big.mark=",")),
      p(strong("Score: "), pk$`Peak Score`, " | ", strong("SNPs: "), pk$SNP_Count),
      hr(),
      h5("PROXIMAL DEGs", style="font-weight:bold;"),
      if(nrow(degs) > 0) {
        tags$ul(lapply(1:nrow(degs), function(i) {
          tags$li(
            span(paste0(degs$Symbol[i], ": "), style=paste0("color:", cat_colors[degs$Category[i]], "; font-weight:bold;")),
            span(sprintf("C57 FC = %.2f | SJL FC = %.2f", degs$log2FC_C57[i], degs$log2FC_SJL[i]), style="font-size:0.85em;")
          )
        }))
      } else { p("No DEGs in window.") }
    )
  })
}

shinyApp(ui, server)


