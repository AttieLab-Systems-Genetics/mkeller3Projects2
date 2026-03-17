library(shiny)
library(data.table)
library(ggplot2)
library(GenomicRanges)

# ── 1. DATA ENGINE ────────────────────────────────────────────────────────
prepare_data <- function() {
  base_path <- "W:/General/Projects2/Jeeyeon Cha/NEW integration MAFA SNPs and DEGs"
  
  chr_map <- data.table(
    Chr = paste0("Chr", c(1:19, "X", "Y")),
    Length = as.numeric(c(195154279, 181755017, 159745316, 156860686, 151754605, 149582041, 
                          144995196, 130127694, 124359700, 130530862, 121973369, 120092757, 
                          120885674, 125139656, 104073947, 98008968, 95294699, 90720763, 
                          61420004, 169476575, 91455967))
  )
  chr_map[, Offset := cumsum(data.table::shift(Length, fill = 0, type = "lag"))]
  chr_map[, Midpoint := Offset + (Length / 2)]
  
  master    <- fread(file.path(base_path, "Master_DEG_Strain_Comparison.csv"))
  mafa      <- fread(file.path(base_path, "MafA_Peaks_with_SNPs_V2.csv"))
  all_genes <- fread(file.path(base_path, "mouse_genes_mm39.csv"))
  
  mafa[chr_map, on = "Chr", GlobalPos := Mid + i.Offset]
  return(list(degs = master, mafa = mafa, genes = all_genes, chr_map = chr_map))
}

# ── 2. UI ──────────────────────────────────────────────────────────────────
ui <- fluidPage(
  tags$head(tags$style(HTML("
    .hover-tooltip { 
      position: absolute; z-index: 1000; background: rgba(255, 255, 255, 0.98); 
      border: 2px solid #333; padding: 10px; border-radius: 6px; 
      box-shadow: 3px 3px 12px rgba(0,0,0,0.2); pointer-events: none; font-size: 12px;
    }
    .well-meta { background: #fffdf5; border: 1.5px solid gold; padding: 12px; margin-top: 15px; }
    .meta-title { font-weight: bold; font-size: 1.1em; border-bottom: 1px solid #ddd; margin-bottom: 8px; color: #856404; }
    .btn-rezoom { background-color: #007bff; color: white; font-weight: bold; width: 100%; margin-bottom: 10px; }
    .deg-item { margin-bottom: 4px; font-weight: bold; font-size: 0.92em; line-height: 1.2; }
  "))),
  
  titlePanel(span("MafA Discovery: Integrated Genomic Explorer", style="font-weight:bold;")),
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      actionButton("reset_view", "↺ Reset to Genome-Wide", class = "btn-rezoom"),
      selectizeInput("search_gene", "Search Gene Symbol:", choices = NULL),
      selectInput("zoom_mode", "View Mode:", choices = c("Genome-Wide", "Chromosome", "Locus Zoom")),
      uiOutput("chr_selector_ui"),
      hr(),
      checkboxGroupInput("show_cat", "Visible Categories:", 
                         choices = c("Shared", "C57_Specific", "SJL_Specific", "Discordant"),
                         selected = c("Shared", "C57_Specific", "SJL_Specific", "Discordant")),
      numericInput("win_kb", "Locus Window (Kbp):", value = 250),
      sliderInput("min_score", "Min Peak Score:", 0, 10000, 0),
      uiOutput("metadata_panel")
    ),
    mainPanel(
      width = 9,
      div(style = "position:relative",
          plotOutput("manhattan", height = "480px", click = "p_click", 
                     hover = hoverOpts("p_hover", delay = 50, delayType = "debounce")),
          uiOutput("hover_info_manhattan")),
      hr(),
      div(style = "position:relative",
          plotOutput("schematic", height = "400px", 
                     hover = hoverOpts("s_hover", delay = 50)),
          uiOutput("hover_info_schematic"))
    )
  )
)

# ── 3. SERVER ──────────────────────────────────────────────────────────────
server <- function(input, output, session) {
  d <- prepare_data()
  v <- reactiveValues(active_pk = NULL)
  cat_colors <- c("Shared"="#228B22", "C57_Specific"="#0000CC", "SJL_Specific"="#CC0000", "Discordant"="#FF8C00", "Non-DE"="#D3D3D3")
  
  updateSelectizeInput(session, "search_gene", choices = sort(unique(d$genes$Symbol)), server = TRUE)
  output$chr_selector_ui <- renderUI({ selectInput("sel_chr", "Select Chromosome:", choices = d$chr_map$Chr) })
  
  observeEvent(input$reset_view, { v$active_pk <- NULL; updateSelectInput(session, "zoom_mode", selected = "Genome-Wide") })
  
  filtered_peaks <- reactive({
    req(input$show_cat)
    pks <- d$mafa[`Peak Score` >= input$min_score]
    p_gr  <- GRanges(pks$Chr, IRanges(as.integer(pks$Mid), width=1))
    d_gr  <- GRanges(d$degs$Chr, IRanges(as.integer(d$degs$Start), width=1))
    hits  <- findOverlaps(p_gr, d_gr, maxgap = input$win_kb * 1000)
    if(length(hits) == 0) return(NULL)
    res <- pks[unique(queryHits(hits))]
    idx <- nearest(GRanges(res$Chr, IRanges(as.integer(res$Mid), width=1)), d_gr)
    res[, Category := d$degs$Category[idx]]
    res[, Direction := d$degs$Direction[idx]]
    res[Category %in% input$show_cat]
  })
  
  observeEvent(input$p_click, {
    req(filtered_peaks())
    sel <- nearPoints(filtered_peaks(), input$p_click, xvar="GlobalPos", yvar="Peak Score", threshold=10, maxpoints=1)
    if(nrow(sel) > 0) v$active_pk <- sel
  })
  
  # --- MANHATTAN HOVER ---
  output$hover_info_manhattan <- renderUI({
    hover <- input$p_hover; req(hover, filtered_peaks())
    pt <- nearPoints(filtered_peaks(), hover, xvar="GlobalPos", yvar="Peak Score", threshold=10, maxpoints=1)
    if(nrow(pt) == 0) return(NULL)
    
    win_size <- input$win_kb * 1000
    local_degs <- d$degs[Chr == pt$Chr & abs(Start - pt$Mid) <= win_size]
    deg_counts <- if(nrow(local_degs) > 0) local_degs[, .N, by=Category] else data.table(Category=character(), N=numeric())
    
    div(class = "hover-tooltip", style = paste0("left:", hover$coords_css$x + 15, "px; top:", hover$coords_css$y + 15, "px; border-color:", cat_colors[pt$Category]),
        p(strong(paste("MafA peak", sub(".*peak_", "", pt$PeakID)))),
        p(strong("SNPs: "), pt$SNP_Count),
        hr(style="margin:5px 0"),
        p(em("Window DEGs:")),
        if(nrow(deg_counts) > 0) {
          tags$ul(style="padding-left:15px; margin:0",
                  lapply(1:nrow(deg_counts), function(i) tags$li(paste0(deg_counts$Category[i], ": ", deg_counts$N[i]))))
        } else { p("None", style="font-size:10px") }
    )
  })
  
  output$manhattan <- renderPlot({
    df <- filtered_peaks()
    
    # Define X range
    x_r <- if(input$zoom_mode == "Genome-Wide") {
      c(0, max(d$chr_map$Offset + d$chr_map$Length))
    } else if(input$zoom_mode == "Chromosome") {
      req(input$sel_chr); r <- d$chr_map[Chr == input$sel_chr]; c(r$Offset, r$Offset + r$Length)
    } else {
      req(v$active_pk); c(v$active_pk$GlobalPos - 5e5, v$active_pk$GlobalPos + 5e5)
    }
    
    # AUTO-ZOOM Y-AXIS
    y_max <- 500 
    if(!is.null(df)) {
      visible_pts <- df[GlobalPos >= x_r[1] & GlobalPos <= x_r[2]]
      if(nrow(visible_pts) > 0) y_max <- max(visible_pts$`Peak Score`) * 1.15
    }
    
    p <- ggplot() +
      geom_rect(data = d$chr_map, aes(xmin=Offset, xmax=Offset+Length, ymin=-Inf, ymax=Inf), 
                fill=rep(c("white", "#f9f9f9"), length.out=21), show.legend=FALSE) +
      scale_x_continuous(limits=x_r, breaks=d$chr_map$Midpoint, labels=gsub("Chr","",d$chr_map$Chr), expand=c(0,0)) +
      scale_y_continuous(limits=c(0, y_max), expand=expansion(mult=c(0.02, 0))) +
      theme_minimal() + labs(x="Chromosome", y="MafA Peak Score") +
      theme(axis.line = element_line(color="black", linewidth=0.3), panel.grid = element_blank())
    
    if(!is.null(df)) {
      p <- p + geom_point(data=df, aes(x=GlobalPos, y=`Peak Score`, fill=Category, shape=Direction, size=SNP_Count), 
                          color="black", stroke=0.3, alpha=0.7) +
        scale_fill_manual(values=cat_colors) + 
        scale_shape_manual(values=c("UP"=24, "DOWN"=25))
    }
    
    if(!is.null(v$active_pk)) {
      p <- p + geom_point(data=v$active_pk, aes(x=GlobalPos, y=`Peak Score`), 
                          shape=23, fill="gold", size=7, stroke=1.2)
    }
    p
  })
  
  output$schematic <- renderPlot({
    req(v$active_pk); pk <- v$active_pk; win <- input$win_kb * 1000
    window_genes <- d$genes[Chr == pk$Chr & abs(Start - pk$Mid) <= win]
    if(nrow(window_genes) == 0) return(NULL)
    
    plot_df <- merge(window_genes, d$degs[, .(GeneId, Category)], by="GeneId", all.x=TRUE)
    plot_df[is.na(Category), Category := "Non-DE"]
    plot_df <- plot_df[Category != "Non-DE" | `Gene type` == "protein_coding"]
    if(nrow(plot_df) == 0) return(NULL)
    
    plot_df <- plot_df[order(Start)][, y_lev := 1:.N]
    plot_df[, tss := ifelse(Strand == 1, Start, End)]
    plot_df[, tip := ifelse(Strand == 1, tss + (win * 0.025), tss - (win * 0.025))]
    
    ggplot(plot_df) +
      geom_rect(aes(xmin=pk$Start, xmax=pk$End, ymin=-Inf, ymax=Inf), fill="gold", alpha=0.4) +
      geom_vline(xintercept = pk$Mid, color="gold4", linetype="dashed") +
      geom_rect(aes(xmin=Start, xmax=End, ymin=y_lev-0.3, ymax=y_lev+0.3, fill=Category), color="black", linewidth=0.2) +
      geom_segment(data=plot_df[Category != "Non-DE"], aes(x=tss, xend=tip, y=y_lev, yend=y_lev), arrow=arrow(length=unit(0.45, "cm"), type="closed"), color="black", linewidth=0.6) +
      geom_text(aes(x=(Start+End)/2, y=y_lev+0.6, label=Symbol, color=Category), fontface="bold", size=5) +
      scale_fill_manual(values=cat_colors) + scale_color_manual(values=cat_colors) +
      scale_x_continuous(labels = function(x) format(x/1e6, digits=5)) +
      labs(x=paste(pk$Chr, "(Mbp)"), y="") +
      theme_minimal() + theme(axis.line.x = element_line(color="black", linewidth=0.3), axis.text.y=element_blank(), panel.grid=element_blank())
  })
  
  output$metadata_panel <- renderUI({
    req(v$active_pk); pk <- v$active_pk
    win_size <- input$win_kb * 1000
    local_degs <- d$degs[Chr == pk$Chr & abs(Start - pk$Mid) <= win_size]
    
    wellPanel(class = "well-meta",
              div(class = "meta-title", paste("Peak", sub(".*peak_", "", pk$PeakID))),
              p(strong("Location: "), round(pk$Mid/1e6, 3), " Mbp"),
              p(strong("SNPs: "), pk$SNP_Count),
              hr(),
              p(strong("Locus DEGs & log2FC:")),
              if(nrow(local_degs) > 0) {
                tagList(lapply(1:nrow(local_degs), function(i) {
                  row <- local_degs[i]
                  fc_text <- if(row$Category == "Shared" | row$Category == "Discordant") {
                    paste0(" (C57:", round(row$log2FC_C57, 2), ", SJL:", round(row$log2FC_SJL, 2), ")")
                  } else if(row$Category == "C57_Specific") {
                    paste0(" (C57:", round(row$log2FC_C57, 2), ")")
                  } else {
                    paste0(" (SJL:", round(row$log2FC_SJL, 2), ")")
                  }
                  div(class="deg-item", style=paste0("color:", cat_colors[row$Category]), paste0("• ", row$Symbol, fc_text))
                }))
              } else { p("No DEGs in window.", style="font-style:italic") }
    )
  })
}

shinyApp(ui, server)

