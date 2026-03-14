# install.packages("ggtext") # Run this if not yet installed
library(data.table)
library(ggplot2)
library(GenomicRanges)
library(patchwork)
library(ggtext) 

# ── 1. RESET GRAPHICS ────────────────────────────────────────────────────────
# This clears any "hung" graphics devices that prevent plots from showing.
while (!is.null(dev.list())) dev.off()

# ── 2. CONFIGURATION ──────────────────────────────────────────────────────────
BASE_DIR    <- "W:/General/Projects2/Jeeyeon Cha/NEW integration MAFA SNPs and DEGs"
SNP_PATH    <- "W:/General/main_directory/variant_tables/B6_and_SJL/B6_SJL_high_confidence_variants_no_csq_info.csv"
MAFA_PATH   <- file.path(BASE_DIR, "MAFA_binding sites with mm10 and mm39 genome coordinates.csv")
GEN_REF_PATH <- file.path(BASE_DIR, "mouse_genes_mm39.csv")

WIN_SIZE <- 100000 
STEP     <- 10000

# mm39 Chromosome lengths
CHR_INFO <- data.table(
  Chr = paste0("Chr", c(1:19, "X")),
  AxisLabel = c(1:19, "X"),
  Length = c(195154279, 181755017, 159745316, 156860686, 151754605, 149582041, 
             144995196, 130127694, 124359700, 130530862, 121973369, 120092757, 
             120885674, 125139656, 104073947, 98008968, 95294699, 90720763, 
             61420004, 169476575)
)

CHR_INFO[, Offset := cumsum(c(0, Length[-.N]))]
CHR_INFO[, EndOffset := cumsum(Length)]
CHR_INFO[, Shade := ifelse(1:.N %% 2 == 0, "#F5F5F5", "#FFFFFF")] 
CHR_INFO[, MidpointGlobal := Offset + (Length / 2)]

GLOBAL_X_LIMITS <- c(0, max(CHR_INFO$EndOffset))

# ── 3. HELPER FUNCTIONS ───────────────────────────────────────────────────────

# Stacks coordinates and preserves LocalPos for zooming
stack_coords <- function(data_dt, pos_col, chr_col="Chr") {
  dt <- copy(data_dt) # Prevent modifying original by reference
  dt[, (chr_col) := as.character(get(chr_col))]
  dt[, (chr_col) := ifelse(grepl("Chr", get(chr_col)), get(chr_col), paste0("Chr", get(chr_col)))]
  
  m <- merge(dt, CHR_INFO[, .(Chr, Offset)], by.x = chr_col, by.y = "Chr")
  if (chr_col != "Chr") setnames(m, chr_col, "Chr")
  
  m[, GlobalPos := get(pos_col) + Offset]
  m[, LocalPos  := get(pos_col)] 
  return(m)
}

# Processes density into a format compatible with both views
prep_density_dt <- function(gr, counts) {
  dt <- as.data.table(gr)[, .(Chr=seqnames, start, end, Count = counts)]
  dt[, Midpoint := (start + end)/2]
  return(stack_coords(dt, "Midpoint"))
}

# Helper for loading DEGs
load_deg_coord <- function(fname, strain_label) {
  p <- file.path(BASE_DIR, fname)
  if(!file.exists(p)) return(NULL)
  d <- fread(p)
  setnames(d, names(d)[1], "GeneId")
  fc_col <- names(d)[ncol(d)]
  m <- merge(d, gene_ref[, .(GeneId, Chr, GlobalPos, LocalPos)], by = "GeneId")
  return(m[, .(Chr, GlobalPos, LocalPos, FC = get(fc_col), Strain = strain_label)])
}

# ── 4. DATA LOADING & PROCESSING ─────────────────────────────────────────────
gene_ref <- fread(GEN_REF_PATH)
setnames(gene_ref, old = c("Gene stable ID", "Chromosome", "Gene start (bp)", "Gene end (bp)", "Strand"), 
         new = c("GeneId", "RawChr", "Start", "End", "StrandDir"))
gene_ref[, TSS := ifelse(StrandDir == 1, Start, End)]
gene_ref <- stack_coords(gene_ref, "TSS", chr_col = "RawChr")

mafa <- fread(MAFA_PATH)
mafa <- stack_coords(mafa, "mm39_start", chr_col = "Chr")

snps <- fread(SNP_PATH, select = c("chr", "pos"))
snps <- stack_coords(snps, "pos", chr_col = "chr")

degs <- rbindlist(list(
  load_deg_coord("DEGs UP MixedSJL Hets.csv", "Mixed"), 
  load_deg_coord("DEGs DOWN MixedSJL Hets.csv", "Mixed"), 
  load_deg_coord("DEGs UP C57 Hets.csv", "C57BL/6J"), 
  load_deg_coord("DEGs DOWN C57 Hets.csv", "C57BL/6J")
))

# Density calculation
seq_lens <- setNames(CHR_INFO$Length, CHR_INFO$Chr)
tiles    <- tileGenome(seq_lens, tilewidth=STEP, cut.last.tile.in.chrom=TRUE)
suppressWarnings({
  start(tiles) <- pmax(1, start(tiles) - 45000)
  end(tiles)   <- end(tiles) + 45000
  tiles        <- trim(tiles)
})

snp_plot_dt  <- prep_density_dt(tiles, countOverlaps(tiles, GRanges(snps$Chr, IRanges(snps$pos, width=1))))
mafa_plot_dt <- prep_density_dt(tiles, countOverlaps(tiles, GRanges(mafa$Chr, IRanges(mafa$mm39_start, width=1))))

# ── 5. PLOTTING THEMES ───────────────────────────────────────────────────────
master_theme <- theme_minimal() + 
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    axis.text.x = element_blank(), 
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 9, color = "black"),
    axis.title.y = element_text(size = 10, angle = 90, face = "bold"),
    plot.subtitle = element_markdown(size = 11, face = "bold", hjust = 0),
    legend.position = "none"
  )

add_background_stripes <- function(plt) {
  for(i in 1:nrow(CHR_INFO)) {
    plt <- plt + annotate("rect", xmin = CHR_INFO$Offset[i], xmax = CHR_INFO$EndOffset[i], 
                          ymin = -Inf, ymax = Inf, fill = CHR_INFO$Shade[i])
  }
  return(plt)
}

# ── 6. GENERATE WHOLE GENOME VIEW ────────────────────────────────────────────
p1 <- add_background_stripes(ggplot()) +
  geom_line(data = snp_plot_dt, aes(x = GlobalPos, y = Count), color = "darkorange", linewidth = 0.3) +
  scale_x_continuous(limits = GLOBAL_X_LIMITS, expand = c(0,0)) +
  labs(y = "SNPs/100kb", subtitle = "A) SNP density") + master_theme

p2 <- add_background_stripes(ggplot()) +
  geom_line(data = mafa_plot_dt, aes(x = GlobalPos, y = Count), color = "forestgreen", linewidth = 0.3) +
  scale_x_continuous(limits = GLOBAL_X_LIMITS, expand = c(0,0)) +
  labs(y = "MafA sites/100Kb", subtitle = "B) Mafa binding site density") + master_theme

p3_sub <- "C) Distribution of DEGs: <span style='color:red;'>Mixed SJL</span> and <span style='color:blue;'>C57BL/6J</span>"
p3 <- add_background_stripes(ggplot()) +
  geom_point(data = degs, aes(x = GlobalPos, y = FC, fill = Strain), 
             size = 2.5, shape = 21, color = "black", stroke = 0.15, alpha = 0.4) +
  scale_fill_manual(values = c("Mixed" = "red", "C57BL/6J" = "blue")) +
  scale_x_continuous(limits = GLOBAL_X_LIMITS, expand = c(0,0),
                     breaks = CHR_INFO$MidpointGlobal, labels = CHR_INFO$AxisLabel) +
  labs(y = "DEGs (log2 FC)", subtitle = p3_sub, x = "Chromosome") +
  master_theme +
  theme(axis.text.x = element_text(size = 10, color = "black"),
        axis.title.x = element_text(face = "bold", size = 12, margin = margin(t=10)))

genome_plot <- (p1 / p2 / p3) + plot_layout(heights = c(1, 1, 1))

# Force output
print(genome_plot)
ggsave(file.path(BASE_DIR, "final_manhattan_GENOME.png"), 
       plot = genome_plot, width = 18, height = 9, dpi = 300, bg = "white")

# ── 7. ZOOM FUNCTION ─────────────────────────────────────────────────────────
create_zoomed_view <- function(target_chr) {
  # Local Filters
  s_zoom <- snp_plot_dt[Chr == target_chr]
  m_zoom <- mafa_plot_dt[Chr == target_chr]
  d_zoom <- degs[Chr == target_chr]
  c_len  <- CHR_INFO[Chr == target_chr, Length]
  
  # Standard Zoom Theme (Mb labels active)
  zoom_theme <- master_theme + 
    theme(axis.text.x = element_text(size = 10, color = "black"),
          axis.title.x = element_text(face = "bold", size = 12))
  
  # Panel A
  z1 <- ggplot() + 
    annotate("rect", xmin=0, xmax=c_len, ymin=-Inf, ymax=Inf, fill="white") +
    geom_line(data = s_zoom, aes(x = LocalPos, y = Count), color = "navy", linewidth = 0.6) +
    scale_x_continuous(limits = c(0, c_len), expand = c(0,0), labels = function(x) x/1e6) +
    labs(y = "SNPs/100kb", subtitle = paste("A) SNP Density -", target_chr)) + zoom_theme
  
  # Panel B
  z2 <- ggplot() + 
    annotate("rect", xmin=0, xmax=c_len, ymin=-Inf, ymax=Inf, fill="white") +
    geom_line(data = m_zoom, aes(x = LocalPos, y = Count), color = "forestgreen", linewidth = 0.6) +
    scale_x_continuous(limits = c(0, c_len), expand = c(0,0), labels = function(x) x/1e6) +
    labs(y = "MafA sites/100Kb", subtitle = "B) MafA Binding Sites") + zoom_theme
  
  # Panel C
  z3 <- ggplot() + 
    annotate("rect", xmin=0, xmax=c_len, ymin=-Inf, ymax=Inf, fill="white") +
    geom_point(data = d_zoom, aes(x = LocalPos, y = FC, fill = Strain), 
               size = 4, shape = 21, color = "black", stroke = 0.3, alpha = 0.7) +
    scale_fill_manual(values = c("Mixed" = "red", "C57BL/6J" = "blue")) +
    scale_x_continuous(limits = c(0, c_len), expand = c(0,0), labels = function(x) x/1e6) +
    labs(y = "log2 FC", x = paste(target_chr, "Position (Mb)"), 
         subtitle = paste0("C) DEGs: <span style='color:red;'>Mixed SJL</span> and <span style='color:blue;'>C57BL/6J</span>")) +
    zoom_theme
  
  zoom_plot <- (z1 / z2 / z3) + plot_layout(heights = c(1, 1, 1))
  
  print(zoom_plot)
  ggsave(file.path(BASE_DIR, paste0("Illustrative_Zoom_", target_chr, ".png")), 
         plot = zoom_plot, width = 12, height = 10, dpi = 300, bg = "white")
}

# ── 8. EXECUTION ─────────────────────────────────────────────────────────────
create_zoomed_view("Chr4")
