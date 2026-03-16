# install.packages("ggtext") # Run this if not yet installed
library(data.table)
library(ggplot2)
library(GenomicRanges)
library(patchwork)
library(ggtext) 

# ── 1. RESET GRAPHICS ────────────────────────────────────────────────────────
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
stack_coords <- function(data_dt, pos_col, chr_col="Chr") {
  dt <- copy(data_dt)
  dt[, (chr_col) := as.character(get(chr_col))]
  dt[, (chr_col) := ifelse(grepl("Chr", get(chr_col)), get(chr_col), paste0("Chr", get(chr_col)))]
  m <- merge(dt, CHR_INFO[, .(Chr, Offset)], by.x = chr_col, by.y = "Chr")
  if (chr_col != "Chr") setnames(m, chr_col, "Chr")
  m[, GlobalPos := get(pos_col) + Offset]
  m[, LocalPos  := get(pos_col)] 
  return(m)
}

prep_density_dt <- function(gr, counts) {
  dt <- as.data.table(gr)[, .(Chr=seqnames, start, end, Count = counts)]
  dt[, Midpoint := (start + end)/2]
  return(stack_coords(dt, "Midpoint"))
}

load_deg_coord <- function(fname, strain_label) {
  p <- file.path(BASE_DIR, fname)
  if(!file.exists(p)) return(NULL)
  d <- fread(p)
  setnames(d, names(d)[1], "GeneId")
  fc_col <- names(d)[ncol(d)]
  m <- merge(d, gene_ref[, .(GeneId, Chr, GlobalPos, LocalPos)], by = "GeneId")
  return(m[, .(Chr, GlobalPos, LocalPos, FC = get(fc_col), Strain = strain_label)])
}

# ── 4. DATA LOADING ──────────────────────────────────────────────────────────
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

seq_lens <- setNames(CHR_INFO$Length, CHR_INFO$Chr)
tiles    <- tileGenome(seq_lens, tilewidth=STEP, cut.last.tile.in.chrom=TRUE)
suppressWarnings({
  start(tiles) <- pmax(1, start(tiles) - 45000); end(tiles) <- end(tiles) + 45000; tiles <- trim(tiles)
})

snp_plot_dt  <- prep_density_dt(tiles, countOverlaps(tiles, GRanges(snps$Chr, IRanges(snps$pos, width=1))))
mafa_plot_dt <- prep_density_dt(tiles, countOverlaps(tiles, GRanges(mafa$Chr, IRanges(mafa$mm39_start, width=1))))

# ── 5. THEMES ────────────────────────────────────────────────────────────────
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

zoom_theme <- master_theme + 
  theme(axis.text.x = element_text(size = 10, color = "black"),
        axis.title.x = element_text(face = "bold", size = 12))

add_background_stripes <- function(plt) {
  for(i in 1:nrow(CHR_INFO)) {
    plt <- plt + annotate("rect", xmin = CHR_INFO$Offset[i], xmax = CHR_INFO$EndOffset[i], 
                          ymin = -Inf, ymax = Inf, fill = CHR_INFO$Shade[i])
  }
  return(plt)
}

# ── 6. PLOTTING GENOME ───────────────────────────────────────────────────────
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
print(genome_plot)
ggsave(file.path(BASE_DIR, "final_manhattan_GENOME.png"), plot = genome_plot, width = 18, height = 9, dpi = 300, bg = "white")

# ── 7. ZOOM FUNCTION (CLEANED X-AXIS) ────────────────────────────────────────
create_zoomed_view <- function(target_chr) {
  s_zoom <- snp_plot_dt[Chr == target_chr]
  m_zoom <- mafa_plot_dt[Chr == target_chr]
  d_zoom <- degs[Chr == target_chr]
  c_len  <- CHR_INFO[Chr == target_chr, Length]
  
  # Top Panel (SNPs) - Uses master_theme to hide X labels
  z1 <- ggplot() + 
    annotate("rect", xmin=0, xmax=c_len, ymin=-Inf, ymax=Inf, fill="white") +
    geom_line(data = s_zoom, aes(x = LocalPos, y = Count), color = "navy", linewidth = 0.6) +
    scale_x_continuous(limits = c(0, c_len), expand = c(0,0)) +
    labs(y = "SNPs/100kb", subtitle = paste("A) SNP Density -", target_chr)) + 
    master_theme
  
  # Middle Panel (MafA) - Uses master_theme to hide X labels
  z2 <- ggplot() + 
    annotate("rect", xmin=0, xmax=c_len, ymin=-Inf, ymax=Inf, fill="white") +
    geom_line(data = m_zoom, aes(x = LocalPos, y = Count), color = "forestgreen", linewidth = 0.6) +
    scale_x_continuous(limits = c(0, c_len), expand = c(0,0)) +
    labs(y = "MafA sites/100Kb", subtitle = "B) MafA Binding Sites") + 
    master_theme
  
  # Bottom Panel (DEGs)
  z3_sub <- paste0("C) DEGs: <span style='color:red;'>Mixed SJL</span> and <span style='color:blue;'>C57BL/6J</span>")
  z3 <- ggplot() + 
    annotate("rect", xmin=0, xmax=c_len, ymin=-Inf, ymax=Inf, fill="white") +
    # Larger dashes for the y=0 baseline
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
    geom_point(data = d_zoom, aes(x = LocalPos, y = FC, fill = Strain), 
               size = 4, shape = 21, color = "black", stroke = 0.3, alpha = 0.7) +
    scale_fill_manual(values = c("Mixed" = "red", "C57BL/6J" = "blue")) +
    scale_x_continuous(limits = c(0, c_len), expand = c(0,0), labels = function(x) x/1e6) +
    labs(y = "log2 FC", x = paste(target_chr, "Position (Mb)"), subtitle = z3_sub) +
    zoom_theme
  
  zoom_plot <- (z1 / z2 / z3) + plot_layout(heights = c(1, 1, 1))
  print(zoom_plot)
  ggsave(file.path(BASE_DIR, paste0("Illustrative_Zoom_", target_chr, ".png")), 
         plot = zoom_plot, width = 18, height = 9, dpi = 600, bg = "white")
}

# ── 8. EXECUTION ─────────────────────────────────────────────────────────────
create_zoomed_view("Chr1")
create_zoomed_view("Chr2")
create_zoomed_view("Chr3")
create_zoomed_view("Chr4")
create_zoomed_view("Chr5")
create_zoomed_view("Chr6")
create_zoomed_view("Chr7")
create_zoomed_view("Chr8")
create_zoomed_view("Chr9")
create_zoomed_view("Chr11")
create_zoomed_view("Chr12")
create_zoomed_view("Chr13")
create_zoomed_view("Chr14")
create_zoomed_view("Chr15")
create_zoomed_view("Chr16")
create_zoomed_view("Chr17")
create_zoomed_view("Chr18")
create_zoomed_view("Chr19")


###----------------------------------------------------
# Next vignette occurs here
###----------------------------------------------------


###----------------------------------------------------
# Vignette: MafA Binding Site Length Distribution
###----------------------------------------------------

# 1. Create the Histogram
# Using the pre-existing 'length' column in the mafa data table
p_mafa_hist <- ggplot(mafa, aes(x = length)) +
  geom_histogram(
    fill = "forestgreen", 
    color = "white", 
    bins = 50, 
    alpha = 0.9
  ) +
  # Use scale_x_log10() if you have a very wide range of lengths, 
  # but standard linear scale is often clearer for ChIP-seq peak widths.
  scale_x_continuous(expand = c(0.02, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(
    title = "Distribution of MafA Binding Site Lengths",
    subtitle = "Vignette: Genomic footprint analysis of MafA binding sites (mm39)",
    x = "Peak Length (bp)",
    y = "Frequency (Count)"
  ) +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(face = "bold", size = 11),
    axis.text = element_text(color = "black")
  )

# 2. Output and Save
print(p_mafa_hist)

ggsave(
  filename = file.path(BASE_DIR, "Vignette_MafA_Site_Lengths.png"),
  plot = p_mafa_hist,
  width = 12, 
  height = 6, 
  dpi = 600, 
  bg = "white"
)

# 3. Summary Statistics Console Log
message("--- MafA Site Length Statistics ---")
print(mafa[, .(
  Mean = mean(length), 
  Median = median(length), 
  Min = min(length), 
  Max = max(length)
)])



###----------------------------------------------------
# Vignette: MafA Distance to TSS (Mbp/Kb Refinement)
###----------------------------------------------------

# 1. Clean data: Ensure "Distance to TSS" is numeric
mafa[, Dist_TSS_Num := as.numeric(`Distance to TSS`)]

# 2. Main Plot: -1 to +1 Mbp Range
p_main_tss <- ggplot(mafa, aes(x = Dist_TSS_Num)) +
  geom_histogram(
    binwidth = 25000, 
    fill = "forestgreen", 
    color = "white",
    alpha = 0.8
  ) +
  # Adding vertical and horizontal lines at the origin
  geom_vline(xintercept = 0, color = "black", linewidth = 0.6) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.6) +
  coord_cartesian(xlim = c(-1000000, 1000000)) +
  scale_x_continuous(labels = function(x) paste0(x/1e6, " Mb")) +
  labs(
    title = "MafA Genomic Positioning Relative to TSS",
    subtitle = "Window: +/- 1 Mbp | Bin Width: 25 kb",
    x = "Distance to TSS",
    y = "Frequency"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(), # Remove grid to focus on the axis lines
    axis.line = element_line(color = "black"), # Explicit axis lines
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(face = "bold", size = 11),
  )

# 3. Inset Plot: -1 to +1 Kb Range
p_inset_tss <- ggplot(mafa, aes(x = Dist_TSS_Num)) +
  geom_histogram(
    binwidth = 25, 
    fill = "forestgreen", 
    color = "white"
  ) +
  # Adding vertical and horizontal lines at the origin
  geom_vline(xintercept = 0, color = "black", linewidth = 0.5) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  coord_cartesian(xlim = c(-1000, 1000)) +
  scale_x_continuous(labels = function(x) paste0(x/1000, " kb")) +
  labs(title = NULL, x = "Distance to TSS", y = "Frequency") +
  theme_minimal() +
  theme(
    plot.background = element_rect(fill = "white", color = "white", linewidth = 0.8),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text = element_text(size = 9, color = "black")
  )

# 4. Combine with Patchwork (Upper Right Placement)
# Using coordinates to anchor the inset in the top-right quadrant
final_tss_vignette <- p_main_tss + 
  inset_element(
    p_inset_tss, 
    left = 0.60,   # Position relative to plot area (0 to 1)
    bottom = 0.25, 
    right = 0.98, 
    top = 0.98,
    align_to = 'panel'
  )

# 5. Output and Save
print(final_tss_vignette)

ggsave(
  filename = file.path(BASE_DIR, "Vignette_MafA_TSS_Structural.png"),
  plot = final_tss_vignette,
  width = 12, 
  height = 6, 
  dpi = 600, 
  bg = "white"
)








###----------------------------------------------------
# Vignette: SNP-Containing MafA Binding Sites
###----------------------------------------------------

# 1. Convert Data to GRanges for High-Speed Overlap
# We use mm39 coordinates for the intersection
mafa_gr <- GRanges(
  seqnames = mafa$Chr,
  ranges   = IRanges(start = mafa$mm39_start, end = mafa$mm39_end),
  mcols    = mafa # Keep all original metadata
)

snp_gr <- GRanges(
  seqnames = snps$Chr,
  ranges   = IRanges(start = snps$LocalPos, width = 1)
)

# 2. Identify Overlaps
# findOverlaps returns a mapping of which SNP is in which Peak
overlaps <- findOverlaps(mafa_gr, snp_gr)

# 3. Count SNPs per Peak
# We create a table of how many SNPs fall into each MafA site
snp_counts <- as.data.table(as.data.frame(overlaps))
snp_counts <- snp_counts[, .(SNP_Count = .N), by = queryHits]

# 4. Filter for SNP-containing Peaks
# Add the count back to the main mafa table
mafa[, SNP_Count := 0]
mafa[snp_counts$queryHits, SNP_Count := snp_counts$SNP_Count]

# Create the subset: Peaks with at least 1 SNP
mafa_with_snps <- mafa[SNP_Count > 0]

# Save the list for external review
fwrite(mafa_with_snps, file.path(BASE_DIR, "MafA_Peaks_with_SNPs.csv"))
message(paste("Identified", nrow(mafa_with_snps), "MafA sites containing at least one SNP."))

# 5. Histogram of Distances to TSS for SNP-Mafa Peaks
# We maintain the "forestgreen" and structural lines from the previous vignette
p_snp_tss <- ggplot(mafa_with_snps, aes(x = Dist_TSS_Num)) +
  geom_histogram(
    binwidth = 25000, 
    fill = "forestgreen", 
    color = "white"
  ) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.6) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.6) +
  coord_cartesian(xlim = c(-1000000, 1000000)) +
  scale_x_continuous(labels = function(x) paste0(x/1e6, " Mb")) +
  labs(
    title = "TSS Proximity of SNP-Containing MafA Sites",
    subtitle = paste("Analysis of", nrow(mafa_with_snps), "filtered sites | +/- 1 Mbp"),
    x = "Distance to TSS",
    y = "Frequency"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black")
  )

# 6. High-Resolution Inset for SNP-Mafa Peaks
p_snp_inset <- ggplot(mafa_with_snps, aes(x = Dist_TSS_Num)) +
  geom_histogram(binwidth = 25, fill = "forestgreen", color = "white") +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.5) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  coord_cartesian(xlim = c(-1000, 1000)) +
  scale_x_continuous(labels = function(x) paste0(x/1000, " kb")) +
  labs(title = NULL, x = NULL, y = NULL) +
  theme_minimal() +
  theme(
    plot.background = element_rect(fill = "white", color = "white", linewidth = 0.8),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text = element_text(size = 9, color = "black")
  )

# 7. Combine and Save
# Using the 30% height increase (bottom = 0.421)
final_snp_vignette <- p_snp_tss + 
  inset_element(p_snp_inset, 
                left = 0.6, 
                bottom = 0.25, 
                right = 0.98, 
                top = 0.98)

print(final_snp_vignette)

ggsave(
  filename = file.path(BASE_DIR, "Vignette_MafA_SNP_TSS_Intersection.png"),
  plot = final_snp_vignette,
  width = 10, height = 7, dpi = 300, bg = "white"
)

