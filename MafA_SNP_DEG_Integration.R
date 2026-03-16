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
  labs(title = NULL, x = "Distance to TSS", y = "Frequency") +
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
  width = 12, 
  height = 6, 
  dpi = 600, 
  bg = "white"
)




###----------------------------------------------------
# Vignette: SNP-Containing MafA Binding Sites
###----------------------------------------------------

# Define the global limit once for all plots
TOTAL_GENOME_LENGTH <- max(CHR_INFO$EndOffset)

# 2. Re-filter for Proximal Sites (1 or more DEGs within 200kb)
PROXIMAL_WINDOW <- 200000

# Create GRanges
mafa_gr <- GRanges(seqnames = mafa$Chr, 
                   ranges = IRanges(start = mafa$mm39_start, end = mafa$mm39_end))
deg_gr  <- GRanges(seqnames = degs$Chr, 
                   ranges = IRanges(start = degs$LocalPos, width = 1))

# Find indices of MafA sites near at least one DEG
prox_overlaps    <- findOverlaps(mafa_gr, deg_gr, maxgap = PROXIMAL_WINDOW)
proximal_indices <- unique(queryHits(prox_overlaps))

# Create the two datasets
mafa_all_snps <- mafa[SNP_Count > 0]
mafa_prox_snps <- mafa[proximal_indices][SNP_Count > 0]

# Ensure GlobalPos is present in both (using the Offset from CHR_INFO)
if(!"Offset" %in% names(mafa_all_snps)) {
  mafa_all_snps <- merge(mafa_all_snps, CHR_INFO[, .(Chr, Offset)], by = "Chr")
  mafa_all_snps[, GlobalPos := mm39_start + Offset]
}
if(!"Offset" %in% names(mafa_prox_snps)) {
  mafa_prox_snps <- merge(mafa_prox_snps, CHR_INFO[, .(Chr, Offset)], by = "Chr")
  mafa_prox_snps[, GlobalPos := mm39_start + Offset]
}

##--------------------------------------------------------------------
# All mafa peaks
##--------------------------------------------------------------------

p_all <- ggplot() +
  # Background stripes
  geom_rect(data = CHR_INFO, aes(xmin = Offset, xmax = EndOffset, ymin = -Inf, ymax = Inf, fill = Shade), 
            alpha = 1, show.legend = FALSE) +
  scale_fill_identity() +
  
  # Data Points
  geom_point(data = mafa_all_snps, 
             aes(x = GlobalPos, y = `Peak Score`, size = SNP_Count), 
             shape = 21, fill = "forestgreen", color = "black", stroke = 0.2, alpha = 0.6) +        
  
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  
  # FIXED X-AXIS
  scale_x_continuous(
    limits = c(0, TOTAL_GENOME_LENGTH),
    breaks = CHR_INFO$MidpointGlobal, 
    labels = CHR_INFO$AxisLabel,
    expand = c(0, 0)
  ) +
  scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.1))) +
  scale_size_continuous(range = c(0.8, 8), name = "SNPs per Peak") +
  
  labs(title = "SNP-Impacted MafA Binding Sites (All)", x = "Chromosome", y = "MafA Peak Score") +
  theme_minimal() + 
  theme(panel.grid = element_blank(), axis.text.x = element_text(size = 10, face = "bold"),
        axis.line.y = element_line(color = "black"), legend.position = "bottom")

print(p_all)

ggsave(file.path(BASE_DIR, "Manhattan_MafA_SNP_ALL_Fixed.png"), 
       plot = p_all, width = 13, height = 6.5, dpi = 600, bg = "white")

##--------------------------------------------------------------------
# Mafa peaks proximal to DEGs
##--------------------------------------------------------------------


p_prox <- ggplot() +
  # Identical background stripes
  geom_rect(data = CHR_INFO, aes(xmin = Offset, xmax = EndOffset, ymin = -Inf, ymax = Inf, fill = Shade), 
            alpha = 1, show.legend = FALSE) +
  scale_fill_identity() +
  
  # Filtered Data Points
  geom_point(data = mafa_prox_snps, 
             aes(x = GlobalPos, y = `Peak Score`, size = SNP_Count), 
             shape = 21, fill = "forestgreen", color = "black", stroke = 0.2, alpha = 0.6) +        
  
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  
  # IDENTICAL FIXED X-AXIS
  scale_x_continuous(
    limits = c(0, TOTAL_GENOME_LENGTH),
    breaks = CHR_INFO$MidpointGlobal, 
    labels = CHR_INFO$AxisLabel,
    expand = c(0, 0)
  ) +
  scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.1))) +
  scale_size_continuous(range = c(0.8, 8), name = "SNPs per Peak") +
  
  labs(title = "SNP-Impacted MafA Binding Sites Proximal to DEGs (+/- 200kb)", 
       x = "Chromosome", y = "MafA Peak Score") +
  theme_minimal() + 
  theme(panel.grid = element_blank(), axis.text.x = element_text(size = 10, face = "bold"),
        axis.line.y = element_line(color = "black"), legend.position = "bottom")

print(p_prox)

ggsave(file.path(BASE_DIR, "Manhattan_MafA_SNP_PROXIMAL_Fixed.png"), 
       plot = p_prox, width = 13, height = 6.5, dpi = 600, bg = "white")





###----------------------------------------------------
# Vignette: Functional Manhattan (Pixel-Perfect Alignment)
###----------------------------------------------------

library(ggtext)
library(ggnewscale)

# 1. Data Mapping
nearest_deg_idx <- nearest(mafa_prox_gr, deg_gr)
mafa_prox_snps[, Strain := degs$Strain[nearest_deg_idx]]
mafa_prox_snps[, FC := degs$FC[nearest_deg_idx]]
mafa_prox_snps[, Direction := ifelse(FC > 0, "UP", "DOWN")]

# 2. Setup Shared Aesthetics
impact_shapes <- c("UP" = 24, "DOWN" = 25) 
strain_colors <- c("Mixed" = "#CC0000", "C57BL/6J" = "#0000CC")

# 3. Build the Plot
p_directional_final <- ggplot() +
  # INFRASTRUCTURE
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), fill = "white") +
  geom_rect(data = CHR_INFO, 
            aes(xmin = Offset, xmax = EndOffset, ymin = -Inf, ymax = Inf, fill = Shade), 
            alpha = 1, show.legend = FALSE) +
  scale_fill_identity() + 
  
  # DATA (Directional Triangles)
  new_scale_fill() + 
  geom_point(data = mafa_prox_snps, 
             aes(x = GlobalPos, y = `Peak Score`, 
                 fill = Strain, 
                 shape = Direction, 
                 size = SNP_Count), 
             color = "black", 
             stroke = 0.25, 
             alpha = 0.6) +       
  
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  
  # SCALES
  scale_fill_manual(values = strain_colors, guide = "none") + 
  scale_shape_manual(values = impact_shapes, name = "DEG Direction") +
  scale_size_continuous(range = c(0.8, 8), name = "SNPs per Peak") +
  
  # AXES
  scale_x_continuous(limits = c(0, TOTAL_GENOME_LENGTH),
                     breaks = CHR_INFO$MidpointGlobal, 
                     labels = CHR_INFO$AxisLabel,
                     expand = c(0, 0)) +
  scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.1))) +
  
  # TITLES (HTML color-coded)
  labs(
    title = paste0("Functional Impact of MafA Sites for ",
                   "<span style='color:#0000CC;'>C57BL/6J</span> and ",
                   "<span style='color:#CC0000;'>SJL/Mixed</span> Strains"),
    x = "Chromosome", 
    y = "MafA Peak Score"
  ) +
  
  # THEME (Force alignment by zeroing margins)
  theme_minimal() + 
  theme(
    panel.grid = element_blank(), 
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.text.x = element_text(size = 10, face = "bold"),
    axis.line.y = element_line(color = "black"), 
    axis.ticks.y = element_blank(),
    legend.position = "bottom",
    legend.key = element_blank(),
    legend.box = "horizontal",
    
    # ALIGNMENT FIX: Explicitly zero out subtitle and set title margin
    plot.title = element_markdown(lineheight = 1.1, margin = margin(t = 5, b = 5)),
    plot.subtitle = element_blank(), # Removes the subtitle object entirely
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10) # Standardizes outer spacing
  ) +
  
  # LEGEND GUIDES
  guides(
    shape = guide_legend(override.aes = list(size = 4, fill = "grey50")),
    size = guide_legend(override.aes = list(shape = 21, fill = "grey80"))
  )

print(p_directional_final)

# 4. Save
ggsave(
  filename = file.path(BASE_DIR, "Manhattan_MafA_Directional_Perfect_Align.png"), 
  plot = p_directional_final,
  width = 13, height = 6.5, dpi = 600, bg = "white"
)


####---------------------------------------------------------
# Harmonize plot layouts for mafa peaks, mafa peaks with SNPs and 
# mafa peaks with SNPs proximal to directional DEGs via triangles
####---------------------------------------------------------




# 1. Determine the absolute max Peak Score across all datasets to lock the Y-axis
global_max_y <- max(mafa_all_snps$`Peak Score`, na.rm = TRUE) * 1.1 

# 2. Create a standardized theme object
shared_manhattan_theme <- theme_minimal() + 
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 10, face = "bold"),
    axis.line.y = element_line(color = "black"),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.key = element_blank(),
    plot.title = ggtext::element_markdown(lineheight = 1.1, margin = margin(t = 5, b = 20)),
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5),
    # Force the plot area to be exactly the same regardless of legend complexity
    aspect.ratio = 6.5 / 13 
  )

# 3. Create a shared Y-scale to prevent "wobble"
shared_y_scale <- scale_y_continuous(
  limits = c(0, global_max_y), 
  labels = scales::comma, 
  expand = c(0, 0)
)


p_all <- ggplot() +
  geom_rect(data = CHR_INFO, aes(xmin = Offset, xmax = EndOffset, ymin = -Inf, ymax = Inf, fill = Shade), alpha = 1, show.legend = FALSE) +
  scale_fill_identity() +
  geom_point(data = mafa_all_snps, aes(x = GlobalPos, y = `Peak Score`, size = SNP_Count), 
             shape = 21, fill = "forestgreen", color = "black", stroke = 0.2, alpha = 0.6) +
  scale_x_continuous(limits = c(0, TOTAL_GENOME_LENGTH), breaks = CHR_INFO$MidpointGlobal, labels = CHR_INFO$AxisLabel, expand = c(0, 0)) +
  shared_y_scale + 
  scale_size_continuous(range = c(0.8, 8), name = "SNPs per Peak", limits = c(1, max(mafa_all_snps$SNP_Count))) +
  labs(title = "SNP-Impacted MafA Binding Sites (All)", x = "Chromosome", y = "MafA Peak Score") +
  shared_manhattan_theme

ggsave(file.path(BASE_DIR, "1_All_SNPs.png"), 
       plot = p_all, 
       width = 13, 
       height = 6.5, 
       dpi = 600)


p_prox <- ggplot() +
  geom_rect(data = CHR_INFO, aes(xmin = Offset, xmax = EndOffset, ymin = -Inf, ymax = Inf, fill = Shade), alpha = 1, show.legend = FALSE) +
  scale_fill_identity() +
  geom_point(data = mafa_prox_snps, aes(x = GlobalPos, y = `Peak Score`, size = SNP_Count), 
             shape = 21, fill = "forestgreen", color = "black", stroke = 0.2, alpha = 0.6) +
  scale_x_continuous(limits = c(0, TOTAL_GENOME_LENGTH), breaks = CHR_INFO$MidpointGlobal, labels = CHR_INFO$AxisLabel, expand = c(0, 0)) +
  shared_y_scale +
  scale_size_continuous(range = c(0.8, 8), name = "SNPs per Peak", limits = c(1, max(mafa_all_snps$SNP_Count))) +
  labs(title = "SNP-Impacted MafA Binding Sites Proximal to DEGs", x = "Chromosome", y = "MafA Peak Score") +
  shared_manhattan_theme

ggsave(file.path(BASE_DIR, "2_Proximal_SNPs.png"), 
       plot = p_prox, 
       width = 13, 
       height = 6.5, 
       dpi = 600)


p_directional_final <- ggplot() +
  geom_rect(data = CHR_INFO, aes(xmin = Offset, xmax = EndOffset, ymin = -Inf, ymax = Inf, fill = Shade), alpha = 1, show.legend = FALSE) +
  scale_fill_identity() +
  new_scale_fill() +
  geom_point(data = mafa_prox_snps, aes(x = GlobalPos, y = `Peak Score`, fill = Strain, shape = Direction, size = SNP_Count), 
             color = "black", stroke = 0.25, alpha = 0.6) +
  scale_fill_manual(values = strain_colors, name = "Strain") + 
  scale_shape_manual(values = impact_shapes, name = "DEG Direction") +
  scale_x_continuous(limits = c(0, TOTAL_GENOME_LENGTH), breaks = CHR_INFO$MidpointGlobal, labels = CHR_INFO$AxisLabel, expand = c(0, 0)) +
  shared_y_scale +
  scale_size_continuous(range = c(0.8, 8), name = "SNPs per Peak", limits = c(1, max(mafa_all_snps$SNP_Count))) +
  labs(title = "Functional Impact of MafA Sites for <span style='color:#0000CC;'>C57BL/6J</span> and <span style='color:#CC0000;'>SJL/Mixed</span>", 
       x = "Chromosome", y = "MafA Peak Score") +
  shared_manhattan_theme +
  guides(shape = guide_legend(override.aes = list(size = 4, fill = "grey50")),
         fill = guide_legend(override.aes = list(shape = 21)),
         size = guide_legend(override.aes = list(shape = 21, fill = "grey80")))

ggsave(file.path(BASE_DIR, "3_Functional_Impact.png"), 
       plot = p_directional_final, 
       width = 13, 
       height = 6.5, 
       dpi = 600)



# 1. Total mafa binding sites with SNPs
total_snp_sites <- nrow(mafa[SNP_Count > 0])

# 2. Subset of those with at least one DEG within 200kb
# This matches your mafa_prox_snps logic
proximal_snp_sites <- nrow(mafa_prox_snps)

# Print results
cat("Total MafA sites with SNPs:", total_snp_sites, "\n")
cat("SNP sites near at least one DEG (200kb):", proximal_snp_sites, "\n")














