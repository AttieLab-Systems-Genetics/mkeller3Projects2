# -------------------------------------------------------------------------
# Part 1: Environment Setup and Data Ingestion
# -------------------------------------------------------------------------

# Load necessary libraries
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("dplyr")) install.packages("dplyr")
if (!require("tidyverse")) install.packages("tidyverse")

library(ggplot2)
library(dplyr)
library(tidyverse)

# 1. Define the local path and file
file_path <- "C:/Users/mkeller3/Desktop/Kahn FAHFA analysis Jan2026"
file_name <- "DO500_adipose_fahfas_all_mice_additive_peaks.csv"

# 2. Load and Map Columns
# We rename the new columns to match the logic in the rest of the script
qtl_data <- read.csv(file.path(file_path, file_name), stringsAsFactors = FALSE) %>%
  rename(
    lodcolumn = phenotype,
    chr = qtl_chr,
    pos = qtl_pos,
    lod = qtl_lod
  )

# ==============================================================================
# MODULE 2: FAHFA QTL MANHATTAN (ADJUSTED SCALING)
# ==============================================================================

# GLOBAL GRCm39 CONFIGURATION
chr_lens <- c("1"=195.15, "2"=181.75, "3"=159.75, "4"=156.86, "5"=151.76, "6"=149.58, 
              "7"=144.99, "8"=130.13, "9"=124.36, "10"=130.53, "11"=121.97, "12"=120.09, 
              "13"=120.88, "14"=125.14, "15"=104.07, "16"=98.01, "17"=95.29, "18"=90.72, 
              "19"=61.42, "X"=169.48)

GLOBAL_MAP <- data.frame(Chr = names(chr_lens), Length = as.numeric(chr_lens)) %>%
  mutate(offset = lag(cumsum(Length), default = 0),
         center = offset + Length / 2,
         chr_col = rep(c("darkorange", "darkgreen"), length.out = 20))

TOTAL_LEN <- sum(GLOBAL_MAP$Length)

# PLOTTING FUNCTION
generate_fahfa_manhattan <- function(data_df, dataset_name, pt_size = 7) {
  
  df_plot <- data_df %>%
    mutate(Chr = as.character(chr)) %>%
    filter(Chr %in% names(chr_lens)) %>%
    inner_join(GLOBAL_MAP, by = "Chr") %>%
    mutate(BP_cum = offset + pos)
  
  p <- ggplot(df_plot, aes(x = BP_cum, y = lod)) +
    geom_hline(yintercept = 7, linetype = "dashed", color = "black", alpha = 0.5) +
    geom_point(aes(fill = chr_col), 
               shape = 21, 
               color = "black", 
               size = pt_size, 
               stroke = 0.4, 
               alpha = 0.85) +
    scale_fill_identity() +
    scale_x_continuous(breaks = GLOBAL_MAP$center, labels = GLOBAL_MAP$Chr, 
                       limits = c(0, TOTAL_LEN), expand = c(0.01, 0.01)) +
    scale_y_continuous(expand = c(0, 0), limits = c(6.5, max(df_plot$lod, na.rm = TRUE) + 2)) +
    theme_classic(base_size = 16) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(size = 9, face = "bold"),
      axis.title = element_text(face = "bold")
    ) +
    labs(
      title = "Manhattan plot",
      subtitle = "Additive QTL with LOD > 7",
      x = "Chromosome", 
      y = "LOD Score"
    )
  
  return(p)
}

# EXECUTION
p_fahfa <- generate_fahfa_manhattan(qtl_data, "FAHFA_Adipose", pt_size = 7)
print(p_fahfa)

# -------------------------------------------------------------------------
# Save the Manhattan Plot to Desktop
# -------------------------------------------------------------------------

# 1. Define the Desktop path
# Using 'Sys.getenv("USERPROFILE")' automatically finds your specific Windows user folder
desktop_path <- file.path(Sys.getenv("USERPROFILE"), "Desktop")

# 2. Save as a high-resolution PNG (300 DPI is standard for print)
ggsave(
  filename = file.path(desktop_path, "FAHFA_Adipose_Manhattan_HighRes.png"), 
  plot = p_fahfa, 
  width = 14,      # Width in inches
  height = 6,      # Height in inches
  dpi = 300        # Dots per inch (Resolution)
)

# 3. Save as a PDF (Vector format - best for publication)
ggsave(
  filename = file.path(desktop_path, "FAHFA_Adipose_Manhattan_Vector.pdf"), 
  plot = p_fahfa, 
  width = 14, 
  height = 6
)

message("Plots have been saved to your Desktop.")





# -------------------------------------------------------------------------
# Part 3: Allele Effect Heatmap (Selection Parameters)
# -------------------------------------------------------------------------

target_chr <- "X"     # Example: Set to 9 to see the major hotspot
start_mb   <- 52
end_mb     <- 57

generate_allele_heatmap <- function(data_df, chr_select, start_select, end_select) {
  
  plot_df <- data_df %>%
    filter(chr == chr_select, 
           pos >= start_select, 
           pos <= end_select,
           lod >= 7) %>% 
    arrange(pos) %>%
    select(lodcolumn, A, B, C, D, E, F, G, H) %>%
    pivot_longer(cols = A:H, names_to = "Strain", values_to = "Effect")
  
  plot_df$Strain <- factor(plot_df$Strain, levels = c("A", "B", "C", "D", "E", "F", "G", "H"))
  plot_df$lodcolumn <- factor(plot_df$lodcolumn, levels = unique(plot_df$lodcolumn))
  
  p <- ggplot(plot_df, aes(x = Strain, y = lodcolumn, fill = Effect)) +
    geom_tile(color = "white", linewidth = 0.3) +
    scale_fill_gradient2(low = "royalblue4", 
                         mid = "white", 
                         high = "firebrick3", 
                         midpoint = 0, 
                         limits = c(-0.7, 0.7), 
                         oob = scales::squish, 
                         name = "Allele effect") +
    theme_minimal() +
    labs(title = paste("Chr", chr_select, ":", start_select, "-", end_select, "Mb"),
         #x = "Founder Strain",
         y = "FAHFA lipid") +
    theme(axis.text.y = element_text(size = 10),
          axis.text.x = element_text(face = "bold", size = 11),
          panel.grid = element_blank(),
          plot.title = element_text(face = "bold"))
  
  return(p)
}

p_heatmap <- generate_allele_heatmap(qtl_data, target_chr, start_mb, end_mb)
print(p_heatmap)




# ==============================================================================
# MODULE 8: SEX-MIRROR MANHATTAN (TRUE 6.8 BASELINE)
# ==============================================================================

if (!require("patchwork")) install.packages("patchwork")
library(tidyverse)
library(patchwork)

# 1. SET PATHS & LOAD DATA
desktop_path <- file.path(Sys.getenv("USERPROFILE"), "Desktop")

# Load and rename columns to standardize
qtl_female <- read.csv(file.path(desktop_path, "DO500_adipose_fahfas_female_mice_additive_peaks.csv")) %>%
  rename(lodcolumn = phenotype, chr = qtl_chr, pos = qtl_pos, lod = qtl_lod)

qtl_male <- read.csv(file.path(desktop_path, "DO500_adipose_fahfas_male_mice_additive_peaks.csv")) %>%
  rename(lodcolumn = phenotype, chr = qtl_chr, pos = qtl_pos, lod = qtl_lod)

# 2. CALCULATE GLOBAL SCALE LIMITS (Ensures both panels are identical)
# We find the highest LOD score across both datasets to set the outer limit
max_lod_val <- max(c(qtl_female$lod, qtl_male$lod), na.rm = TRUE)
y_max_limit <- ceiling(max_lod_val) + 0.5
y_min_limit <- 6.0  # Your requested baseline

# 3. GLOBAL GRCm39 MAP CONFIGURATION
chr_lens <- c("1"=195.15, "2"=181.75, "3"=159.75, "4"=156.86, "5"=151.76, "6"=149.58, 
              "7"=144.99, "8"=130.13, "9"=124.36, "10"=130.53, "11"=121.97, "12"=120.09, 
              "13"=120.88, "14"=125.14, "15"=104.07, "16"=98.01, "17"=95.29, "18"=90.72, 
              "19"=61.42, "X"=169.48)

GLOBAL_MAP <- data.frame(Chr = names(chr_lens), Length = as.numeric(chr_lens)) %>%
  mutate(offset = lag(cumsum(Length), default = 0),
         center = offset + Length / 2,
         chr_col = rep(c("darkorange", "darkgreen"), length.out = 20))

TOTAL_LEN <- sum(GLOBAL_MAP$Length)

# Helper function to map genomic positions
prep_mirror_data <- function(df) {
  df %>%
    mutate(Chr = as.character(chr)) %>%
    filter(Chr %in% names(chr_lens)) %>%
    inner_join(GLOBAL_MAP, by = "Chr") %>%
    mutate(BP_cum = offset + pos)
}

df_female_plot <- prep_mirror_data(qtl_female)
df_male_plot   <- prep_mirror_data(qtl_male)

# 4. CONSTRUCT THE TOP PANEL (FEMALE)
p_female <- ggplot(df_female_plot, aes(x = BP_cum, y = lod)) +
  geom_hline(yintercept = 7, linetype = "dashed", color = "red", alpha = 0.5) +
  geom_point(aes(fill = chr_col), shape = 21, color = "black", size = 6, stroke = 0.4, alpha = 0.8) +
  scale_fill_identity() +
  scale_x_continuous(limits = c(0, TOTAL_LEN), expand = c(0.01, 0.01)) +
  # Standard axis: 6.8 (bottom) to Max (top)
  scale_y_continuous(limits = c(y_min_limit, y_max_limit), expand = c(0, 0)) +
  theme_classic(base_size = 16) +
  labs(title = "Sex-Mirror Manhattan: Adipose FAHFAs", 
       subtitle = "Mirror point at LOD 6.0",
       y = "Female LOD") +
  theme(
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin  = margin(b = 0) # Remove bottom margin to touch the lower plot
  )

# 5. CONSTRUCT THE BOTTOM PANEL (MALE)
p_male <- ggplot(df_male_plot, aes(x = BP_cum, y = lod)) +
  geom_hline(yintercept = 7, linetype = "dashed", color = "red", alpha = 0.5) +
  geom_point(aes(fill = chr_col), shape = 21, color = "black", size = 6, stroke = 0.4, alpha = 0.8) +
  scale_fill_identity() +
  scale_x_continuous(breaks = GLOBAL_MAP$center, labels = GLOBAL_MAP$Chr, 
                     limits = c(0, TOTAL_LEN), expand = c(0.01, 0.01)) +
  # REVERSED axis: 6.8 (top) to Max (bottom)
  scale_y_reverse(limits = c(y_max_limit, y_min_limit), expand = c(0, 0)) +
  theme_classic(base_size = 16) +
  labs(x = "Chromosome (GRCm39)", y = "Male LOD") +
  theme(
    plot.margin = margin(t = 0), # Remove top margin to touch the upper plot
    axis.text.x = element_text(size = 12, face = "bold")
  )

# 6. COMBINE AND SAVE
# The '/' operator from patchwork stacks plots vertically
p_mirror_final <- p_female / p_male

print(p_mirror_final)

ggsave(file.path(desktop_path, "FAHFA_Mirror_Manhattan_6.8_Baseline.png"), 
       p_mirror_final, width = 12.5, height = 6, dpi = 300)








# ==============================================================================
# MODULE 2B: SEX DIFFERENCE QTL MANHATTAN (GRCm39 SANITIZED)
# ==============================================================================

library(tidyverse)

# 1. LOAD THE SEX DIFFERENCE DATA
file_path <- "C:/Users/mkeller3/Desktop/Kahn FAHFA analysis Jan2026"
file_diff <- "Sex difference QTL for FAHFA lipids.csv"

qtl_diff_data <- read.csv(file.path(file_path, file_diff), stringsAsFactors = FALSE)

# 2. GLOBAL GRCm39 CONFIGURATION (Ensures all Chrs 1-X are shown)
chr_lens <- c("1"=195.15, "2"=181.75, "3"=159.75, "4"=156.86, "5"=151.76, "6"=149.58, 
              "7"=144.99, "8"=130.13, "9"=124.36, "10"=130.53, "11"=121.97, "12"=120.09, 
              "13"=120.88, "14"=125.14, "15"=104.07, "16"=98.01, "17"=95.29, "18"=90.72, 
              "19"=61.42, "X"=169.48)

GLOBAL_MAP <- data.frame(Chr = names(chr_lens), Length = as.numeric(chr_lens)) %>%
  mutate(offset = lag(cumsum(Length), default = 0),
         center = offset + Length / 2,
         chr_col = rep(c("royalblue4", "skyblue3"), length.out = 20)) # Specific blue palette

TOTAL_LEN <- sum(GLOBAL_MAP$Length)

# 3. PLOTTING FUNCTION
generate_diff_manhattan <- function(data_df, dataset_name, pt_size = 5.0) {
  
  # A. Sanitize and Align Coordinates
  df_plot <- data_df %>%
    mutate(Chr = as.character(chr)) %>%
    filter(Chr %in% names(chr_lens)) %>%
    inner_join(GLOBAL_MAP, by = "Chr") %>%
    mutate(BP_cum = offset + pos)
  
  # B. The Plot
  p <- ggplot(df_plot, aes(x = BP_cum, y = lod)) +
    # Significance Line
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", alpha = 0.5) +
    
    # --- OPTIMIZATION: CIRCLE STYLE ---
    geom_point(aes(fill = chr_col), 
               shape = 21,          # Circle with border
               color = "black",     # Border color
               size = pt_size,      # <--- CHANGE POINT SIZE HERE
               stroke = 0.3,        # Border thickness
               alpha = 0.8) + 
    # ----------------------------------
  
  scale_fill_identity() +
    scale_x_continuous(breaks = GLOBAL_MAP$center, labels = GLOBAL_MAP$Chr, 
                       limits = c(0, TOTAL_LEN), expand = c(0.01, 0.01)) +
    
    # Y-axis starts at 7 (filters out low-LOD noise visually)
    scale_y_continuous(expand = c(0, 0), limits = c(5.5, max(df_plot$lod) + 1)) +
    
    theme_classic(base_size = 14) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(size = 9, angle = 0, face = "bold"),
      axis.title = element_text(face = "bold")
    ) +
    labs(
      title = paste("Manhattan Plot:", dataset_name),
      subtitle = "Interactive Scan (Difference Profiles) - Threshold LOD > 7",
      x = "Chromosome (GRCm39)", 
      y = "LOD Score (Sex Interaction)"
    )
  
  return(p)
}

# 4. EXECUTION
# You can adjust pt_size here to optimize the visual density
p_diff <- generate_diff_manhattan(qtl_diff_data, "Sex Difference QTL", pt_size = 4.5)

# Display
print(p_diff)

# 5. SAVE THE PLOT
ggsave(file.path(file_path, "Sex_Difference_Manhattan_GRCm39.png"), 
       p_diff, width = 12, height = 5, dpi = 300)





# ==============================================================================
# MODULE: DUAL MANHATTAN (ADDITIVE TOP | SEX DIFFERENCE BOTTOM)
# ==============================================================================

library(tidyverse)
if (!require("patchwork")) install.packages("patchwork")
library(patchwork)

# 1. SETUP PATHS AND DATA
file_path <- "C:/Users/mkeller3/Desktop/Kahn FAHFA analysis Jan2026"
file_add  <- "Additive_QTL_alleles_FAHFA_DO_adipose.csv"
file_diff <- "Sex difference QTL for FAHFA lipids.csv"

qtl_add  <- read.csv(file.path(file_path, file_add), stringsAsFactors = FALSE)
qtl_diff <- read.csv(file.path(file_path, file_diff), stringsAsFactors = FALSE)

# 2. GLOBAL GRCm39 CONFIGURATION
chr_lens <- c("1"=195.15, "2"=181.75, "3"=159.75, "4"=156.86, "5"=151.76, "6"=149.58, 
              "7"=144.99, "8"=130.13, "9"=124.36, "10"=130.53, "11"=121.97, "12"=120.09, 
              "13"=120.88, "14"=125.14, "15"=104.07, "16"=98.01, "17"=95.29, "18"=90.72, 
              "19"=61.42, "X"=169.48)

GLOBAL_MAP <- data.frame(Chr = names(chr_lens), Length = as.numeric(chr_lens)) %>%
  mutate(offset = lag(cumsum(Length), default = 0),
         center = offset + Length / 2,
         chr_col = rep(c("darkorange", "darkgreen"), length.out = 20))

TOTAL_LEN <- sum(GLOBAL_MAP$Length)

# ==============================================================================
# 3. TOP PANEL: ADDITIVE QTL
# ==============================================================================
# --- MANUALLY OPTIMIZE TOP PANEL HERE ---
y_min_top <- 6.8
y_max_top <- 15 # Adjust based on your highest LOD
pt_size_top <- 5
# ----------------------------------------

df_top <- qtl_add %>%
  mutate(Chr = as.character(chr)) %>%
  filter(Chr %in% names(chr_lens)) %>%
  inner_join(GLOBAL_MAP, by = "Chr") %>%
  mutate(BP_cum = offset + pos)

p_top <- ggplot(df_top, aes(x = BP_cum, y = lod)) +
  geom_hline(yintercept = y_min_top, linetype = "dashed", color = "red", alpha = 0.5) +
  geom_point(aes(fill = chr_col), shape = 21, color = "black", size = pt_size_top, stroke = 0.3, alpha = 0.8) +
  scale_fill_identity() +
  scale_x_continuous(breaks = GLOBAL_MAP$center, labels = GLOBAL_MAP$Chr, limits = c(0, TOTAL_LEN), expand = c(0.01, 0.01)) +
  scale_y_continuous(limits = c(y_min_top, y_max_top), expand = c(0, 0)) +
  theme_classic() +
  labs(title = "Additive QTL Manhattan, LOD>7", x = NULL, y = "Additive LOD") +
  theme(panel.grid = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())

# ==============================================================================
# 4. BOTTOM PANEL: SEX DIFFERENCE QTL (Interaction - Additive)
# ==============================================================================
# --- MANUALLY OPTIMIZE BOTTOM PANEL HERE ---
y_min_bot <- 6
y_max_bot <- 10 # Adjust based on interaction peaks
pt_size_bot <- 5
# -------------------------------------------

df_bot <- qtl_diff %>%
  mutate(Chr = as.character(chr)) %>%
  filter(Chr %in% names(chr_lens)) %>%
  inner_join(GLOBAL_MAP, by = "Chr") %>%
  mutate(BP_cum = offset + pos)

p_bot <- ggplot(df_bot, aes(x = BP_cum, y = lod)) +
  # Add the dotted line at y=7
  geom_hline(yintercept = 7, linetype = "dashed", color = "black", linewidth = 0.8) +
  
  geom_point(aes(fill = chr_col), shape = 21, color = "black", size = 5, stroke = 0.3, alpha = 0.8) +
  scale_fill_identity() +
  scale_x_continuous(breaks = GLOBAL_MAP$center, labels = GLOBAL_MAP$Chr, limits = c(0, TOTAL_LEN), expand = c(0.01, 0.01)) +
  scale_y_continuous(limits = c(y_min_bot, y_max_bot), expand = c(0, 0)) +
  theme_classic() +
  labs(title = "Sex-interactive QTL Manhattan, LOD Diff >6", 
       x = "Chromosome (GRCm39)", 
       y = "Sex Diff LOD") +
  theme(panel.grid = element_blank(), 
        axis.text.x = element_text(size = 9, face = "bold"))

# ==============================================================================
# 5. COMBINE AND SAVE
# ==============================================================================
# Using patchwork to stack (/) and align axes
p_combined <- p_top / p_bot

print(p_combined)

ggsave(file.path(file_path, "FAHFA_Dual_Manhattan_Additive_vs_SexDiff.png"), 
       p_combined, width = 12, height = 6, dpi = 300)




# ==============================================================================
# MODULE 7: QTL FOOTPRINT (WITH LOD THRESHOLDING)
# ==============================================================================

library(tidyverse)
library(viridis)
library(tidyr)

# --- USER OPTIMIZATION PANEL ---
# SELECT DATASET: qtl_add or qtl_diff
data_to_plot <- qtl_diff  

# FILTERING
lod_threshold <- 7.0             # <--- SET YOUR MINIMUM LOD HERE

# AESTHETICS
genome_bg_color <- "grey95"      
chr_line_color  <- "grey50"      
chr_line_width  <- 0.05          
y_label_size    <- 8             
point_width     <- 12            
# -------------------------------

# 1. Prepare Data & Apply LOD Threshold
footprint_data <- data_to_plot %>%
  filter(lod >= lod_threshold) %>% # <--- APPLY FILTER HERE
  mutate(Chr = as.character(chr)) %>%
  filter(Chr %in% names(chr_lens)) %>%
  inner_join(GLOBAL_MAP, by = "Chr") %>%
  mutate(BP_cum = offset + pos)

# Check if data exists after filtering
if(nrow(footprint_data) == 0) stop("No QTLs found above the current LOD threshold.")

# Set Title and Filename
plot_title <- if(identical(data_to_plot, qtl_add)) "Additive QTL Heatmap" else "Sex Difference QTL Heatmap"
save_name  <- if(identical(data_to_plot, qtl_add)) "FAHFA_Additive_Filtered.png" else "FAHFA_SexDiff_Filtered.png"

# 2. Hierarchical Clustering (Updated based on filtered hits)
bin_size <- 5 
cluster_matrix <- footprint_data %>%
  mutate(bin = floor(BP_cum / bin_size)) %>%
  group_by(lodcolumn, bin) %>%
  summarize(max_lod = max(lod), .groups = 'drop') %>%
  pivot_wider(names_from = bin, values_from = max_lod, values_fill = 0)

mat <- as.matrix(cluster_matrix[,-1])
rownames(mat) <- cluster_matrix$lodcolumn
hc <- hclust(dist(mat), method = "complete")
lipid_clustered_order <- rownames(mat)[hc$order]
footprint_data$lodcolumn <- factor(footprint_data$lodcolumn, levels = lipid_clustered_order)

# 3. Generate the Heatmap
p_footprint_unified <- ggplot(footprint_data, aes(x = BP_cum, y = lodcolumn)) +
  geom_tile(aes(fill = lod), width = point_width, height = 0.85) + 
  scale_fill_viridis_c(option = "magma", direction = -1, name = "LOD Score", end = 0.9) +
  geom_vline(xintercept = GLOBAL_MAP$offset, color = chr_line_color, linewidth = chr_line_width) +
  scale_x_continuous(breaks = GLOBAL_MAP$center, labels = GLOBAL_MAP$Chr, 
                     limits = c(0, TOTAL_LEN), expand = c(0, 0)) +
  theme_minimal() +
  labs(title = paste("Genome-Wide", plot_title),
       subtitle = paste("Filtered for LOD >=", lod_threshold),
       x = "Chromosome (GRCm39)",
       y = "FAHFA Lipid Phenotype") +
  theme(
    panel.background = element_rect(fill = genome_bg_color, color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.text.y      = element_text(size = y_label_size, family = "mono"), 
    axis.text.x      = element_text(face = "bold", size = 10),
    panel.grid       = element_blank(),
    legend.position  = "right",
    plot.title       = element_text(face = "bold", size = 16)
  )

# 4. Save and Display
print(p_footprint_unified)


ggsave(file.path(file_path, save_name), p_footprint_unified, width = 18, height = 14, dpi = 300)


