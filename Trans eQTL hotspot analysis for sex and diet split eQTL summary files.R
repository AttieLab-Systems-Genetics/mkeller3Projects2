# Load necessary libraries
library(tidyverse)
library(data.table)

# ==========================================
# 1. Configuration and Directory Setup
# ==========================================
input_dir <- "W:/General/main_directory/annotated_peak_summaries"
base_output_dir <- "C:/Users/mkeller3/Desktop/eQTL_Analysis"
task_name <- "Hotspots"
output_dir <- file.path(base_output_dir, task_name)

# Create directories if they don't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Define filenames based on your provided list
files <- list(
  HF     = "DO1200_liver_genes_HF_mice_additive_peaks.csv",
  HC     = "DO1200_liver_genes_HC_mice_additive_peaks.csv",
  Male   = "DO1200_liver_genes_male_mice_additive_peaks.csv",
  Female = "DO1200_liver_genes_female_mice_additive_peaks.csv"
)

# ==========================================
# 2. Module: Data Loading and Filtering
# ==========================================
load_trans_eqtls <- function(dir, filename, group_label) {
  path <- file.path(dir, filename)
  if (!file.exists(path)) {
    warning(paste("File not found:", path))
    return(NULL)
  }
  
  df <- fread(path) %>%
    as_tibble() %>%
    filter(cis == FALSE) %>%
    # Ensure positions are in Mb if they look like bp
    mutate(qtl_pos = if_else(max(qtl_pos, na.rm = TRUE) > 10000, qtl_pos / 1e6, qtl_pos)) %>%
    select(qtl_chr, qtl_pos, gene_symbol) %>%
    mutate(Group = group_label)
  
  return(df)
}

# ==========================================
# 3. Module: Sliding Window Calculation
# ==========================================
# This function counts eQTLs in a 1Mb window sliding by 0.5Mb increments
calculate_hotspot_density <- function(data, window_size_mb = 1, step_mb = 0.5) {
  results <- list()
  
  # Ensure chromosomes are handled as factors or specific order
  chrs <- unique(data$qtl_chr)
  
  for (chr in chrs) {
    chr_data <- data %>% filter(qtl_chr == chr)
    if (nrow(chr_data) == 0) next
    
    max_pos <- max(chr_data$qtl_pos, na.rm = TRUE)
    # Create sliding window centers
    windows <- seq(0, max_pos, by = step_mb)
    
    counts <- sapply(windows, function(w_center) {
      sum(chr_data$qtl_pos >= (w_center - window_size_mb/2) & 
            chr_data$qtl_pos <  (w_center + window_size_mb/2))
    })
    
    results[[as.character(chr)]] <- tibble(
      qtl_chr = chr,
      pos_mb = windows,
      count = counts
    )
  }
  
  return(bind_rows(results))
}

# ==========================================
# 4. Module: Plotting
# ==========================================
plot_hotspots <- function(density_data, title, filename) {
  # Clean up chromosome levels for plotting
  density_data$qtl_chr <- factor(density_data$qtl_chr, 
                                 levels = c(1:19, "X", "Y", "M"))
  
  p <- ggplot(density_data, aes(x = pos_mb, y = count, color = Group)) +
    geom_line(size = 0.8) +
    facet_wrap(~qtl_chr, scales = "free_x", nrow = 2) +
    theme_minimal() +
    labs(title = title,
         x = "Genomic Position (Mb)",
         y = "Trans-eQTL Count (1Mb Window)",
         color = "Group") +
    theme(legend.position = "bottom",
          strip.background = element_rect(fill = "grey90"))
  
  ggsave(file.path(output_dir, filename), plot = p, width = 14, height = 7)
  return(p)
}

# ==========================================
# 5. Execution Logic
# ==========================================

# --- Process DIET ---
cat("Processing Diet Hotspots...\n")
diet_data <- bind_rows(
  load_trans_eqtls(input_dir, files$HF, "High Fat"),
  load_trans_eqtls(input_dir, files$HC, "High Carbohydrate")
)

diet_density <- diet_data %>%
  group_by(Group) %>%
  do(calculate_hotspot_density(.))

plot_hotspots(diet_density, "Trans-eQTL Hotspots by Diet", "Hotspots_by_Diet.png")

# --- Process SEX ---
cat("Processing Sex Hotspots...\n")
sex_data <- bind_rows(
  load_trans_eqtls(input_dir, files$Male, "Male"),
  load_trans_eqtls(input_dir, files$Female, "Female")
)

sex_density <- sex_data %>%
  group_by(Group) %>%
  do(calculate_hotspot_density(.))

plot_hotspots(sex_density, "Trans-eQTL Hotspots by Sex", "Hotspots_by_Sex.png")

cat(paste("Analysis complete. Files saved to:", output_dir))


# ==========================================
# 6. Module: Differential Hotspot Identification
# ==========================================

find_differential_hotspots <- function(density_data, raw_data, threshold = 10) {
  # 1. Pivot data to compare groups side-by-side
  comparison <- density_data %>%
    pivot_wider(names_from = Group, values_from = count, values_fill = 0)
  
  # Get group names dynamically
  groups <- unique(density_data$Group)
  col1 <- groups[1]
  col2 <- groups[2]
  
  # 2. Calculate difference
  diff_hotspots <- comparison %>%
    mutate(diff = abs(.data[[col1]] - .data[[col2]])) %>%
    filter(diff >= threshold) %>%
    arrange(desc(diff))
  
  # 3. Extract traits mapping to these windows
  mapping_traits <- list()
  
  for(i in seq_len(nrow(diff_hotspots))) {
    row <- diff_hotspots[i,]
    
    traits <- raw_data %>%
      filter(qtl_chr == row$qtl_chr,
             qtl_pos >= (row$pos_mb - 0.5), # matches 1Mb window logic
             qtl_pos <= (row$pos_mb + 0.5)) %>%
      mutate(hotspot_region = paste0("Chr", row$qtl_chr, ":", row$pos_mb, "Mb"),
             diff_score = row$diff)
    
    mapping_traits[[i]] <- traits
  }
  
  return(bind_rows(mapping_traits) %>% distinct())
}

# ==========================================
# 7. Execution: Differential Analysis
# ==========================================

# Create new sub-directory
diff_dir <- file.path(output_dir, "Differential_Hotspots")
if (!dir.exists(diff_dir)) dir.create(diff_dir)

# --- DIET DIFFERENCES ---
diet_diff_traits <- find_differential_hotspots(diet_density, diet_data, threshold = 15)
write_csv(diet_diff_traits, file.path(diff_dir, "Diet_Specific_Hotspot_Traits.csv"))

# --- SEX DIFFERENCES ---
sex_diff_traits <- find_differential_hotspots(sex_density, sex_data, threshold = 15)
write_csv(sex_diff_traits, file.path(diff_dir, "Sex_Specific_Hotspot_Traits.csv"))

cat(paste("Differential analysis complete. Lists saved to:", diff_dir))

plot_differential_hotspots <- function(density_data, group_names, title, filename) {
  # Pivot to calculate the difference
  diff_plot_data <- density_data %>%
    pivot_wider(names_from = Group, values_from = count, values_fill = 0) %>%
    mutate(diff = .data[[group_names[1]]] - .data[[group_names[2]]],
           Dominant_Group = ifelse(diff > 0, group_names[1], group_names[2]))
  
  # Clean up chromosome levels
  diff_plot_data$qtl_chr <- factor(diff_plot_data$qtl_chr, 
                                   levels = c(1:19, "X", "Y", "M"))
  
  p <- ggplot(diff_plot_data, aes(x = pos_mb, y = diff, fill = Dominant_Group)) +
    geom_col(width = 0.5) + # Bars make the "magnitude" of difference very clear
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    facet_wrap(~qtl_chr, scales = "free_x", nrow = 2) +
    theme_minimal() +
    scale_fill_manual(values = c("Female" = "#e41a1c", "Male" = "#377eb8", 
                                 "High Carbohydrate" = "#4daf4a", "High Fat" = "#984ea3")) +
    labs(title = title,
         subtitle = paste(group_names[1], "(+) vs", group_names[2], "(-)"),
         x = "Genomic Position (Mb)",
         y = "Difference in Trans-eQTL Count",
         fill = "Higher Density In:") +
    theme(legend.position = "bottom",
          strip.background = element_rect(fill = "grey95"))
  
  ggsave(file.path(output_dir, filename), plot = p, width = 14, height = 7)
  return(p)
}

# --- Execute Visualization ---
# For Diet
plot_differential_hotspots(diet_density, 
                           group_names = c("High Carbohydrate", "High Fat"), 
                           title = "Differential eQTL Hotspots: HC vs HF", 
                           filename = "Differential_Hotspots_Diet.png")

# For Sex
plot_differential_hotspots(sex_density, 
                           group_names = c("Female", "Male"), 
                           title = "Differential eQTL Hotspots: Female vs Male", 
                           filename = "Differential_Hotspots_Sex.png")






