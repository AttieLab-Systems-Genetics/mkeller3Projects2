# ===============================================================
# Module 1: Data Preparation Logic
# ===============================================================

library(tidyverse)
library(data.table)


# ==============================================================================
# 1. SETUP PATHS & LOAD DATA (Updated to use full split-scan peaks)
# ==============================================================================
# New directory for broader additive peaks
input_dir <- "W:/General/main_directory/annotated_peak_summaries"
output_dir <- "C:/Users/mkeller3/Desktop/QTL_Differential_Analysis"
if(!dir.exists(output_dir)) dir.create(output_dir)

# Define the NEW files
files <- list(
  HC     = "DO1200_liver_metabolites_labeled_HC_mice_additive_peaks.csv",
  HF     = "DO1200_liver_metabolites_labeled_HF_mice_additive_peaks.csv",
  Female = "DO1200_liver_metabolites_labeled_female_mice_additive_peaks.csv",
  Male   = "DO1200_liver_metabolites_labeled_male_mice_additive_peaks.csv"
)

# BACKUP: Interaction-filtered files (Commented out)
# input_dir_bak <- "W:/General/main_directory/annotated_peak_summaries/qtlxcovar_peaks/qtlxcovar_peaks_in_splitby_additive_scans"
# files_bak <- list(
#   HC     = "DO1200_liver_metabolites_labeled_qtlxdiet_peaks_in_HF_mice_additive.csv",
#   HF     = "DO1200_liver_metabolites_labeled_qtlxdiet_peaks_in_HC_mice_additive.csv",
#   Female = "DO1200_liver_metabolites_labeled_qtlxsex_peaks_in_female_mice_additive.csv",
#   Male   = "DO1200_liver_metabolites_labeled_qtlxsex_peaks_in_male_mice_additive.csv"
# )

data_list <- lapply(files, function(f) fread(file.path(input_dir, f)))

# 3. DEFINE SPECIFICITY (5Mb Window)
# A QTL is "Shared" if the same metabolite has a peak on the same Chr within 5Mb in the comparison group
# Note, this definition does not factor allele signatures, only co-mapping to the same locus
# A further refinement of what truly constitutes "shared" archtiecture wouild be to assess allele signatures

get_specificity <- function(target_df, comp_df) {
  target_df <- target_df %>% mutate(qtl_status = "Specific")
  
  for(i in 1:nrow(target_df)) {
    shared <- comp_df %>%
      filter(phenotype == target_df$phenotype[i],
             qtl_chr == target_df$qtl_chr[i],
             abs(qtl_pos - target_df$qtl_pos[i]) <= 5)
    
    if(nrow(shared) > 0) {
      target_df$qtl_status[i] <- "Shared"
    }
  }
  return(target_df)
}

# Run the comparisons
hc_final <- get_specificity(data_list$HC, data_list$HF)
hf_final <- get_specificity(data_list$HF, data_list$HC)
female_final <- get_specificity(data_list$Female, data_list$Male)
male_final <- get_specificity(data_list$Male, data_list$Female)

# 1. Create a function to summarize the Specificity results
summarize_qtl_specificity <- function(df, group_label) {
  summary <- df %>%
    group_by(qtl_status) %>%
    summarise(
      Count = n(),
      Unique_Metabolites = n_distinct(phenotype),
      Avg_LOD = mean(qtl_lod, na.rm = TRUE)
    ) %>%
    mutate(Group = group_label) %>%
    dplyr::select(Group, qtl_status, Count, Unique_Metabolites, Avg_LOD)
  
  return(summary)
}

# 2. Compile the tables
qtl_summary_table <- rbind(
  summarize_qtl_specificity(hc_final, "HC Diet"),
  summarize_qtl_specificity(hf_final, "HF Diet"),
  summarize_qtl_specificity(female_final, "Female Sex"),
  summarize_qtl_specificity(male_final, "Male Sex")
)

# 3. View the table in the console
print(qtl_summary_table)

# 4. (Optional) Export this summary to your Desktop
write.csv(qtl_summary_table, 
          file.path(output_dir, "QTL_Specificity_Summary_Table.csv"), 
          row.names = FALSE)



# ==============================================================================
# MODULE 2: DIFFERENTIAL MANHATTAN (GRCm39 SANITIZED)
# ==============================================================================

library(tidyverse)
library(data.table)

# ==============================================================================
# 1. GLOBAL GRCm39 CONFIGURATION (Run once to lock coordinates)
# ==============================================================================
chr_lens <- c("1"=195.15, "2"=181.75, "3"=159.75, "4"=156.86, "5"=151.76, "6"=149.58, 
              "7"=144.99, "8"=130.13, "9"=124.36, "10"=130.53, "11"=121.97, "12"=120.09, 
              "13"=120.88, "14"=125.14, "15"=104.07, "16"=98.01, "17"=95.29, "18"=90.72, 
              "19"=61.42, "X"=169.48)

GLOBAL_MAP <- data.frame(Chr = names(chr_lens), Length = as.numeric(chr_lens)) %>%
  mutate(offset = lag(cumsum(Length), default = 0),
         center = offset + Length / 2,
         chr_col = rep(c("royalblue3", "orange2"), length.out = 20))

TOTAL_LEN <- sum(GLOBAL_MAP$Length)

# ==============================================================================
# 2. THE MASTER PLOTTING FUNCTION
# ==============================================================================
generate_clean_manhattan <- function(data_df, dataset_name, output_folder, pt_size = 7.0) {
  
  # A. Sanitize Data
  df_plot <- data_df %>%
    mutate(Chr = as.character(qtl_chr),
           Chr = ifelse(Chr == "20", "X", Chr)) %>%
    filter(Chr %in% names(chr_lens)) %>%
    inner_join(GLOBAL_MAP, by = "Chr") %>%
    mutate(BP_cum = offset + qtl_pos)
  
  # B. Calculate counts for subtitle only (No red text on plot)
  n_cropped <- sum(df_plot$qtl_lod > 100)
  
  # C. The Plot
  p <- ggplot(df_plot, aes(x = BP_cum, y = qtl_lod)) +
    # Layer 1: SHARED (Hollow Circles)
    geom_point(data = filter(df_plot, qtl_status == "Shared"),
               shape = 1, color = "black", size = pt_size, stroke = 1.0, alpha = 0.5) +
    # Layer 2: UNIQUE (Filled Circles)
    geom_point(data = filter(df_plot, qtl_status == "Specific"),
               aes(fill = chr_col), shape = 21, color = "black", 
               size = pt_size, stroke = 0.3, alpha = 0.9) +
    # Formatting
    scale_fill_identity() +
    scale_x_continuous(breaks = GLOBAL_MAP$center, labels = GLOBAL_MAP$Chr, 
                       limits = c(0, TOTAL_LEN), expand = c(0.005, 0.005)) +
    coord_cartesian(ylim = c(0, 100)) + # Crops at 100
    theme_classic(base_size = 18) +
    theme(
      legend.position = "top",
      panel.grid = element_blank(),
      axis.text.x = element_text(size = 10, angle = 0),
      axis.title = element_text(face = "bold")
    ) +
    labs(
      title = paste("QTL Manhattan:", dataset_name),
      subtitle = paste("Solid = Unique | Open = Shared |", n_cropped, "peaks with LOD > 100 hidden"),
      x = "Chromosome", 
      y = "LOD Score"
    )
  
  # Save and display
  ggsave(file.path(output_folder, paste0("Manhattan_", dataset_name, ".png")), 
         p, width = 16, height = 6, dpi = 300)
  return(p)
}

# ==============================================================================
# 3. EXECUTION
# ==============================================================================
# (Assuming hc_final, hf_final, female_final, male_final are already in memory)

output_dir <- "C:/Users/mkeller3/Desktop/QTL_Differential_Analysis"
if(!dir.exists(output_dir)) dir.create(output_dir)

p_hc  <- generate_clean_manhattan(hc_final, "HC_Diet", output_dir, pt_size = 7)
p_hf  <- generate_clean_manhattan(hf_final, "HF_Diet", output_dir, pt_size = 7)
p_fem <- generate_clean_manhattan(female_final, "Female_Sex", output_dir, pt_size = 7)
p_mal <- generate_clean_manhattan(male_final, "Male_Sex", output_dir, pt_size = 7)

# To see them all at once in the RStudio Plot Window:
print(p_hc)
print(p_hf)
print(p_fem)
print(p_mal)



# ==============================================================================
# Module 3. Multi-Class Specificity & Manhattan Analysis
# ==============================================================================

# 1. SETUP & PATHS

library(tidyverse)
library(data.table)

input_dir  <- "W:/General/main_directory/annotated_peak_summaries"
output_dir <- "C:/Users/mkeller3/Desktop/QTL_Differential_Analysis"
if(!dir.exists(output_dir)) dir.create(output_dir)

# Define the trait classes and groups
trait_classes <- c("liver_metabolites_labeled", "plasma_metabolites", "liver_lipids", "clinical_traits")
groups <- c("HC", "HF", "female", "male")

# Map files to a nested list: trait_class -> group
file_map <- list(
  liver_metabolites_labeled = list(
    HC="DO1200_liver_metabolites_labeled_HC_mice_additive_peaks.csv",
    HF="DO1200_liver_metabolites_labeled_HF_mice_additive_peaks.csv",
    female="DO1200_liver_metabolites_labeled_female_mice_additive_peaks.csv",
    male="DO1200_liver_metabolites_labeled_male_mice_additive_peaks.csv"
  ),
  plasma_metabolites = list(
    HC="DO1200_plasma_metabolites_HC_mice_additive_peaks.csv",
    HF="DO1200_plasma_metabolites_HF_mice_additive_peaks.csv",
    female="DO1200_plasma_metabolites_female_mice_additive_peaks.csv",
    male="DO1200_plasma_metabolites_male_mice_additive_peaks.csv"
  ),
  liver_lipids = list(
    HC="DO1200_liver_lipids_HC_mice_additive_peaks.csv",
    HF="DO1200_liver_lipids_HF_mice_additive_peaks.csv",
    female="DO1200_liver_lipids_female_mice_additive_peaks.csv",
    male="DO1200_liver_lipids_male_mice_additive_peaks.csv"
  ),
  clinical_traits = list(
    HC="DO1200_clinical_traits_HC_mice_additive_peaks.csv",
    HF="DO1200_clinical_traits_HF_mice_additive_peaks.csv",
    female="DO1200_clinical_traits_female_mice_additive_peaks.csv",
    male="DO1200_clinical_traits_male_mice_additive_peaks.csv"
  )
)

# ==============================================================================
# 2. SPECIFICITY LOGIC (Shared vs Unique)
# ==============================================================================
get_specificity <- function(df_a, df_b, dist_kb = 4000) {
  # Define Shared: Same phenotype on same Chr within 4Mb
  find_shared <- function(target, reference) {
    target %>%
      rowwise() %>%
      mutate(is_shared = any(reference$phenotype == phenotype & 
                               reference$qtl_chr == qtl_chr & 
                               abs(reference$qtl_pos - qtl_pos) <= (dist_kb/1000))) %>%
      ungroup() %>%
      mutate(qtl_status = ifelse(is_shared, "Shared", "Specific"))
  }
  
  list(a_final = find_shared(df_a, df_b), b_final = find_shared(df_b, df_a))
}

# ==============================================================================
# 3. PROCESSING LOOP
# ==============================================================================
for (t_class in trait_classes) {
  message("Processing Trait Class: ", t_class)
  
  # Load Data
  dats <- lapply(file_map[[t_class]], function(f) fread(file.path(input_dir, f)))
  
  # Perform Comparisons
  diet_comp <- get_specificity(dats$HC, dats$HF)
  sex_comp  <- get_specificity(dats$female, dats$male)
  
  final_list <- list(
    HC = diet_comp$a_final, HF = diet_comp$b_final,
    Female = sex_comp$a_final, Male = sex_comp$b_final
  )
  
  # Generate Plots for each group
  for (grp_name in names(final_list)) {
    # We use the 'generate_clean_manhattan' function defined in previous step
    generate_clean_manhattan(
      data_df = final_list[[grp_name]], 
      dataset_name = paste0(t_class, "_", grp_name), 
      output_folder = output_dir, 
      pt_size = 7
    )
  }
}


# ==============================================================================
# 1. SETUP NEW DIRECTORIES
# ==============================================================================
hotspot_dir <- file.path(output_dir, "Hotspot_Analysis")
if(!dir.exists(hotspot_dir)) dir.create(hotspot_dir)

# ==============================================================================
# 2. DATA INTEGRATION & CLUSTERING
# ==============================================================================
# 4Mb window for hotspot clustering
WINDOW_SIZE <- 4.0 

all_peaks <- list()

# Load and Combine all files into one master data frame
for (t_class in names(file_map)) {
  for (grp in names(file_map[[t_class]])) {
    file_path <- file.path(input_dir, file_map[[t_class]][[grp]])
    if(file.exists(file_path)) {
      df <- fread(file_path) %>%
        select(phenotype, qtl_chr, qtl_pos, qtl_lod) %>%
        mutate(trait_class = t_class, 
               group = grp,
               # Ensure X chromosome is standardized
               qtl_chr = as.character(qtl_chr),
               qtl_chr = ifelse(qtl_chr == "20", "X", qtl_chr))
      all_peaks[[paste0(t_class, "_", grp)]] <- df
    }
  }
}

master_peaks <- bind_rows(all_peaks)

# Cluster identification logic
identify_hotspots <- function(df, gap_limit) {
  df %>%
    arrange(group, qtl_chr, qtl_pos) %>%
    group_by(group, qtl_chr) %>%
    mutate(
      # Identify gaps between peaks > gap_limit
      diff_pos = qtl_pos - lag(qtl_pos, default = first(qtl_pos)),
      new_cluster = ifelse(diff_pos > gap_limit, 1, 0),
      cluster_id = cumsum(new_cluster)
    ) %>%
    group_by(group, qtl_chr, cluster_id) %>%
    summarise(
      Start_Mb = min(qtl_pos),
      End_Mb = max(qtl_pos),
      Total_Traits = n(),
      Traits_List = paste(unique(phenotype), collapse = "; "),
      # Count unique traits per class
      n_liver_metab = sum(trait_class == "liver_metabolites_labeled"),
      n_plasma_metab = sum(trait_class == "plasma_metabolites"),
      n_liver_lipids = sum(trait_class == "liver_lipids"),
      n_clinical = sum(trait_class == "clinical_traits"),
      .groups = "drop"
    )
}

# Apply clustering
all_clusters <- identify_hotspots(master_peaks, WINDOW_SIZE)

# ==============================================================================
# 3. APPLY HOTSPOT CRITERIA & RANKING
# ==============================================================================
# Criteria: 2+ traits in any class AND 1+ in any other class
final_hotspots <- all_clusters %>%
  rowwise() %>%
  mutate(
    counts = list(c(n_liver_metab, n_plasma_metab, n_liver_lipids, n_clinical)),
    has_2_plus = any(unlist(counts) >= 2),
    has_other = sum(unlist(counts) > 0) >= 2
  ) %>%
  filter(has_2_plus & has_other) %>%
  select(-counts, -has_2_plus, -has_other) %>%
  # Rank by total mapping traits
  arrange(desc(Total_Traits))

# ==============================================================================
# 4. EXPORT SUMMARY TABLE
# ==============================================================================
write.csv(final_hotspots, 
          file.path(hotspot_dir, "Global_Hotspot_Summary_Ranked.csv"), 
          row.names = FALSE)

message("Hotspot analysis complete. File saved to: ", hotspot_dir)


# ==============================================================================
# MODULE 4: INTEGRATED MULTI-CLASS MANHATTAN PLOT (RANKING & GRID FIX)
# ==============================================================================
library(ggplot2)
library(dplyr)
library(readr)
library(stringr)

# ==============================================================================
# 1. USER CONTROLS (THE TOGGLE)
# ==============================================================================
TARGET_PATH    <- "W:/General/main_directory/annotated_peak_summaries"

# --- TOGGLE RANKING HERE ---

RANK_BY        <- "Density" # Options: "Diversity" or "Density"
#RANK_BY        <- "Diversity" # Options: "Diversity" or "Density"

# ---------------------------

# Define filename based on ranking method
FILE_NAME <- paste0("Global_Hotspot_Manhattan_Ranked_by_", RANK_BY, ".png")
OUT_FILE  <- file.path(hotspot_dir, FILE_NAME)


BASE_FONT_SIZE <- 14    
PT_SIZE        <- 6     
LABEL_Y        <- 95    

# ==============================================================================
# 2. DATA LOADING & COORDINATES
# ==============================================================================
all_files <- list.files(path = TARGET_PATH, pattern = "additive_peaks.csv", full.names = TRUE)
target_files <- all_files[grepl("clinical|liver_lipids|liver_metabolites_labeled|plasma_metabolites", all_files)]
target_files <- target_files[!grepl("all_mice|genes|isoforms|splice", target_files)]

master_peaks <- do.call(rbind, lapply(target_files, function(f) {
  ctx <- case_when(
    str_detect(f, "_HC_mice_")     ~ "HC Diet",
    str_detect(f, "_HF_mice_")     ~ "HF Diet",
    str_detect(f, "_female_mice_") ~ "Female",
    str_detect(f, "_male_mice_")   ~ "Male",
    TRUE                           ~ "Other"
  )
  t_class <- case_when(
    str_detect(f, "clinical")      ~ "Clinical",
    str_detect(f, "liver_lipids")  ~ "Liver Lipids",
    str_detect(f, "labeled")       ~ "Liver Metab",
    str_detect(f, "plasma")        ~ "Plasma Metab"
  )
  read_csv(f, show_col_types = FALSE) %>%
    select(phenotype, qtl_chr, qtl_pos, qtl_lod) %>%
    mutate(trait_class = t_class, context = ctx, 
           qtl_chr = as.character(qtl_chr),
           qtl_chr = ifelse(qtl_chr == "20", "X", qtl_chr))
}))

chr_info <- data.frame(
  qtl_chr = c(as.character(1:19), "X"),
  length = c(195.1, 181.7, 159.5, 155.8, 152.0, 149.6, 144.9, 130.1, 124.3, 
             130.5, 121.9, 120.0, 120.9, 124.9, 104.0, 98.0, 95.0, 90.7, 61.3, 169.3)
) %>% mutate(offset = lag(cumsum(as.numeric(length)), default = 0))

plot_data <- master_peaks %>%
  inner_join(chr_info, by = "qtl_chr") %>%
  mutate(cum_pos = qtl_pos + offset)

class_colors   <- c("Clinical"="#984EA3", "Liver Lipids"="#4DAF4A", 
                    "Liver Metab"="#E41A1C", "Plasma Metab"="#377EB8")
context_shapes <- c("HC Diet"=21, "HF Diet"=24, "Female"=22, "Male"=23)

# ==============================================================================
# 3. RANKING LOGIC (STRICT TOGGLE)
# ==============================================================================
hotspots <- plot_data %>%
  arrange(qtl_chr, qtl_pos) %>%
  group_by(qtl_chr) %>%
  mutate(group = cumsum(c(1, diff(qtl_pos) > 2))) %>% 
  group_by(qtl_chr, group) %>%
  summarise(
    Total_Traits = n(),
    n_classes = n_distinct(trait_class),
    cum_center = cum_pos[which.max(qtl_lod)],
    pos_at_max = qtl_pos[which.max(qtl_lod)],
    max_lod    = max(qtl_lod),
    .groups = 'drop'
  ) 

top_10_hs <- if(RANK_BY == "Diversity") {
  hotspots %>% arrange(desc(n_classes), desc(Total_Traits)) %>% head(10)
} else {
  hotspots %>% arrange(desc(Total_Traits)) %>% head(10)
}

top_10_hs$label <- paste0("Chr", top_10_hs$qtl_chr, ":", round(top_10_hs$pos_at_max, 1), "Mb")

# ==============================================================================
# 4. THE PLOT (GRID AND LINE FIX)
# ==============================================================================
p <- ggplot(plot_data, aes(x = cum_pos, y = qtl_lod)) +
  # 1. Background Rectangles (Grey for Even Chrs)
  geom_rect(data = chr_info %>% filter(row_number() %% 2 == 0),
            aes(xmin = offset, xmax = offset + length, ymin = 0, ymax = Inf),
            fill = "grey95", color = NA, inherit.aes = FALSE) +
  
  # 2. Leader Lines (FOR ALL 10 LOCI)
  # Uniform solid black lines for all hotspots to ensure clear visual connection
  geom_segment(data = top_10_hs,
               aes(x = cum_center, xend = cum_center, 
                   y = max_lod + 2, yend = LABEL_Y - 3),
               color = "black", 
               linetype = "dotted", 
               linewidth = 0.75, 
               alpha = 1.0, # Increased alpha to 1.0 for true solid black
               inherit.aes = FALSE) + 
  
  # 3. Points
  geom_point(aes(fill = trait_class, shape = context), 
             size = PT_SIZE, color = "black", stroke = 0.3, alpha = 0.7) +
  
  # 4. Labels
  geom_text(data = top_10_hs, aes(x = cum_center, y = LABEL_Y, label = label), 
            angle = 90, vjust = -0.5, size = (BASE_FONT_SIZE * 0.3), 
            color = "white", fontface = "bold", inherit.aes = FALSE) + 
  geom_text(data = top_10_hs, aes(x = cum_center, y = LABEL_Y, label = label), 
            angle = 90, vjust = -0.5, size = (BASE_FONT_SIZE * 0.3), 
            color = "black", fontface = "bold", inherit.aes = FALSE) +
  
  # Scales & Legend
  scale_fill_manual(values = class_colors, name = "Trait Class:  ") +
  scale_shape_manual(values = context_shapes, name = "Context (Shape):  ") +
  scale_x_continuous(label = chr_info$qtl_chr, breaks = chr_info$offset + (chr_info$length/2), expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 105), breaks = seq(0, 100, 20)) +
  guides(fill = guide_legend(override.aes = list(shape = 21, size = 4)),
         shape = guide_legend(override.aes = list(fill = "grey70", size = 4))) +
  labs(title = "Integrated Multi-Class Manhattan Plot",
       subtitle = paste("Points shaped by Context | Ranked by", RANK_BY),
       x = "Chromosome (GRCm39 Mb)", y = "LOD Score") +
  theme_minimal(base_size = BASE_FONT_SIZE) +
  theme(
    panel.grid.major.x = element_blank(), # REMOVES THE VERTICAL LINES YOU SAW
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.box = "vertical"
  )

print(p)


# Save to the hotspot_dir
ggsave(OUT_FILE, p, width = 16, height = 8, dpi = 600)

message("Manhattan plot saved to: ", OUT_FILE)




# ==============================================================================
# MODULE 6: eQTL HOTSPOT ANALYSIS (DENSITY CALCULATION & SEPARATE PLOTTING)
# ==============================================================================
library(tidyverse)
library(data.table)

# 1. SETUP DIRECTORIES
eqtl_input_dir <- "W:/General/main_directory/annotated_peak_summaries"
base_output_dir <- "C:/Users/mkeller3/Desktop/QTL_Differential_Analysis"
eqtl_outdir <- file.path(base_output_dir, "eQTL_Hotspots")

if (!dir.exists(eqtl_outdir)) {
  dir.create(eqtl_outdir, recursive = TRUE)
}

# 2. DEFINE GRCm39 CHROMOSOME LENGTHS
grcm39_lengths <- tibble(
  qtl_chr = factor(c(1:19, "X"), levels = c(1:19, "X")),
  length_mb = c(195.15, 181.75, 159.62, 156.45, 151.19, 149.37, 144.99, 130.12, 124.36, 
                130.53, 121.97, 120.09, 120.88, 125.13, 104.07, 98.01, 95.29, 90.67, 61.42, 169.47)
)

# 3. DEFINE FILES
eqtl_files <- list(
  gene_HF = "DO1200_liver_genes_HF_mice_additive_peaks.csv",
  gene_HC = "DO1200_liver_genes_HC_mice_additive_peaks.csv",
  gene_Male = "DO1200_liver_genes_male_mice_additive_peaks.csv",
  gene_Female = "DO1200_liver_genes_female_mice_additive_peaks.csv",
  iso_HF = "DO1200_liver_isoforms_HF_mice_additive_peaks.csv",
  iso_HC = "DO1200_liver_isoforms_HC_mice_additive_peaks.csv",
  iso_Male = "DO1200_liver_isoforms_male_mice_additive_peaks.csv",
  iso_Female = "DO1200_liver_isoforms_female_mice_additive_peaks.csv"
)

# 4. LOAD eQTL DATA
message("Loading eQTL summary files...")
master_eqtl <- map_df(names(eqtl_files), function(f_key) {
  group_label <- str_split(f_key, "_")[[1]][2] 
  type_label  <- str_split(f_key, "_")[[1]][1] 
  file_path <- file.path(eqtl_input_dir, eqtl_files[[f_key]])
  if(!file.exists(file_path)) return(NULL)
  
  df <- fread(file_path) %>% as_tibble()
  
  possible_ids <- c("trait", "phenotype", "gene_id", "gene", "isoform", "id")
  found_id_col <- intersect(possible_ids, colnames(df))
  if(length(found_id_col) > 0) df <- df %>% rename_with(~ "trait", all_of(found_id_col[1]))
  
  if(!"qtl_chr" %in% colnames(df) && "chr" %in% colnames(df)) df <- df %>% rename(qtl_chr = chr)
  if(!"qtl_pos" %in% colnames(df) && "pos" %in% colnames(df)) df <- df %>% rename(qtl_pos = pos)
  
  df %>% mutate(Source_Group = group_label, 
                Source_Type = ifelse(type_label == "gene", "Gene", "Isoform"),
                qtl_chr = factor(qtl_chr, levels = c(1:19, "X")))
})

# 5. CALCULATE eQTL DENSITY
window_size <- 4; step_size <- 1
window_grid <- grcm39_lengths %>% rowwise() %>%
  do(data.frame(qtl_chr = .$qtl_chr, pos_mb = seq(0, ceiling(.$length_mb), by = step_size))) %>% ungroup()

eqtl_density <- master_eqtl %>%
  mutate(pos_mb = floor(qtl_pos)) %>%
  count(Source_Group, Source_Type, qtl_chr, pos_mb, cis, name = "count") %>%
  right_join(expand_grid(window_grid, Source_Group = unique(master_eqtl$Source_Group), 
                         Source_Type = unique(master_eqtl$Source_Type), cis = c(TRUE, FALSE)), 
             by = c("Source_Group", "Source_Type", "qtl_chr", "pos_mb", "cis")) %>%
  replace_na(list(count = 0))

# 6. PLOTTING FUNCTION WITH Y-LIMIT SUPPORT
plot_differential_by_type <- function(df, group_pair, target_type, is_cis = FALSE, y_limit = NULL) {
  type_prefix <- ifelse(is_cis, "Cis", "Trans")
  plot_df <- df %>% filter(cis == is_cis, Source_Group %in% group_pair, Source_Type == target_type)
  pair_data <- plot_df %>% pivot_wider(names_from = Source_Group, values_from = count, values_fill = 0)
  g1 <- group_pair[1]; g2 <- group_pair[2]
  
  plot_data <- pair_data %>%
    mutate(diff = .data[[g1]] - .data[[g2]],
           Dominant_Group = ifelse(diff >= 0, g1, g2),
           Fill_Key = paste(Dominant_Group, Source_Type))
  
  color_values <- c("Female Gene"="#D55E00", "Female Isoform"="#D55E00", "Male Gene"="#0072B2", 
                    "Male Isoform"="#0072B2", "HC Gene"="#E6AB02", "HC Isoform"="#E6AB02",
                    "HF Gene"="#66A61E", "HF Isoform"="#66A61E")
  
  p <- ggplot(plot_data, aes(x = pos_mb, y = diff, fill = Fill_Key)) +
    geom_col(width = step_size, color = "black", linewidth = 0.02) + 
    geom_hline(yintercept = 0, color = "grey20") +
    facet_grid(. ~ qtl_chr, scales = "free_x", space = "free_x") +
    theme_minimal() + scale_fill_manual(values = color_values) +
    labs(title = paste(target_type, type_prefix, "Hotspot Differences"),
         subtitle = paste(g1, "(+) vs", g2, "(-)"), x = "Mb", y = "Count Diff")
  
  if (!is.null(y_limit)) {
    p <- p + coord_cartesian(ylim = c(-y_limit, y_limit))
  }
  return(p)
}

# 7. CALCULATE GLOBAL MAX DIFFS FOR STANDARDIZED AXES
# Find max diff across Gene AND Isoform for Diet
max_diet <- eqtl_density %>%
  filter(cis == FALSE, Source_Group %in% c("HC", "HF")) %>%
  pivot_wider(names_from = Source_Group, values_from = count, values_fill = 0) %>%
  mutate(abs_diff = abs(HC - HF)) %>%
  pull(abs_diff) %>% max(na.rm = TRUE)

# Find max diff across Gene AND Isoform for Sex
max_sex <- eqtl_density %>%
  filter(cis == FALSE, Source_Group %in% c("Female", "Male")) %>%
  pivot_wider(names_from = Source_Group, values_from = count, values_fill = 0) %>%
  mutate(abs_diff = abs(Female - Male)) %>%
  pull(abs_diff) %>% max(na.rm = TRUE)

# SAVE eQTL PLOTS WITH STANDARDIZED Y-AXIS
for (stype in c("Gene", "Isoform")) {
  # Plot Diet with shared max_diet limit
  p_diet <- plot_differential_by_type(eqtl_density, c("HC", "HF"), stype, y_limit = max_diet)
  print(p_diet)
  ggsave(file.path(eqtl_outdir, paste0("Split_Trans_Diet_", stype, ".png")), p_diet, width = 20, height = 5)
  
  # Plot Sex with shared max_sex limit
  p_sex <- plot_differential_by_type(eqtl_density, c("Female", "Male"), stype, y_limit = max_sex)
  print(p_sex)
  ggsave(file.path(eqtl_outdir, paste0("Split_Trans_Sex_", stype, ".png")), p_sex, width = 20, height = 5)
}


# ==============================================================================
# MODULE 7: DIRECTIONAL CO-MAPPING & VISUALIZATION (COMPREHENSIVE LOOP)
# ==============================================================================
# Dependency: Requires 'eqtl_density', 'master_eqtl', and 'master_peaks' in environment.

library(tidyverse)
library(data.table)
library(ggrepel)
library(ggnewscale)

# 1. PREPARE PHYSIOLOGICAL DATA
if (exists("master_peaks")) {
  message("Starting Comprehensive Directional Co-mapping Analysis...")
  phys_data <- master_peaks %>% as_tibble()
  
  # Standardize trait and position columns
  possible_id_cols <- c("trait", "phenotype", "pheno", "id", "trait_id")
  found_id_col <- intersect(possible_id_cols, colnames(phys_data))
  if(length(found_id_col) > 0) phys_data <- phys_data %>% rename_with(~ "trait", all_of(found_id_col[1]))
  
  if("chr" %in% colnames(phys_data)) phys_data <- phys_data %>% rename(qtl_chr = chr)
  if("pos" %in% colnames(phys_data)) phys_data <- phys_data %>% rename(qtl_pos = pos)
  
  # 1.1 CLASSIFY TRAITS INTO 4 KEY CATEGORIES
  phys_data <- phys_data %>%
    mutate(trait_class = case_when(
      str_detect(trait, "(?i)pl_|plasma|serum|_pl|_P_|_P$|\\.P$|_plas|HDL|LDL|VLDL|ApoB|ApoA|Chol_pl|TG_pl|_CHOL") ~ "Plasma Metabolites",
      str_detect(trait, "(?i)TG|TAG|lipid|cholesterol|fatty|DG|CE|PC|PE|PS|PI|Sph|Cer|SM") ~ "Liver Lipids",
      str_detect(trait, "(?i)lv_|liver_|_lv|_L_|_L$|\\.L$") ~ "Liver Metabolites",
      str_detect(trait, "(?i)metab") & str_detect(trait, "(?i)pl|plasma|serum") ~ "Plasma Metabolites",
      str_detect(trait, "(?i)metab") & str_detect(trait, "(?i)lv|liver") ~ "Liver Metabolites",
      str_detect(trait, "(?i)glucose|insulin|weight|adipose|leptin|height|length|gluc|flow") ~ "Clinical",
      TRUE ~ "Clinical" 
    ))
  
  # Identify grouping column
  available_groups <- intersect(c("group", "Group", "source", "Source", "context"), colnames(phys_data))
  phys_group_col <- if(length(available_groups) > 0) available_groups[1] else NULL
  
  # 2. CALCULATE PHYSIOLOGICAL DENSITY
  phys_density_wide <- phys_data %>%
    mutate(!!sym(phys_group_col) := str_replace_all(!!sym(phys_group_col), " Diet", ""),
           pos_mb = floor(qtl_pos),
           qtl_chr = factor(qtl_chr, levels = levels(master_eqtl$qtl_chr))) %>%
    group_by(qtl_chr, pos_mb, !!sym(phys_group_col)) %>%
    summarise(phys_count = n(), 
              dominant_class = {
                counts <- table(trait_class)
                specific_classes <- counts[!names(counts) %in% c("Clinical", "None")]
                if(length(specific_classes) > 0) {
                  names(which.max(specific_classes))
                } else {
                  names(which.max(counts))
                }
              },
              .groups = 'drop') %>%
    pivot_wider(names_from = !!sym(phys_group_col), 
                values_from = c(phys_count, dominant_class), 
                values_fill = list(phys_count = 0, dominant_class = "None"))
  
  # 3. PIVOT eQTL DENSITY WIDE
  eqtl_density_wide <- eqtl_density %>%
    filter(cis == FALSE) %>%
    pivot_wider(names_from = Source_Group, values_from = count, 
                values_fill = 0, names_prefix = "eQTL_count_")
  
  # 4. JOIN
  final_alignment <- phys_density_wide %>%
    inner_join(eqtl_density_wide, by = c("qtl_chr", "pos_mb"))
  
  # GLOBAL SCALING CONSTANTS
  global_chr_limits <- master_eqtl %>%
    group_by(qtl_chr) %>%
    summarise(max_mb = max(qtl_pos, na.rm = TRUE) + 2, .groups = "drop")
  
  # Find global max eQTL count to fix the bubble size across all 4 plots
  global_max_eqtl <- max(eqtl_density$count[eqtl_density$cis == FALSE], na.rm = TRUE)
  
  # 5. GENERALIZED PLOTTING FUNCTION
  plot_directional_alignment <- function(df, group_a, group_b) {
    
    phys_a <- paste0("phys_count_", group_a)
    phys_b <- paste0("phys_count_", group_b)
    eqtl_a <- paste0("eQTL_count_", group_a)
    eqtl_b <- paste0("eQTL_count_", group_b)
    class_a <- paste0("dominant_class_", group_a)
    
    plot_df <- df %>%
      filter(!!sym(phys_a) > 0 | !!sym(eqtl_a) > 0) %>% 
      mutate(
        Agreement_Score = ((!!sym(phys_a) + 1) / (!!sym(phys_b) + 1)) * ((!!sym(eqtl_a) + 1) / (!!sym(eqtl_b) + 1)),
        # Clip score at 39.2 to sit clearly inside the 0-40 boundary
        Score_Clipped = pmin(Agreement_Score, 39.2),
        eQTL_Size = !!sym(eqtl_a),
        Trait_Class = factor(!!sym(class_a), levels = c("Clinical", "Liver Lipids", "Plasma Metabolites", "Liver Metabolites"))
      ) %>%
      filter(Agreement_Score > 1.2, !is.na(Trait_Class))
    
    if(nrow(plot_df) == 0) return(NULL)
    
    chr_bg <- global_chr_limits %>%
      mutate(fill_col = as.numeric(qtl_chr) %% 2 == 0)
    
    # Mimic Module 5 Color Palette
    class_colors <- c(
      "Clinical" = "#377EB8", 
      "Liver Lipids" = "#E41A1C", 
      "Plasma Metabolites" = "#4DAF4A", 
      "Liver Metabolites" = "#984EA3"
    )
    
    p <- ggplot(plot_df, aes(x = pos_mb, y = Score_Clipped)) +
      # Module 5 style background shading
      geom_rect(data = chr_bg, aes(xmin = 0, xmax = max_mb, ymin = -Inf, ymax = Inf, fill = fill_col), 
                inherit.aes = FALSE) +
      scale_fill_manual(values = c("TRUE" = "grey95", "FALSE" = "white"), guide = "none") +
      
      new_scale_fill() +
      # Shapes: 21 (Circle) for Gene, 24 (Triangle) for Isoform to match Module 5
      geom_point(aes(fill = Trait_Class, shape = Source_Type, size = eQTL_Size), 
                 color = "black", stroke = 0.3, alpha = 0.85) +
      
      facet_grid(. ~ qtl_chr, scales = "free_x", space = "free_x") +
      
      scale_fill_manual(values = class_colors, name = "Physiological Class", drop = FALSE) +
      scale_shape_manual(values = c("Gene" = 21, "Isoform" = 24), name = "Expression Source") +
      
      # Standardized bubble scaling
      scale_size_continuous(range = c(1.5, 10), limits = c(1, global_max_eqtl), name = "eQTL Density") +
      
      # Standardized Y-axis
      scale_y_continuous(limits = c(0, 40), breaks = seq(0, 40, 10), expand = c(0,0)) +
      
      theme_minimal() +
      labs(title = paste(group_a, "Dominant Co-mapping Hotspots"),
           subtitle = paste("Agreement Score relative to", group_b, "| Visual design matched to Module 5"),
           y = "Combined Agreement Score", x = "Genomic Position (Mb)") +
      theme(panel.spacing = unit(0.1, "lines"), 
            strip.text.x = element_text(size = 9, face = "bold"),
            strip.background = element_rect(fill = "grey90", color = NA),
            legend.position = "bottom", 
            legend.box = "vertical",
            legend.text = element_text(size = 9),
            axis.text.x = element_blank(), 
            axis.ticks.x = element_blank(),
            panel.grid.major.x = element_blank(), 
            panel.grid.minor.x = element_blank(),
            panel.grid.minor.y = element_blank(),
            panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.5)) +
      guides(fill = guide_legend(override.aes = list(shape = 21, size = 4, alpha = 1), order = 1),
             shape = guide_legend(override.aes = list(size = 4), order = 2),
             size = guide_legend(order = 3))
    
    return(p)
  }
  
  comparisons <- list(c("Female", "Male"), c("Male", "Female"), c("HF", "HC"), c("HC", "HF"))
  
  for(comp in comparisons) {
    grp_a <- comp[1]; grp_b <- comp[2]
    p_out <- plot_directional_alignment(final_alignment, grp_a, grp_b)
    if(!is.null(p_out)) {
      file_name <- paste0("Alignment_Map_", grp_a, "_Dominant.png")
      ggsave(file.path(eqtl_outdir, file_name), p_out, width = 18, height = 7)
      if(grp_a == "Female") print(p_out)
    }
  }
}


# ==============================================================================
# MODULE 8: GLOBAL HOTSPOT DISCOVERY (TOP 20) & TARGETED VISUALIZATION
# ==============================================================================
# Customized for environment with: phenotype, gene_symbol, transcript_symbol.
# This script identifies top hotspots and generates high-res zoom/BLUP plots.

library(tidyverse)
library(data.table)
library(ggrepel)
library(patchwork)

if (exists("master_peaks") && exists("final_alignment") && exists("master_eqtl")) {
  
  # ----------------------------------------------------------------------------
  # PART A: GLOBAL DISCOVERY (Runs every time to provide the menu)
  # ----------------------------------------------------------------------------
  message("Calculating Global Top 20 Agreement Hotspots...")
  
  discovery_base <- final_alignment %>%
    mutate(
      score_diet = ((phys_count_HC + 1) / (phys_count_HF + 1)) * ((eQTL_count_HC + 1) / (eQTL_count_HF + 1)),
      score_sex  = ((phys_count_Female + 1) / (phys_count_Male + 1)) * ((eQTL_count_Female + 1) / (eQTL_count_Male + 1))
    )
  
  top_diet <- discovery_base %>%
    group_by(qtl_chr, pos_mb) %>%
    summarise(Max_Score = max(score_diet), Top_Class = first(dominant_class_HC), .groups = "drop") %>%
    mutate(Effect_Type = "Diet")
  
  top_sex <- discovery_base %>%
    group_by(qtl_chr, pos_mb) %>%
    summarise(Max_Score = max(score_sex), Top_Class = first(dominant_class_Female), .groups = "drop") %>%
    mutate(Effect_Type = "Sex")
  
  global_top_20 <- bind_rows(top_diet, top_sex) %>%
    arrange(desc(Max_Score)) %>%
    mutate(Rank = row_number()) %>%
    select(Rank, qtl_chr, pos_mb, Max_Score, Top_Class, Effect_Type) %>%
    head(20)
  
  message("\n--- GLOBAL TOP 20 HOTSPOT DISCOVERY TABLE ---")
  print(as.data.frame(global_top_20), row.names = FALSE)
  
  # ----------------------------------------------------------------------------
  # PART B: TARGETED ZOOM PLOT FUNCTION
  # ----------------------------------------------------------------------------
  plot_targeted_hotspot <- function(chr, effect, mb_pos, buffer = 5) {
    
    if (effect == "Sex") { g_a <- "Female"; g_b <- "Male" } else { g_a <- "HC"; g_b <- "HF" }
    
    # Diagnostic confirmed 'phenotype' is the column name for phys traits
    phys_targeted <- master_peaks %>%
      as_tibble() %>%
      filter(qtl_chr == chr, qtl_pos >= (mb_pos - buffer), qtl_pos <= (mb_pos + buffer)) %>%
      rename(trait = phenotype) 
    
    all_classes <- c("Clinical", "Liver Lipids", "Plasma Metabolites", "Liver Metabolites")
    
    phys_targeted <- phys_targeted %>%
      mutate(trait_class = factor(case_when(
        str_detect(trait, "(?i)pl_|plasma|serum|_pl|_P_|_P$|\\.P$|_plas|HDL|LDL|VLDL|ApoB|ApoA|Chol_pl|TG_pl") ~ "Plasma Metabolites",
        str_detect(trait, "(?i)TG|TAG|lipid|cholesterol|fatty|DG|CE|PC|PE|PS|PI|Sph|Cer|SM") ~ "Liver Lipids",
        str_detect(trait, "(?i)lv_|liver_|_lv|_L_|_L$|\\.L$") ~ "Liver Metabolites",
        TRUE ~ "Clinical" 
      ), levels = all_classes))
    
    plot_data <- final_alignment %>%
      filter(qtl_chr == chr, pos_mb >= (mb_pos - buffer), pos_mb <= (mb_pos + buffer)) %>%
      mutate(
        p_a = !!sym(paste0("phys_count_", g_a)), p_b = !!sym(paste0("phys_count_", g_b)),
        e_a = !!sym(paste0("eQTL_count_", g_a)), e_b = !!sym(paste0("eQTL_count_", g_b)),
        Agreement_Score = ((p_a + 1) / (p_b + 1)) * ((e_a + 1) / (e_b + 1)),
        Score_Clipped = pmin(Agreement_Score, 39.2)
      )
    
    labeled_traits <- phys_targeted %>%
      mutate(pos_mb_bin = floor(qtl_pos)) %>%
      inner_join(plot_data %>% select(pos_mb, Score_Clipped, e_a, Source_Type), 
                 by = c("pos_mb_bin" = "pos_mb"), relationship = "many-to-many")
    
    hotspot_labels <- labeled_traits %>% filter(pos_mb_bin == mb_pos)
    class_colors <- c("Clinical" = "#377EB8", "Liver Lipids" = "#E41A1C", 
                      "Plasma Metabolites" = "#4DAF4A", "Liver Metabolites" = "#984EA3")
    
    p_zoom <- ggplot(labeled_traits, aes(x = qtl_pos, y = Score_Clipped)) +
      annotate("rect", xmin = mb_pos, xmax = mb_pos + 1, ymin = -Inf, ymax = Inf, fill = "yellow", alpha = 0.15) + 
      geom_point(aes(fill = trait_class, size = e_a, shape = Source_Type), color = "black", stroke = 0.5, alpha = 0.9) +
      geom_text_repel(data = hotspot_labels, aes(label = trait), size = 2.8, force = 8, box.padding = 0.8) +
      scale_fill_manual(values = class_colors, name = "Physiological Class", drop = FALSE) +
      scale_shape_manual(values = c("Gene" = 21, "Isoform" = 24)) +
      scale_size_continuous(range = c(2, 12), name = paste(g_a, "eQTL Density")) +
      scale_x_continuous(limits = c(mb_pos - buffer, mb_pos + buffer)) +
      scale_y_continuous(limits = c(0, 40)) +
      theme_minimal() +
      labs(title = paste0("Hotspot Zoom: Chr ", chr, " @ ", mb_pos, "Mb"),
           subtitle = paste("Condition:", effect),
           y = "Agreement Score", x = "Position (Mb)")
    
    return(p_zoom)
  }
  
  # ----------------------------------------------------------------------------
  # PART C: ALLELE EFFECT (BLUP) HEATMAP FUNCTION
  # ----------------------------------------------------------------------------
  plot_hotspot_blup_heatmap <- function(chr, effect, mb_pos) {
    
    if (effect == "Sex") { g_a <- "Female" } else { g_a <- "HC" }
    
    # 1. Physiological Extraction (Diagnostic confirms no A-H columns here)
    phys_subset <- master_peaks %>%
      filter(qtl_chr == chr, floor(qtl_pos) == mb_pos) %>%
      mutate(trait = phenotype, Type = "Physiological", Strain = "No Data", Effect = 0) %>%
      select(trait, Strain, Effect, Type)
    
    # 2. eQTL Extraction (Diagnostic confirmed gene_symbol/transcript_symbol/A-H exist)
    eqtl_subset_raw <- master_eqtl %>%
      filter(Source_Group == g_a, qtl_chr == chr, floor(qtl_pos) == mb_pos, cis == FALSE)
    
    founder_cols <- c("A", "B", "C", "D", "E", "F", "G", "H")
    eqtl_cols <- colnames(eqtl_subset_raw)
    
    if(nrow(eqtl_subset_raw) > 0) {
      eqtl_subset <- eqtl_subset_raw %>%
        mutate(trait = case_when(
          "gene_symbol" %in% eqtl_cols ~ as.character(gene_symbol),
          "transcript_symbol" %in% eqtl_cols ~ as.character(transcript_symbol),
          TRUE ~ as.character(gene_id)
        )) %>%
        select(trait, all_of(founder_cols)) %>%
        pivot_longer(cols = all_of(founder_cols), names_to = "Strain", values_to = "Effect") %>%
        mutate(Type = "eQTL")
    } else {
      eqtl_subset <- tibble()
    }
    
    heatmap_data <- bind_rows(phys_subset, eqtl_subset)
    
    if(nrow(heatmap_data) == 0) return(ggplot() + labs(title = "No traits found at hotspot center"))
    
    # Z-score normalization per trait (where data exists)
    heatmap_data <- heatmap_data %>%
      group_by(trait) %>%
      mutate(Effect_Z = if(all(Effect == 0)) 0 else as.vector(scale(Effect))) %>%
      ungroup()
    
    p_heat <- ggplot(heatmap_data, aes(x = Strain, y = trait, fill = Effect_Z)) +
      geom_tile(color = "white") +
      facet_grid(Type ~ ., scales = "free_y", space = "free_y") +
      scale_fill_gradient2(low = "#053061", mid = "white", high = "#67001f", midpoint = 0, name = "Z-Score") +
      theme_minimal() +
      labs(title = paste("Allele Effects @ Chr", chr, mb_pos, "Mb"),
           subtitle = "Physiological traits listed (allele effects missing in master_peaks)",
           x = "Founder Strain", y = "") +
      theme(axis.text.y = element_text(size = 7), 
            strip.text.y = element_text(face = "bold", angle = 0),
            panel.grid = element_blank())
    
    return(p_heat)
  }
  
  # ----------------------------------------------------------------------------
  # PART D: EXECUTION
  # ----------------------------------------------------------------------------
  # Define target (Example: Top Hotspot)
  target_chr <- "8"
  target_mb  <- 43
  target_eff <- "Diet"
  
  zoom_p <- plot_targeted_hotspot(target_chr, target_eff, target_mb)
  heat_p <- plot_hotspot_blup_heatmap(target_chr, target_eff, target_mb)
  
  # Final display using patchwork
  final_combined <- zoom_p / heat_p + plot_layout(heights = c(1, 1.2))
  print(final_combined)
  
} else {
  message("Error: Required data frames (master_peaks, final_alignment, master_eqtl) not found in R environment.")
}














# ==============================================================================
# MODULE 9: FUNCTIONAL ENRICHMENT PIPELINE (GO ANALYSIS)
# ==============================================================================
library(clusterProfiler)
library(org.Mm.eg.db)

clean_and_map_genes <- function(gene_list) {
  clean_ids <- gsub("^liver_", "", gene_list)
  clean_ids <- gsub("\\..*$", "", clean_ids)
  
  is_transcript <- any(grepl("ENSMUST", clean_ids))
  id_type <- if(is_transcript) "ENSEMBLTRANS" else "ENSEMBL"
  
  mapped <- bitr(clean_ids, 
                 fromType = id_type, 
                 toType = c("ENTREZID", "SYMBOL"), 
                 OrgDb = org.Mm.eg.db)
  return(mapped)
}

# Run for target hotspot
target_h = "Chr10_76"
h_info = phys_hotspots %>% filter(Hotspot_Label == target_h)
hotspot_egenes <- master_eqtl %>%
  filter(qtl_chr == h_info$qtl_chr, 
         qtl_pos >= h_info$peak_mb - 5, 
         qtl_pos <= h_info$peak_mb + 5,
         cis == FALSE)

enrichment_results <- map_df(unique(hotspot_egenes$Source_Group), function(grp) {
  map_df(unique(hotspot_egenes$Source_Type), function(typ) {
    raw_genes <- hotspot_egenes %>% filter(Source_Group == grp, Source_Type == typ) %>% pull(trait) %>% unique()
    if(length(raw_genes) < 5) return(NULL)
    
    mapped_df <- tryCatch({ clean_and_map_genes(raw_genes) }, error = function(e) NULL)
    if(is.null(mapped_df) || nrow(mapped_df) == 0) return(NULL)
    
    ego <- enrichGO(gene = mapped_df$ENTREZID, OrgDb = org.Mm.eg.db, ont = "BP", pvalueCutoff = 0.05, readable = TRUE)
    if(is.null(ego) || nrow(as.data.frame(ego)) == 0) return(NULL)
    
    as.data.frame(ego) %>% mutate(Source_Group = grp, Source_Type = typ)
  })
})

if(nrow(enrichment_results) > 0) {
  write_csv(enrichment_results, file.path(eqtl_outdir, paste0(target_h, "_GO_Enrichment_Results.csv")))
}















