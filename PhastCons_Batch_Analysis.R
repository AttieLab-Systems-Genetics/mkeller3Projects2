# Dependencies
library(data.table)
library(ggplot2)
library(GenomicScores)
library(phastCons35way.UCSC.mm39)
library(GenomicRanges)

message("Initializing PhastCons mapping script...")

# Directories
input_dir <- "/mnt/rdrive/mkeller3/General/main_directory/top_snps_all_impact/one_row_per_csq"
output_dir <- "/mnt/rdrive/mkeller3/General/Projects2/Diet_Sex dependent_ciseQTL_DO1200"

# Files
files <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE)

# PhastCons Obj
phast_obj <- getGScores("phastCons35way.UCSC.mm39")

all_results <- list()

for (f in files) {
  message("Processing file: ", basename(f))
  dt <- fread(f)
  
  if (!"variant_chr" %in% names(dt)) {
    stop("Missing variant_chr column in ", basename(f))
  }
  
  # Base coordinates directly mapped
  pos_bp <- as.numeric(dt$variant_pos) * 1e6
  clean_chr <- gsub("chr", "", as.character(dt$variant_chr), ignore.case = TRUE)
  
  gr <- GRanges(seqnames = paste0("chr", clean_chr), ranges = IRanges(start = pos_bp, width = 1))
  
  message(" -> Executing query against genomic reference...")
  dt$phastCons <- gscores(phast_obj, gr)$default
  dt[is.na(phastCons), phastCons := 0] # NA mapping
  
  # Remove the old SIFT score columns to avoid confusion
  if ("sift_score" %in% names(dt)) dt[, sift_score := NULL]
  if ("sift_prediction" %in% names(dt)) dt[, sift_prediction := NULL]
  
  # Extract plot data
  if ("impact" %in% names(dt)) {
    dt_plot <- dt[, .(file = basename(f), variant_id, impact, phastCons)]
  } else {
    dt_plot <- dt[, .(file = basename(f), variant_id, impact = NA_character_, phastCons)]
  }
  all_results[[basename(f)]] <- dt_plot
  
  # Output updated file 
  out_name <- file.path(output_dir, paste0(gsub("\\.csv", "", basename(f)), "_PhastConsScored.csv"))
  fwrite(dt, out_name)
  message(" -> Saved ", basename(out_name))
}

message("Compiling visual summaries...")
combined_dt <- rbindlist(all_results, fill=TRUE)

# Relevel impacts
impact_levels <- c("HIGH", "MODERATE", "MODIFIER", "LOW")
combined_dt[, impact := factor(impact, levels = impact_levels)]

# Visuals
pdf_path <- file.path(output_dir, "PhastCons_Summary_Plots.pdf")
pdf(pdf_path, width = 10, height = 8)

# 1. Histogram
p1 <- ggplot(combined_dt, aes(x = phastCons)) +
  geom_histogram(binwidth = 0.05, fill = "steelblue", color = "black") +
  theme_minimal() +
  labs(title = "Global Distribution of PhastCons Scores",
       subtitle = "Ranging from 0 (un-conserved) to 1 (perfectly conserved)",
       x = "PhastCons Score",
       y = "Frequency / Count")
print(p1)

# 2. Boxplot by Impact
p2 <- ggplot(combined_dt[!is.na(impact)], aes(x = impact, y = phastCons, fill = impact)) +
  geom_boxplot(outlier.alpha = 0.2) +
  scale_fill_manual(values = c("HIGH"="#ca0020", "MODERATE"="#f4a582", "MODIFIER"="#92c5de", "LOW"="#0571b0")) +
  theme_minimal() +
  labs(title = "PhastCons Conservation by GWAS Functional Impact",
       x = "Consequence Impact Level",
       y = "PhastCons Score") +
  theme(legend.position = "none")
print(p2)

# 3. Density mapping by dataset
p3 <- ggplot(combined_dt, aes(x = phastCons, fill = file)) +
  geom_density(alpha = 0.4) +
  theme_minimal() +
  labs(title = "PhastCons Distribution Comparison Across Datasets",
       x = "PhastCons Score",
       y = "Density") +
  theme(legend.position = "bottom", legend.text = element_text(size=6), legend.title = element_blank()) +
  guides(fill = guide_legend(ncol = 1))
print(p3)

dev.off()

message("All tasks completed. Visuals generated at: ", pdf_path)
