library(dplyr)
library(readr)

# Path to your data
folder_path <- "W:/General/Projects2/Jeeyeon Cha/NEW integration MAFA SNPs and DEGs/"

# 1. Load Original Files
c57_up <- read_csv(paste0(folder_path, "DEGs UP C57 Hets.csv"))
c57_dn <- read_csv(paste0(folder_path, "DEGs DOWN C57 Hets.csv"))
sjl_up <- read_csv(paste0(folder_path, "DEGs UP MixedSJL Hets.csv"))
sjl_dn <- read_csv(paste0(folder_path, "DEGs DOWN MixedSJL Hets.csv"))
gene_coords <- read_csv(paste0(folder_path, "mouse_genes_mm39.csv")) # THE MISSING PIECE

# 2. Standardize and Combine
c57_all <- bind_rows(c57_up %>% mutate(Dir_C57 = "UP"), c57_dn %>% mutate(Dir_C57 = "DOWN")) %>%
  rename(log2FC_C57 = starts_with("diffExpr"))
sjl_all <- bind_rows(sjl_up %>% mutate(Dir_SJL = "UP"), sjl_dn %>% mutate(Dir_SJL = "DOWN")) %>%
  rename(log2FC_SJL = starts_with("diffExpr"))

# 3. Join Strains
master <- full_join(c57_all, sjl_all, by = c("GeneId", "Symbol", "Gene Name"))

# 4. Refined Categorization Logic (The 4 Categories)
master <- master %>%
  mutate(
    Category = case_when(
      !is.na(Dir_C57) & !is.na(Dir_SJL) & Dir_C57 == Dir_SJL ~ "Shared",
      !is.na(Dir_C57) & !is.na(Dir_SJL) & Dir_C57 != Dir_SJL ~ "Discordant",
      !is.na(Dir_C57) & is.na(Dir_SJL)                       ~ "C57_Specific",
      is.na(Dir_C57) & !is.na(Dir_SJL)                       ~ "SJL_Specific",
      TRUE ~ "Unknown"
    ),
    # Unified Direction for the Plot Shape
    Direction = ifelse(!is.na(Dir_C57), Dir_C57, Dir_SJL)
  )

# 5. Merge Genomic Coordinates (CRITICAL STEP)
# Harmonize names to match what the coordinates file typically uses
gene_coords_sub <- gene_coords %>% 
  select(GeneId = `Gene stable ID`, Chr = `Chromosome`, Start = `Gene start (bp)`) %>%
  mutate(Chr = ifelse(grepl("Chr", Chr), Chr, paste0("Chr", Chr)))

final_master <- inner_join(master, gene_coords_sub, by = "GeneId")

# Save this as the new source for your Shiny App
write_csv(final_master, paste0(folder_path, "Master_DEGs_with_Coordinates.csv"))

print("Master file updated with 4 categories and coordinates.")