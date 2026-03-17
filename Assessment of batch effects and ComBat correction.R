# -------------------------------------------------------------------------
# Part 1: Load Data and Perform ANOVA for Batch
# -------------------------------------------------------------------------

# Load necessary libraries
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("tidyr")) install.packages("tidyr")
if (!require("dplyr")) install.packages("dplyr")

library(ggplot2)
library(tidyr)
library(dplyr)

# 1. Load the dataset
# Update the path below to your local file location
file_path <- "C:/Users/mkeller3/Desktop"
# For the purpose of this script, ensure the file name matches your upload
data <- read.csv("20221206_DO_Batch_1to9_MPK_dietdays.csv", stringsAsFactors = FALSE)

# 2. Identify lipid columns
metadata_cols <- c("mouse", "Sex", "DOwave", "Batch", "diet_days")
lipid_cols <- setdiff(names(data), metadata_cols)

# Ensure Batch is treated as a factor for ANOVA
data$Batch <- as.factor(data$Batch)

# 3. Perform ANOVA for each lipid (Batch effect)
anova_results <- data.frame(Lipid = character(), P_Value = numeric(), stringsAsFactors = FALSE)

for (lipid in lipid_cols) {
  if (sd(data[[lipid]], na.rm = TRUE) == 0 || all(is.na(data[[lipid]]))) next
  formula_str <- paste0("`", lipid, "` ~ Batch")
  res_aov <- aov(as.formula(formula_str), data = data)
  p_val <- summary(res_aov)[[1]][["Pr(>F)"]][1]
  anova_results <- rbind(anova_results, data.frame(Lipid = lipid, P_Value = p_val))
}

# 4. Rank lipids by significance and SAVE RESULTS
anova_results <- anova_results %>%
  arrange(P_Value) %>%
  mutate(Adj_P_Value = p.adjust(P_Value, method = "BH"))

# Save the Batch ANOVA results to CSV
write.csv(anova_results, "Lipid_Batch_ANOVA_Results.csv", row.names = FALSE)

# 5. Identify the top 12 lipids to act as the MASTER LIST for both plots
top_lipids <- head(anova_results$Lipid, 12)

# 6. Generate Plot 1 (Before Correction)
plot_data_raw <- data %>%
  select(Batch, all_of(top_lipids)) %>%
  pivot_longer(cols = -Batch, names_to = "Lipid", values_to = "Intensity")

# CRITICAL: Force the order to match top_lipids significance
plot_data_raw$Lipid <- factor(plot_data_raw$Lipid, levels = top_lipids)

ggplot(plot_data_raw, aes(x = Batch, y = Intensity, fill = Batch)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.5) +
  facet_wrap(~ Lipid, scales = "free_y", ncol = 3) +
  theme_minimal() +
  labs(title = "Top 12 Lipids Most Influenced by Batch Effect (Pre-Correction)",
       subtitle = "Ordered by significance: Most significant in top-left",
       x = "Batch", y = "Measured Value") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")


# -------------------------------------------------------------------------
# Part 2: ComBat Batch Correction & Verification Plot
# -------------------------------------------------------------------------

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!require("sva")) BiocManager::install("sva")
library(sva)
library(tibble)

# 1. Prepare data for ComBat
lipid_matrix_raw <- t(as.matrix(data[, lipid_cols]))
colnames(lipid_matrix_raw) <- data$mouse

# 2. Run ComBat Correction
combat_edata <- ComBat(dat = lipid_matrix_raw, 
                       batch = data$Batch, 
                       par.prior = TRUE, 
                       prior.plots = FALSE)

# 3. Reconstruct the corrected dataframe
data_corrected <- as.data.frame(t(combat_edata)) %>%
  rownames_to_column("mouse") %>%
  left_join(data %>% select(mouse, Sex, DOwave, Batch), by = "mouse")

# 4. Prepare data for the Post-Correction plot
plot_data_corrected <- data_corrected %>%
  select(Batch, all_of(top_lipids)) %>%
  pivot_longer(cols = -Batch, names_to = "Lipid", values_to = "Intensity")

# CRITICAL: Force the order to match Plot 1 EXACTLY
plot_data_corrected$Lipid <- factor(plot_data_corrected$Lipid, levels = top_lipids)

# 5. Generate Plot 2 (After Correction)
ggplot(plot_data_corrected, aes(x = Batch, y = Intensity, fill = Batch)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.5) +
  facet_wrap(~ Lipid, scales = "free_y", ncol = 3) + 
  theme_minimal() +
  labs(title = "Top 12 Lipids Following ComBat Correction (Post-Correction)",
       subtitle = "Verified: Exact same lipids in the exact same positions as Plot 1",
       x = "Batch", y = "Corrected Intensity Value") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

# 6. Save the corrected data for QTL analysis
write.csv(data_corrected, "20221206_DO_Batch_1to9_MPK_ComBat_Corrected.csv", row.names = FALSE)



# -------------------------------------------------------------------------
# Part 3: ANOVA for Diet Duration (diet_days) Effect
# -------------------------------------------------------------------------

# This script uses the batch-corrected data to identify which lipids are 
# significantly influenced by the number of days on diet.

library(ggplot2)
library(tidyr)
library(dplyr)

# 1. Load the corrected data and original metadata
# We need to rejoin diet_days because it might not have been in the corrected file
desktop_path <- "C:/Users/mkeller3/Desktop"
data_corrected <- read.csv(file.path(desktop_path, "20221206_DO_Batch_1to9_MPK_ComBat_Corrected.csv"))
original_data <- read.csv(file.path(desktop_path, "20221206_DO_Batch_1to9_MPK_dietdays.csv"))

# Merge diet_days into the corrected dataset
data_eval <- data_corrected %>%
  left_join(original_data %>% select(mouse, diet_days), by = "mouse")

# 2. Identify lipid columns (all columns except metadata)
metadata_cols <- c("mouse", "Sex", "DOwave", "Batch", "diet_days")
lipid_cols <- setdiff(names(data_eval), metadata_cols)

# Ensure Batch is treated as a factor for ANOVA
data_eval$diet_days <- as.factor(data_eval$diet_days)

# 3. Perform ANOVA/Linear Regression for diet_days
diet_results <- data.frame(Lipid = character(), 
                           P_Value = numeric(), 
                           R_Squared = numeric(), 
                           stringsAsFactors = FALSE)

for (lipid in lipid_cols) {
  # Skip columns with no variance or all NAs
  if (sd(data_eval[[lipid]], na.rm = TRUE) == 0 || all(is.na(data_eval[[lipid]]))) next
  
  # Linear model: Intensity ~ diet_days
  # Since diet_days is continuous, this is effectively a regression
  formula_str <- paste0("`", lipid, "` ~ diet_days")
  fit <- lm(as.formula(formula_str), data = data_eval)
  
  # Extract P-value for the diet_days coefficient
  p_val <- summary(fit)$coefficients[2, 4]
  r_sq  <- summary(fit)$r.squared
  
  diet_results <- rbind(diet_results, data.frame(Lipid = lipid, P_Value = p_val, R_Squared = r_sq))
}

# 4. Adjust P-Values and Rank
diet_results <- diet_results %>%
  arrange(P_Value) %>%
  mutate(Adj_P_Value = p.adjust(P_Value, method = "BH"))

# Save the results
write.csv(diet_results, file.path(desktop_path, "Lipid_DietDays_ANOVA_Results_FACTOR.csv"), row.names = FALSE)

# 5. Visualize the Top 6 (Greatest Influence) and Bottom 6 (Least Influence)
# Greatest influence = smallest P-values
most_influenced <- head(diet_results$Lipid, 6)
# Least influence = largest P-values
least_influenced <- tail(diet_results$Lipid, 6)

# Combine for plotting
target_lipids <- c(most_influenced, least_influenced)

plot_data_diet <- data_eval %>%
  select(diet_days, all_of(target_lipids)) %>%
  pivot_longer(cols = -diet_days, names_to = "Lipid", values_to = "Intensity") %>%
  mutate(Group = ifelse(Lipid %in% most_influenced, "Greatest Influence", "Least Influence"))

# Set levels to maintain specific order (Most then Least)
plot_data_diet$Lipid <- factor(plot_data_diet$Lipid, levels = target_lipids)
plot_data_diet$Group <- factor(plot_data_diet$Group, levels = c("Greatest Influence", "Least Influence"))

ggplot(plot_data_diet, aes(x = diet_days, y = Intensity)) +
  geom_point(alpha = 0.4, color = "steelblue") +
  geom_smooth(method = "lm", color = "darkred", se = TRUE) +
  facet_wrap(Group ~ Lipid, scales = "free_y", ncol = 3) +
  theme_minimal() +
  labs(title = "Lipids with Greatest vs. Least Influence from Diet Duration",
       subtitle = "Comparing most significant regressions to the least significant",
       x = "Days on High Fat Diet", 
       y = "Corrected Intensity") +
  theme(strip.text = element_text(face = "bold"))

message("Analysis complete. Results saved to 'Lipid_DietDays_ANOVA_Results.csv'.")



# -------------------------------------------------------------------------
# Part 6 (OPTIMIZED): Manual Comparison with Clean Separation
# -------------------------------------------------------------------------

# 1. USER INPUT
target_lipid <- "FL36_1_SAHOA_peak_7" 

# 2. Data Preparation (Remains the same)
df_before <- data %>%
  select(Batch, Sex, Value = all_of(target_lipid)) %>%
  mutate(Status = "Raw (Before ComBat)")

df_after <- data_corrected %>%
  select(Batch, Sex, Value = all_of(target_lipid)) %>%
  mutate(Status = "Corrected (After ComBat)")

comparison_df <- rbind(df_before, df_after)
comparison_df$Status <- factor(comparison_df$Status, 
                               levels = c("Raw (Before ComBat)", "Corrected (After ComBat)"))

# 3. Generate Plot
p_compare <- ggplot(comparison_df, aes(x = Batch, y = Value, fill = Sex, color = Sex)) +
  
  # --- OPTIMIZATION: BOXPLOT BORDER ---
  # 'color = "black"' forces a thin black border; 'linewidth' controls the thickness
  geom_boxplot(outlier.shape = NA, alpha = 0.5, 
               color = "black", linewidth = 0.3, # <--- CHANGE BOX BORDER HERE
               position = position_dodge(width = 0.8)) +
  
  # --- OPTIMIZATION: DOT SIZE ---
  # 'size' controls the diameter of the jittered points
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), 
             size = 1.5, # <--- CHANGE DOT SIZE HERE
             alpha = 0.6) +
  
  facet_wrap(~ Status, scales = "fixed") + 
  
  scale_fill_manual(values = c("F" = "#F8766D", "M" = "#00BFC4")) +
  scale_color_manual(values = c("F" = "#F8766D", "M" = "#00BFC4")) +
  
  theme_minimal() +
  labs(title = paste("ComBat Correction Impact:", target_lipid),
       x = "Experimental Batch",
       y = "Intensity (Log Scale)") +
  
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "bottom",
    # --- OPTIMIZATION: PANEL SPACING ---
    panel.spacing = unit(2, "lines"), # <--- ADDS SPACE BETWEEN PLOTS
    axis.title = element_text(face = "bold")
  )

# 4. Display and Save
print(p_compare)
ggsave(filename = file.path(file_path, paste0("ComBat_Check_", target_lipid, ".png")), 
       plot = p_compare, width = 12, height = 6, dpi = 300)













# -------------------------------------------------------------------------
# Part 3: Evaluate Sex Effects on Corrected Data
# -------------------------------------------------------------------------

# 1. Perform ANOVA for Sex (Using Corrected Data)
anova_sex <- data.frame(Lipid = character(), P_Value = numeric(), stringsAsFactors = FALSE)

for (lipid in lipid_cols) {
  formula_sex <- paste0("`", lipid, "` ~ Sex")
  res_sex <- aov(as.formula(formula_sex), data = data_corrected)
  p_val <- summary(res_sex)[[1]][["Pr(>F)"]][1]
  anova_sex <- rbind(anova_sex, data.frame(Lipid = lipid, P_Value = p_val))
}

# 2. Save and Print Sex Results
anova_sex <- anova_sex %>% arrange(P_Value) %>% mutate(Adj_P_Value = p.adjust(P_Value, method = "BH"))
write.csv(anova_sex, "Lipid_Sex_ANOVA_Results.csv", row.names = FALSE)

print("Sex ANOVA complete. Results saved to Lipid_Sex_ANOVA_Results.csv")

# -------------------------------------------------------------------------
# Part 4: Plot Top 12 Lipids with Greatest Sex Effect
# -------------------------------------------------------------------------

# 1. Identify the top 12 lipids from the Sex ANOVA performed in Part 3
top_lipids_sex <- head(anova_sex$Lipid, 18)

# 2. Prepare the data for plotting
plot_data_sex <- data_corrected %>%
  select(Sex, all_of(top_lipids_sex)) %>%
  pivot_longer(cols = -Sex, names_to = "Lipid", values_to = "Intensity")

# 3. CRITICAL: Lock the order so the most significant sex effects appear first
plot_data_sex$Lipid <- factor(plot_data_sex$Lipid, levels = top_lipids_sex)

# 4. Generate the Plot
ggplot(plot_data_sex, aes(x = Sex, y = Intensity, fill = Sex)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.4) +
  facet_wrap(~ Lipid, scales = "free_y", ncol = 3) +
  # Using distinct colors for Sex (Female = Pinkish, Male = Teal)
  scale_fill_manual(values = c("F" = "#F8766D", "M" = "#00BFC4")) +
  theme_minimal() +
  labs(title = "Top 18 Lipids with Strongest Sexual Dimorphism",
       subtitle = "Analysis performed on batch-corrected data",
       x = "Sex",
       y = "Corrected Intensity Value") +
  theme(strip.text = element_text(face = "bold"),
        legend.position = "none")

# Optional: Save this specific plot
# ggsave("Top_12_Sex_Effect_Lipids.png", width = 10, height = 8)

# -------------------------------------------------------------------------
# Part 5: Identify and Plot Top 9 Female-High and Top 9 Male-High Lipids
# -------------------------------------------------------------------------

# 1. Calculate the mean difference for all lipids to determine direction
# (Female Mean - Male Mean)
sex_direction <- data_corrected %>%
  pivot_longer(cols = all_of(lipid_cols), names_to = "Lipid", values_to = "Intensity") %>%
  group_by(Lipid, Sex) %>%
  summarise(mean_val = mean(Intensity, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = Sex, values_from = mean_val) %>%
  mutate(diff = F - M)

# 2. Merge direction with ANOVA significance results
anova_sex_directed <- anova_sex %>%
  left_join(sex_direction %>% select(Lipid, diff), by = "Lipid")

# 3. Pull Top 9 where Female > Male (diff is positive, sorted by p-value)
top_9_female_high <- anova_sex_directed %>%
  filter(diff > 0) %>%
  arrange(P_Value) %>%
  head(9) %>%
  pull(Lipid)

# 4. Pull Top 9 where Male > Female (diff is negative, sorted by p-value)
top_9_male_high <- anova_sex_directed %>%
  filter(diff < 0) %>%
  arrange(P_Value) %>%
  head(9) %>%
  pull(Lipid)

# 5. Combine for a final list of 18
top_18_split <- c(top_9_female_high, top_9_male_high)

# 6. Prepare Plot Data
plot_data_sex_split <- data_corrected %>%
  select(Sex, all_of(top_18_split)) %>%
  pivot_longer(cols = -Sex, names_to = "Lipid", values_to = "Intensity")

# Lock the order: first 9 are Female-High, next 9 are Male-High
plot_data_sex_split$Lipid <- factor(plot_data_sex_split$Lipid, levels = top_18_split)

# 7. Generate the Plot
ggplot(plot_data_sex_split, aes(x = Sex, y = Intensity, fill = Sex)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 0.8, alpha = 0.3) +
  facet_wrap(~ Lipid, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c("F" = "#F8766D", "M" = "#00BFC4")) +
  theme_minimal() +
  labs(title = "Top 18 Sex-Influenced Lipids (Split by Direction)",
       subtitle = "Rows 1-3: Top 9 Higher in Females | Rows 4-6: Top 9 Higher in Males",
       x = "Sex", y = "Corrected Intensity") +
  theme(strip.text = element_text(size = 7, face = "bold"),
        legend.position = "none")

# 8. Save the directional results for reference
write.csv(anova_sex_directed, "Lipid_Sex_ANOVA_with_Direction.csv", row.names = FALSE)


# -------------------------------------------------------------------------
# Part 6: Volcano Plot of Sex Effects
# -------------------------------------------------------------------------

# 1. Use the 'anova_sex_directed' dataframe created in Part 5
# This dataframe already contains Lipid, P_Value, and diff (F-M)
volcano_data <- anova_sex_directed %>%
  mutate(log_p = -log10(P_Value),
         Significance = case_when(
           P_Value < 0.05 & diff > 0 ~ "Higher in Female",
           P_Value < 0.05 & diff < 0 ~ "Higher in Male",
           TRUE ~ "Not Significant"
         ))

# 2. Generate the Volcano Plot
ggplot(volcano_data, aes(x = diff, y = log_p, color = Significance)) +
  geom_point(alpha = 0.6, size = 3) +
  # Use the same color scheme as your boxplots
  scale_color_manual(values = c("Higher in Female" = "#F8766D", 
                                "Higher in Male" = "#00BFC4", 
                                "Not Significant" = "grey70")) +
  # Add labels for the top 5 most significant lipids in each direction
  geom_text(data = subset(volcano_data, Lipid %in% c(top_9_female_high[1:5], top_9_male_high[1:5])),
            aes(label = Lipid), 
            vjust = 1.5, size = 3, check_overlap = TRUE, color = "black") +
  theme_minimal() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  labs(title = "Volcano Plot: Sex Effects on Lipid Abundance",
       subtitle = "Calculated from ComBat-corrected data",
       x = "Effect Size (Female Mean - Male Mean)",
       y = "-log10(P-value)") +
  theme(legend.position = "bottom")

# Optional: Save the volcano plot
# ggsave("Lipid_Sex_Volcano_Plot.png", width = 8, height = 6)




