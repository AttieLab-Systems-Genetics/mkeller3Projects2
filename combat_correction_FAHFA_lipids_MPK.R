### Purpose: batch correction of phenotype
### Created: 2024-09-12

# sva::ComBat()
# https://rdrr.io/bioc/sva/man/ComBat.html


library(sva)
library(dplyr)
library(tidyr)
library(tibble)
library(data.table)

# load phenotype data
phe_raw <- fread("20221206_DO_Batch_1to9_MPK.csv", data.table=F)


## ComBat parametric adjustment for raw values
batch_data <- phe_raw %>% 
    mutate(Batch = factor(Batch)) %>%
    select(-DOwave, -Sex) %>% 
    remove_rownames() %>%
    column_to_rownames("Mouse")


combat_out <- sva::ComBat(t(batch_data[,-1,drop=F]), batch=batch_data$Batch) %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("Mouse")

combat_out <- phe_raw %>% 
    select(Mouse, DOwave, Sex, Batch) %>%
    left_join(combat_out, by = "Mouse")


fwrite(combat_out, file = "20221206_DO_Batch_1to9_MPK_ComBat.csv", sep = ",")

