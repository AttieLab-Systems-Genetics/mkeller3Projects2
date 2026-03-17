
library(tidyverse)
library(qtl2)
library(rstudioapi)

#markers <- read.csv(file.choose())
#chr_sep <- read.csv(file.choose())
#map <- readRDS(file.choose())
#combined_qtl <- readRDS(file.choose())



map=readRDS("W:\\General\\mapping_files_for_mark\\pmap_DO_diet_grcm39.rds")
markers <- read.csv("W:\\General\\mapping_files_for_mark\\marker_annotations_DO_diet_grcm39.csv")
chr_sep <- read_csv("W:\\General\\Projects2\\Holland500\\mapping_files\\Grcm39\\chromosomal_sep_mm11.csv")

clinical = readRDS("W:\\General\\main_directory\\annotated_peak_summaries\\qtlxcovar_interaction\\DO1200_clinical_traits_all_mice_qtlxdiet_lod_profiles.rds")



thr=10.5
peaks1 <- find_peaks(scan1_output = clinical, map = map, threshold = thr, prob=0.95)
print(peaks1)


########## Use this section to exclude any traits from the heatmap
subset_modules <- unique(peaks1$lodcolumn)
sub_add_female <- data.frame(add_female)[subset_modules]

excludeme1 <- c("MEdarkmagenta",
                "MEdarkolivegreen",
                "MEmagenta",
                "MEorange",
                "MEpaleturquoise",
                "MEroyalblue",
                "MEsienna3",
                "MEviolet",
                "MEyellow",
                "MEyellowgreen",
                "MEgrey")
sub_add_female2 <- sub_add_female[,-which(colnames(sub_add_female) %in% excludeme1)]



#### Draw heatmap

heatmap <- qtl_heatmap(qtl.temp = sub_add_female2, 
                       mrkrs = markers, 
                       low.thr = 6,
                       mid.thr = 8,
                       high.thr = 9, 
                       chr_breaks = chr_sep,
                       set2 = 0,
                       lowcolor = "white",
                       midcolor = 'darkred',
                       highcolor = 'black')

plot(heatmap)


pdf(file = "heatmap_additive_female.pdf", width=9, height=4)
plot(heatmap)
dev.off()






