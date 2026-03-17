
library(tidyverse)
library(qtl2)
library(rstudioapi)

#markers <- read.csv(file.choose())
#chr_sep <- read.csv(file.choose())
#map <- readRDS(file.choose())
#combined_qtl <- readRDS(file.choose())


map=readRDS("W:\\General\\Projects2\\Holland500\\mapping_files\\Grcm39\\cleaned_DO500_pmap_grcm39.rds")
markers <- read_csv("W:\\General\\Projects2\\Holland500\\mapping_files\\Grcm39\\marker_annotations_DO500_grcm39.csv")
chr_sep <- read_csv("W:\\General\\Projects2\\Holland500\\mapping_files\\Grcm39\\chromosomal_sep_mm11.csv")

add_female = readRDS("W:\\General\\Projects2\\Klein\\WGCNA 500 DO cohort\\SkeletalMuscle\\Sex_specific_modules\\QTL_additive_SM_MEs_FEMALE.rds")
add_male = readRDS("W:\\General\\Projects2\\Klein\\WGCNA 500 DO cohort\\SkeletalMuscle\\Sex_specific_modules\\QTL_additive_SM_MEs_MALE.rds")


thr=6
peaks1 <- find_peaks(scan1_output = add_female, map = map, threshold = thr, prob=0.95)
print(peaks1)

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






thr=6
peaks1 <- find_peaks(scan1_output = add_male, map = map, threshold = thr, prob=0.95)
print(peaks1)

subset_modules <- unique(peaks1$lodcolumn)
sub_add_male <- data.frame(add_male)[subset_modules]

excludeme2 <- c("MEdarkgrey",
               "MEdarkorange",
               "MElightyellow",
               "MEorange",
               "MEskyblue",
               "MEwhite",
               "MEgrey")

sub_add_male2 <- sub_add_male[,-which(colnames(sub_add_male) %in% excludeme2)]


heatmap <- qtl_heatmap(qtl.temp = sub_add_male2, 
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

pdf(file = "heatmap_additive_male.pdf", width=9, height=4)
plot(heatmap)
dev.off()




