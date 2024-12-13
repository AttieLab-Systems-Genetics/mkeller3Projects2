# load libraries ======================================================================================================
require("rstudioapi")
require("bmediatR")
require("dplyr")
require("stringr")
require("tidyverse")
require("BiocManager")
require("intermediate")
require("ggplot2")
require("cowplot")
require("qtl2")
require("grid")
require("ggrepel")
require("gridGraphics")
require("ggpubr")
#require("shiny")
#require("shinyFiles")
require("bslib")
require("spsComps")
require("DT")
#require("shinyjs")
#require("shinycssloaders")
require("qtl2fst")
library("YandellIntermediate")


source("C:\\Lipid_mediations\\data\\source_data_loading_and_functionset_v2.R")

# set working directory
setwd(selectDirectory())




# CLINICAL TRAITS
#phenotypes <- read.csv("W:\\General\\Projects2\\Clinical_traits_DO1200\\Final in vivo DO data for analysis rankz.csv")
#rownames(phenotypes)<-phenotypes$mouse

# set QTL_list
#QTL_list <- read.csv("C:/Users/mkeller3/Desktop/Clinical_traits_QTL_additive.csv")

# prepare covariate matrix:
#covariate_statement = "~Sex*Diet+Gen_Lit"
#covar_matrix = model.matrix(as.formula(covariate_statement), data=phenotypes)[,-1]




# LIVER LIPIDS
#phenotypes <- read.csv("W:\\General\\Projects2\\Holland\\Liver_lipids_combat_badmiceremoved_allratios_rankz.csv")
#rownames(phenotypes)<-phenotypes$mouse

# set QTL_list
#QTL_list <- read.csv("C:/Users/mkeller3/Desktop/Liver_lipid_traits_QTL_additive.csv")

# prepare covariate matrix:
#covariate_statement = "~Sex*Diet+GenLit+Batch"
#covar_matrix = model.matrix(as.formula(covariate_statement), data=phenotypes)[,-1]




# PLASMA METABOLITES 13C-LABEL
phenotypes <- read.csv("W:\\General\\Projects2\\Kibbey\\D and 13C abeled molecules in plasma from 1200 DO.csv")
rownames(phenotypes)<-phenotypes$mouse

# set QTL_list
QTL_list <- read.csv("C:/Users/mkeller3/Desktop/13C_Labeled_metabolites_plasma_QTL_additive.csv")

# prepare covariate matrix:
covariate_statement = "~Sex*Diet+GenLit+C_batch"
covar_matrix = model.matrix(as.formula(covariate_statement), data=phenotypes)[,-1]






# mediator choice
mediator_object_choice <- "Isoforms"
#mediator_object_choice <- "Genes"







# set the phenotype list
pheno_list <- data.frame(matrix(nrow = nrow(QTL_list), ncol=1))
colnames(pheno_list) <- "Phenotype"
pheno_list$Phenotype <- QTL_list$lodcolumn

# mediator compartment
mediator_compartment <- "Liver"

# select mediators
chosen_mediators <- .GlobalEnv[[mediator_list_object[[mediator_object_choice]]]]
mediator_list <- make_mediator_list(mediator_object_choice, mediator_compartment)
if (mediator_object_choice=="RNA"){
  mediator_type <- "RNA/Protein"
}
if (mediator_object_choice=="Genes"){
  mediator_type <- "RNA/Protein"
  RNA_type <- "genes"
}
if (mediator_object_choice=="Isoforms"){
  mediator_type <- "RNA/Protein"
  RNA_type <- "isoforms"
}
if (mediator_object_choice=="Proteins"){
  mediator_type <- "RNA/Protein"
}
if (mediator_object_choice=="Phenotypes"){
  mediator_type <- "Phenotypes"
}

# set window
window = 4
if (window == 0){
  usewindow <- FALSE
}
if (window > 0){
  usewindow <- TRUE
  window <- as.numeric(window)
}

# set LOD thr
LOD_thr <- 7

# set lod drop and prob thresholds
probability_effect_thr <- 0.4
lod_drop_thr <- 0.4

# set DO set
DO_set <- "second"

# set rankz option- this determines whether to re-rankz your data (recommended)
rankz_set = FALSE

# set counter
cntr <- 0

# start loop

repeat{
  #increment counter
  cntr <- cntr+1
  
  # find the next phenotype
  pheno_list_subset <- subset(pheno_list, Phenotype == QTL_list$lodcolumn[cntr])
  pheno_list_subset <- data.frame(pheno_list_subset)
  pheno_list_subset <- unique(pheno_list_subset)
  colnames(pheno_list_subset) <- "Phenotype"
  
  
  # set counter 2
  cntr2 <- 0
  
  # start second loop 
  QTL_list2 <- subset(QTL_list, lodcolumn %in% pheno_list_subset$Phenotype)
  
  repeat{
    # increment second counter
    cntr2 <- cntr2+1
    
  # find the next peak
  QTL_list3 <- subset(QTL_list2, marker == QTL_list2$marker[cntr2])
  
  # run the function:
  completed_mediation <- Emfinger_mediation(
    wr_dir = getwd(),
    K = kinship,
    genoprobs = allele_probabilities,
    map = map,
    pheno = phenotypes,
    covar = covar_matrix,
    test_these = pheno_list_subset,
    mediator = chosen_mediators,
    cc_colors = "new",
    QTL_list = QTL_list3,
    mediator_list = mediator_list,
    plot_mediation = TRUE,
    save_mediation = FALSE,
    show_fits = FALSE,
    lod_drop_thr = lod_drop_thr,
    ln_prior_c = ln_prior_c,
    diagrams = diagrams,
    save_rescan=FALSE,
    probability_effect_thr = probability_effect_thr,
    save_plot_obj =FALSE,
    covar_pheno = NULL,
    covar_mediator = NULL,
    weights_all = NULL,
    weights_pheno = NULL,
    weights_mediator = NULL,
    pval = TRUE,
    pval_model_type = "mediation",
    LOD_thr = LOD_thr,
    mediator_type = mediator_type,
    markers = markers,
    chr_breaks = chrom_sep,
    DO_set = DO_set,
    pdf_name_prompt = FALSE,
    ncors = ncors,
    RNA_type = RNA_type,
    usewindow = usewindow,
    window = window,
    save_PDF_file = FALSE,
    perm_p = FALSE,
    combine_graphs = FALSE,
    rescan = FALSE,
    rankz_set = rankz_set
  )
  cat("running Yandell mediation\n")
  annot_done <- completed_mediation$annotation_selected
  if (mediator_object_choice=="Genes"){
    colnames(annot_done)[grep("gene.id",colnames(annot_done))]<-"id"
  }
  if (mediator_object_choice=="Isoforms"){
    colnames(annot_done)[grep("transcript.id",colnames(annot_done))]<-"id"
  }
  if (mediator_object_choice=="phenotypes"){
    colnames(annot_done)[grep("ID",colnames(annot_done))]<-"id"
  }
  if (mediator_object_choice=="Proteins"){
    colnames(annot_done)[grep("protein.id",colnames(annot_done))]<-"id"
  }
  #run Brian's function version
  yandell_mediation <- YandellIntermediate::mediation_scan(
    target = completed_mediation$pheno_subset,
    driver = completed_mediation$temporary_genoprobs,
    mediator = completed_mediation$selected_mediator_matrix,
    annotation = annot_done,
    covar = completed_mediation$subset_covariates,
    verbose = TRUE
  )
  
  out_table <- completed_mediation$tabular_results
  
  yandell_mediation$LR_median <- median(yandell_mediation$LR)
  
  if (mediator_object_choice=="Isoforms"){
    out_table <- out_table[,-(grep("[.]y",colnames(out_table)))]
    the_xs <- colnames(out_table)[grep("[.]x",colnames(out_table))]
    for (i in 1:length(the_xs)){
      the_xs[i]<-str_split(the_xs[i], pattern="[.]x")[[1]][1]
    }
    colnames(out_table)[grep("[.]x",colnames(out_table))]<-the_xs
    
    #join the two tables
    out_table <- yandell_mediation[c("ID","LR","LR_median")] %>% inner_join(.,out_table,by=c("ID"="ID"))
    assign("out_table_test", out_table, envir = .GlobalEnv)
    out_table <- out_table %>%
      relocate("symbol","ID","LR","original_drop","complete","colocal","chr","pos","median_drop")
  }
  if (mediator_object_choice=="Genes"){
    out_table <- out_table[,-(grep("ID.y",colnames(out_table)))]
    colnames(out_table)[grep("ID",colnames(out_table))]<-"ID"
    the_ys <- colnames(out_table)[grep("[.]y",colnames(out_table))]
    for (i in 1:length(the_ys)){
      the_ys[i]<-str_split(the_ys[i], pattern="[.]y")[[1]][1]
    }
    colnames(out_table)[grep("[.]y",colnames(out_table))]<-the_ys
    out_table <- out_table[,-(grep("[.]x",colnames(out_table)))] %>%
      relocate("original_drop","complete","colocal","chr","pos","median_drop")
    out_table <- yandell_mediation[c("ID","LR","LR_median")] %>% inner_join(.,out_table,by=c("ID"="ID"))
    out_table <- out_table %>%
      relocate("symbol","LR","original_drop","complete","colocal","chr","pos","median_drop")
  }
  
  out_table$trait <- pheno_list_subset$Phenotype
  out_table$marker <- QTL_list3$marker
  out_table$covars <- covariate_statement
  
  flname <- paste0(mediator_compartment,"_",mediator_object_choice,"_",pheno_list_subset$Phenotype[1],"_peak",cntr2,".csv")
  
  # assign to reactive object
  write.csv(out_table, file=flname)
  # end second loop
  if (cntr2 >= nrow(QTL_list2)){break}
  }
  
  # break loop
  if (cntr >= nrow(QTL_list)){break}
}
