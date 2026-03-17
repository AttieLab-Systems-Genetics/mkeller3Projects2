#' this is the loading data function required by the mediation app for the new DO
#' it also loads most of the microfunctions

# load data==================================================================================
# scans
mapping_paths <- read.csv("C:/Lipid_mediations/data/mapping_objects_locations.csv")
allele_probabilities <- readRDS(mapping_paths$object_path[which(mapping_paths$object_name=="allele_probabilities")])
map <- readRDS(mapping_paths$object_path[which(mapping_paths$object_name=="map")])
markers <- readRDS(mapping_paths$object_path[which(mapping_paths$object_name=="markers")])
kinship <- readRDS(mapping_paths$object_path[which(mapping_paths$object_name=="kinship")])
chrom_sep <- read.csv(mapping_paths$object_path[which(mapping_paths$object_name=="chrom_sep")])
diagrams <- readRDS(mapping_paths$object_path[which(mapping_paths$object_name=="diagrams")])

# datasets
data_paths <- read.csv("C:/Lipid_mediations/data/data_object_locations.csv")
all_liver_genes <- readRDS(data_paths$object_file_location[which(data_paths$object_name=="all_liver_genes")])
all_liver_isoforms <-  readRDS(data_paths$object_file_location[which(data_paths$object_name=="all_liver_isoforms")])
HC_diet_genes <-  readRDS(data_paths$object_file_location[which(data_paths$object_name=="HC_diet_genes")])
HF_diet_genes <- readRDS(data_paths$object_file_location[which(data_paths$object_name=="HF_diet_genes")])
clinical <- readRDS(data_paths$object_file_location[which(data_paths$object_name=="clinical")])
#Liver_genes <- readRDS("C:/Users/chris/Downloads/temp_mediation_set/temp_liver_genes.rds")

# process data==============================================================================
# mediation phenotypes
mediation_phenotypes <- list()
mediation_phenotypes$Data$Clinical <- clinical$data$raw
clnms_med_pheno <- data.frame(colnames(clinical$data$raw))
colnames(clnms_med_pheno)<-"data_name"
infomedpheno <- clinical$annot_phenotype[1:3]
infomedpheno <- infomedpheno[-c(1:5),]
clnms_med_pheno <- clnms_med_pheno %>% inner_join(., infomedpheno, by=c("data_name"="data_name"))
colnames(mediation_phenotypes$Data$Clinical)<-clnms_med_pheno$short_name
infomedpheno <- infomedpheno %>% arrange(., category)
mediation_phenotypes$Info$Clinical <- infomedpheno[2:3]
colnames(mediation_phenotypes$Info$Clinical)<-c("ID","Class")
mediation_phenotypes$Info$Clinical$Position <- rownames(mediation_phenotypes$Info$Clinical)
mediation_phenotypes$Info$Clinical$Position <- as.numeric(mediation_phenotypes$Info$Clinical$Position)

# mediation genes 
mediation_genes <- list()
mediation_genes$RNA_values$Liver <- all_liver_genes$data$norm
mediation_genes$RNA_info$Liver <- all_liver_genes$annot.mrna
colnames(mediation_genes$RNA_info$Liver)[1] <- "gene.id"
mediation_genes$RNA_info$Liver <- subset(mediation_genes$RNA_info$Liver, gene.id %in% colnames(mediation_genes$RNA_values$Liver))

# mediation isoforms
mediation_isoforms <- list()
mediation_isoforms$RNA_values$Liver <- all_liver_isoforms$data$norm
mediation_isoforms$RNA_info$Liver <- all_liver_isoforms$annot.mrna
colnames(mediation_isoforms$RNA_info$Liver)[3]<-"transcript.id"
colnames(mediation_isoforms$RNA_info$Liver)[1]<-"gene.id"
mediation_isoforms$RNA_info$Liver <- subset(mediation_isoforms$RNA_info$Liver, transcript.id %in% colnames(mediation_isoforms$RNA_values$Liver))

# set mediator list information
null_obj <- c()
mediator_list_object <- c( "mediation_genes","mediation_isoforms", "mediation_phenotypes")
names(mediator_list_object) <- c("Genes","Isoforms", "Phenotypes")
data_list <- c(data_paths$object_name)
names(data_list) <- data_paths$object_alias
categor_covar <- c("Sex","Wave","Batch","Generation","GenLit","Diet","DietName", "BirthDate", "Generation","sex","DOgen","DOwave","Gen_Lit", "Diet.name","diet")
current_env <- environment()
completed_mediation <- c()
covariates <- "~"

# set standards========================================================================
#DO_set
DO_set <- "second"
# peak type
peak_type <- "additive"
# ln_prior_c
ln_prior_c <- "complete"
# set cores
ncors = 1
# set whether to rerun scan
rescan = TRUE
# newscan
new_scan = NULL
# covar_statement
covar_statement <- "Sex"
# newscan_data 
newscans <- NULL
# newscan points
scan_data <- NULL
# set rescan rankz function
rankz_set <- TRUE
assign("rankz_test",rankz_set, envir = .GlobalEnv)

# source functions=====================================================================
source("C:/Lipid_mediations/data/Additive_mediator_script_for_shiny_v3.1_NewDO.R")
# reformat the phenotypes
pheno_reformat <- function(selected_dataset){
  # get the data to be mediated
  if (selected_dataset$datatype == "pheno"){
    pheno <- data.frame(selected_dataset$data$raw)
    ms_info <- data.frame(selected_dataset$annot.samples)
  }
  if (selected_dataset$datatype == "phenotype"){
    pheno <- data.frame(selected_dataset$data$raw)
    ms_info <- data.frame(selected_dataset$annot_samples)
  }
  if (selected_dataset$datatype != "pheno" && selected_dataset$datatype != "phenotype"){
    pheno <- data.frame(selected_dataset$data$norm)
    ms_info <- data.frame(selected_dataset$annot.samples)
  }
  pheno$Mouse <- rownames(pheno)
  # note- the RDS objects consider the covariate factors separately so these must be 
  # gathered from the sample information in the dataset
  colnames(ms_info)[1]<-"Mouse"
  rownames(ms_info)<-ms_info$Mouse
  pheno <- ms_info %>% inner_join(., pheno, by=c("Mouse"="Mouse"))
  rownames(pheno)<-pheno$Mouse
  return(pheno)
}

# generate QTL peak lists
make_QTL_list <- function(selected_dataset, chosen_trait){
  # populate lists
  QTL_list <- data.frame(selected_dataset$lod.peaks[["additive"]])
  # make sure the columns are named appropriately
  if (selected_dataset$datatype == "mrna"){
    if (selected_dataset$data_subtype == "genes"){
      QTL_list$lodcolumn <- QTL_list[,grep(pattern="^gene",x=colnames(QTL_list))]
    }
  }
  if (selected_dataset$datatype == "pheno"){
    QTL_list$lodcolumn <- QTL_list[,grep(pattern="^data",x=colnames(QTL_list))]
  }
  if (selected_dataset$datatype == "phenotype"){
    QTL_list$lodcolumn <- QTL_list[,grep(pattern="^data",x=colnames(QTL_list))]
  }
  if (selected_dataset$datatype == "protein"){
    QTL_list$lodcolumn <- QTL_list[,grep(pattern="^protein",x=colnames(QTL_list))]
  }
  # select the phenotype to scan within the peaks
  pheno_list <- data.frame(unique(QTL_list$lodcolumn))
  colnames(pheno_list)<-"Phenotype"
  # option selected is 'chosen_type'
  pheno_list <- subset(pheno_list, Phenotype == chosen_trait)
  QTL_list <- subset(QTL_list, lodcolumn == chosen_trait)
  # add the necessary information on the peak location
  if (DO_set == "first"){
    QTL_list <- markers[c("marker.id","chr","pos")] %>% inner_join(., QTL_list, by=c("marker.id"="marker.id"))
  }
  if (DO_set == "second"){
    markers <- markers %>% mutate(pos = bp_mm10/(10^6))
    colnames(markers)[1]<-"marker.id"
    if (selected_dataset$datatype == "phenotype"){
      colnames(QTL_list)[grep("^marker",colnames(QTL_list))]<-"marker.id"
      QTL_list <- markers[c("marker.id","chr")] %>% inner_join(., QTL_list, by=c("marker.id"="marker.id", "chr"="chr"))
    }
    if (selected_dataset$datatype != "phenotype"){
      colnames(QTL_list)[grep(pattern="^data",x=colnames(QTL_list))]<-"marker.id"
      colnames(QTL_list)[grep(pattern="^marker",x=colnames(QTL_list))]<-"marker.id"
      QTL_list <- markers[c("marker.id","chr","pos")] %>% inner_join(., QTL_list, by=c("marker.id"="marker.id"))
    }
  }
  QTL_list <- QTL_list %>% relocate("marker.id","chr","pos","lod")
  return(QTL_list)
}

# generate the mediation list object
make_mediator_list <- function(mediator_object_choice, mediator_compartment){
  chosen_set <- .GlobalEnv[[mediator_list_object[[mediator_object_choice]]]]
  if (mediator_object_choice=="Genes"){
    compartments <- names(chosen_set$RNA_values)
    major_object <- "RNA_values"
    mediator_type <- "RNA"
    major_info_object <- "RNA_info"
    overall_mediator_type <- "RNA/Protein"
    object_name <- mediator_compartment
    tissue_compartment <- mediator_compartment
    info_object <- mediator_compartment
    # select RNA type
    RNA_type <- "genes"
  }
  if (mediator_object_choice=="Isoforms"){
    compartments <- names(chosen_set$RNA_values)
    major_object <- "RNA_values"
    mediator_type <- "RNA"
    major_info_object <- "RNA_info"
    overall_mediator_type <- "RNA/Protein"
    object_name <- mediator_compartment
    tissue_compartment <- mediator_compartment
    info_object <- mediator_compartment
    # select RNA type
    RNA_type <- "isoforms"
  }
  if (mediator_object_choice=="Proteins"){
    compartments <- names(chosen_set$Protein_values)
    major_object <- "Protein_values"
    mediator_type <- "Proteins"
    major_info_object <- "Protein_info"
    overall_mediator_type <- "RNA/Protein"
    object_name <- mediator_compartment
    tissue_compartment <- mediator_compartment
    info_object <- mediator_compartment
  }
  if (mediator_object_choice=="Phenotypes"){
    compartments <- names(chosen_set$Data)
    major_object <- "Data"
    mediator_type <- "Phenotypes"
    major_info_object <- "Info"
    overall_mediator_type <- "Phenotypes"
    object_name <- mediator_compartment
    tissue_compartment <- mediator_compartment
    info_object <- mediator_compartment
  }
  mediator_list <- data.frame(matrix(nrow=1, ncol=6))
  colnames(mediator_list) <- c("Major_object","Object.names", "Tissue.compartment", "Type", "Major_info_object", "Info_object")
  mediator_list$Major_object[1] <- major_object
  mediator_list$Object.names[1] <- object_name
  mediator_list$Tissue.compartment[1] <- tissue_compartment
  mediator_list$Type[1] <- mediator_type
  mediator_list$Major_info_object[1] <- major_info_object
  mediator_list$Info_object[1] <- info_object
  return(mediator_list)
}

# plot the QTL if doing a new scan
QTL_plot_shiny <- function(qtl.temp, phenotype, covar_statement, LOD_thr, mrkrs){
  
  #get the rowname/colname info
  clnms <- colnames(qtl.temp)
  markers_ID <- rownames(qtl.temp)
  qtl.temp <- data.frame(qtl.temp)
  
  #derive the chromosomal location info from the markers
  #add markers
  qtl.temp$markers <- markers_ID
  clnms2 <- colnames(mrkrs)
  which_mrkr_clmn <-grep(pattern="marker", clnms2)
  which_pos_clnm <- grep(pattern="bp", clnms2)
  clnms2[which_mrkr_clmn]<-"markers"
  clnms2[which_pos_clnm]<-"position"
  colnames(mrkrs)<-clnms2
  mrkrs2 <- mrkrs[c("markers","chr","position")]
  mrkrs2$position <- mrkrs2$position/(10^6)
  
  qtl.temp <- qtl.temp %>%
    left_join(., mrkrs2, by=c("markers"="markers"))
  qtl.temp$order <- as.numeric(qtl.temp$chr)
  
  qtl.temp[which(qtl.temp$chr=="X"),"order"]<-20
  qtl.temp[which(qtl.temp$chr=="Y"),"order"]<-21
  qtl.temp[which(qtl.temp$chr=="M"),"order"]<-22
  qtl.temp$chr <- qtl.temp$order
  
  qtl_plot_obj <- qtl.temp %>% 
    
    # Compute chromosome size
    group_by(chr) %>% 
    summarise(chr_len=max(position)) %>% 
    
    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(chr_len)-chr_len) %>%
    select(-chr_len) %>%
    
    
    # Add this info to the initial dataset
    left_join(qtl.temp, ., by=c("chr"="chr")) %>%
    
    # Add a cumulative position of each SNP
    arrange(order, position) %>%
    mutate( BPcum=position+tot)
  
  #x axis mods: we do not want to display the cumulative 
  #position of SNP in bp, but just show the chromosome name instead.
  axisdf = qtl_plot_obj %>% group_by(order) %>% 
    summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
  axisdf$order[which(axisdf$order==20)]<-"X"
  axisdf$order[which(axisdf$order==21)]<-"Y"
  axisdf$order[which(axisdf$order==22)]<-"M"
  axisdf$order <- factor(axisdf$order, levels=c(as.character(1:19),"X","Y","M"))
  
  colnames(qtl_plot_obj)[1]<-"LOD"
  qtl_plot_obj$chr[which(qtl_plot_obj$chr==20)]<-"X"
  qtl_plot_obj$chr[which(qtl_plot_obj$chr==21)]<-"Y"
  qtl_plot_obj$chr[which(qtl_plot_obj$chr==22)]<-"M"
  
  #create the ggplot object
  #add label
  grob <- grobTree(textGrob(paste0(phenotype), x=0.1,  y=.95, hjust=0,
                            gp=gpar(col="red", fontsize=20, fontface="italic")))
  grob2 <- grobTree(textGrob(paste0("Covariates: ", paste0(covar_statement)), x=0.1,  y=.9, hjust=0,
                             gp=gpar(col="orange", fontsize=20, fontface="italic")))
  plot_QTL <- ggplot(qtl_plot_obj, aes(x=BPcum, y=LOD))+
    # Show all points
    geom_line(aes(color=as.factor(chr)), alpha=0.8, size=.5) +
    scale_color_manual(values = rep(c("black", "darkgrey"), 22 )) +
    
    # custom X axis:
    scale_x_continuous( label = axisdf$order, breaks= axisdf$center ) +
    #scale_y_continuous(expand = c(0, 0) ) +     
    # remove space between plot area and x axis
    ylim(0,c(max(qtl_plot_obj$LOD)+.25*max(qtl_plot_obj$LOD)))+
    
    # Custom the theme:
    theme_bw() +
    theme( 
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )+
    
    #change the axis lines
    ggplot2::theme(axis.line = element_line(colour = "black"))+
    
    #change the axis labels
    xlab("Chromosomal position")+
    ylab("LOD")+
    theme(axis.text=element_text(size=18),
          axis.title=element_text(size=20))+
    
    #add significant LOD
    geom_hline(aes(yintercept=LOD_thr,
                   linetype="Chosen LOD"), color="black")+
    
    #add annotations
    annotation_custom(grob)+
    annotation_custom(grob2)
  
  #show(plot_QTL)
  plots <- list(plot_QTL, qtl_plot_obj)
  return(plots)
}

# get allele effects
get_alleleEffects <- function(peaks, 
                              markers, 
                              geno, 
                              kinship, 
                              pheno, 
                              covar,
                              cores){
  if (DO_set == "second"){
    whch_mrkrid <- grep(pattern="marker",x=colnames(markers))
    colnames(markers)[whch_mrkrid]<-"marker.id"
    markers$pos <- markers$bp_mm10/(10^6)
  }
  
  peaks <- peaks %>% inner_join(., markers, by=c("chr"="chr", "pos"="pos"))
  
  peaks[,LETTERS[1:8]] <- 0
  for(i in 1:nrow(peaks)){
    gp <- geno[,peaks$chr[i]]
    #gp <- fst_extract(gp)
    gp[[1]] <- gp[[1]][,,peaks$marker.id[i], drop = FALSE]
    peaks[i, LETTERS[1:8]] <- scan1blup(genoprobs = gp,
                                        kinship   = kinship[[peaks$chr[i]]],
                                        pheno     = pheno[,peaks$lodcolumn[i], drop = FALSE],
                                        addcovar  = covar,
                                        cores = cores)[1, LETTERS[1:8]]
  }
  rm(gp)
  return(peaks)
} 

# rank transform the datasets
rankz = function(x) {
  x = rank(x, na.last = "keep", ties.method = "average") / (sum(!is.na(x)) + 1)
  return(qnorm(x))
}
