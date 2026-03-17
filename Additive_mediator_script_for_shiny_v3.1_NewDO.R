#' This wrapper is for mediating DO data using either the initial intermediate package or
#' the revised paradigm presented in the W. Crouse et al. 
#' "A Bayesian model selection approach to mediation analysis" 
#' 2023 PLoS Genetics paper. This requires the bmediatR package and the intermediate 
#' package
#' 
#' @author Chris Emfinger, PhD. A semi-accomplished mediator. probably
#' 
#' @usage Note that this requires the rownames of the variables analyzed (the
#' mediator, the QTL you're mediating, and the phenotypes, etc) must match
#' 
#' @param wr_dir is the working directory. Default ("choose") prompts the user
#' to select that directory. Otherwise, it is the character element showing the 
#' path to the working directory
#' 
#' @param K is the kinship matrix. If the attie_obj file is FALSE,
#' this assumes that the K is not loaded as part of it. If so, parameter
#' choices are "choose" for the 'file.choose()' command or 
#' the path is specified (as in for CHTC jobs).
#' 
#' @param genoprobs is the genome probabilities file. If the attie_obj file 
#' is FALSE, this assumes that the genoprobs is not loaded as part of it. If so, 
#' parameter choices are "choose" for the 'file.choose()' command or the source path must
#' be a character string
#' 
#' @param map is the map file. If the attie_obj file 
#' is FALSE, this assumes that the map is not loaded as part of it. If so, 
#' parameter choices are "choose" for the 'file.choose()' command or the source path must
#' be a character string
#' 
#' @param pheno is the csv file with the phenotyping information. It should have columns
#' matching the covariates to be used. If "choose" (default) it prompts the user to 
#' select the file. Otherwise, it is the path to the csv file. 
#' 
#' @param covar is the set of covariates used in the analysis. If "choose" it prompts the
#' user for the file. Otherwise it is a character vector defining the path to that file
#' 
#' @param mediator is the Rdata file containing the parameter for use as the mediator. 
#' it also selects a file that lists the mediators of interest. It will assume that the 
#' RDS file contains the entire set of desired mediators. This means if you want specific
#' genes you need a RDS file that keeps only the genes of interest
#' "choose" (default) it prompts the user to select the files. Otherwise it is a character
#' vector defining the path to that file. There's a specific structure to these
#' files. It is a list-type object where one major list 
#' 
#' @param cc_colors defines the colors to graph with ("old" is the classic DO, "new" is
#' the new color set we used for the eLife paper)
#' 
#' @param test_these defines the variables to mediate. If "choose" (default), the user
#' selects the csv files with the phenotypes to mediate. otherwise, this defines the 
#' path to the file listing these phenotypes
#' 
#' @param QTL_list defines the QTL for analysis. Should have the chromosme, LOD, and position
#' "choose" (default) allows the user to select the file with this info. Otherwise it assumes
#' it is a path to that variable
#' 
#' @param mediator_list defines the mediators within the mediator object to look at
#' if "choose" it allows the user to select this .csv file. If not, it is the path
#' to that file. 
#' 
#' @param plot_mediation plots the mediation. if "choose" the user selects the option. otherwise
#' TRUE will set the script to plot the mediation scan. FALSE will set the script not to.
#' 
#' @param save_mediation saves the mediation. if "choose" the user selects the option. otherwise
#' TRUE will set the script to save the mediation scan. FALSE will set the script not to. 
#' 
#' @param show_fits enables showing the QTL effects on target vs QTL effects on mediation
#' it is FALSE by default. 
#' 
#' @param lod_drop_thr is a number setting the fractional drop from the
#' MEDIAN LOD during the mediation analysis. I chose the MEDIAN because
#' it sometimes happens that the overall median LOD is significantly lower
#' than the original one. 
#' 
#' @param ln_prior_c defines the type of bmediatR comparison to use
#' for the moment it is defaulted to "reactive" to show everything. I'll
#' fix that in newer versions
#' 
#' @param diagrams is the location of the illustrations for the different
#' types of mediation. it is a RDS file and is used in making the graphs
#' 
#' @param probability_effect_thr is the parameter describing the threshold for
#' the probability effects to be considered of interest
#' 
#' @param save_rescan TRUE saves the rescanned QTL with auto names. "choose" lets
#' the user select a name
#' 
#' @param save_plot_obj TRUE saves the plot object for the pdf plot
#' 
#' @param covar_pheno NULL (default). If TRUE, it will assume the covariates
#' for the phenotypes are set by covar and will prompt the user to set the 
#' covariates for the mediator through the covar_mediator variable
#' 
#' @param covar_mediator NULL (default). if covar_pheno isn't NULL it will prompt
#' the user to select the file for the covariates for the mediator ("choose"), or is
#' a path to that covariate matrix file. 
#' 
#' @param weights_all NULL(default). this defines the matrix of weights if desired. 
#' If "choose", it will assume the weight assignments
#' are the same for both mediator and the phenotype variables
#' 
#' @param weights_pheno NULL(default). This defines the weights for the phenotype 
#' variable. if TRUE, it will assume that the "weights" variable is the weights for
#' the phenotype variable and will prompt the user to load the weights_mediator variable
#' 
#' @param weights_mediator NULL (default). if weights_pheno isn't NULL it will prompt
#' the user to select the file for the weights for the mediator ("choose"), or is
#' a path to that weights matrix file.  
#' 
#' @param pval TRUE (default). if TRUE and if ln_prior_c is not "reactive" it will
#' calculate the approximate and permuted pvalues for the bmediatR object
#' 
#' @param pval_model_type sets the parameter model_type for the get_approx_pval and
#' get_perm_pval functions of bmediatR. default is "mediation" see options for those
#' functions for more info. 
#' 
#' @param LOD_thr sets the LOD threshold for checking QTL rescans. Default is 6
#' 
#' @param perm_p determines whether to use permuted p-vals. this is only practical on
#' the CHTC because the size of the generated matrices and the time it takes to run each 
#' trial. It may also be buggy due to my not having been able to test it properly for these 
#' reasons
#' 
#' @param mediator_type either "RNA/Protein" or "Phenotypes"
#'
#' @param markers is the markers file. Used when converting QTL markers
#' to genomic coordinates. This assumes either "choose" where the user chooses
#' or that it is the character path to the file name. 
#' 
#' @param chr_breaks is the file with the chromosomal breaks in the genome
#' build you're using. Used to set the axes. Assumes either "choose" where
#' the user selects the file or otherwise assumes the variable is a character
#' string identifying the file's path. 
#' 
#' @param DO_set determines whether you're using the "first" DO (500 mice) or
#' the "second" DO (~1200 mice)
#'
#' @param pdf_name_prompt is a TRUE/FALSE. Because some terms can get really long,
#' the PDF maker will break if the file name is too long. I added this in to
#' let the user determine the pdf name if desired.
#' 
#' @param usewindow TRUE or FALSE. Limits the genes or proteins to those within a specific
#' window around the peak of interest in Mbp. default is FALSE
#' 
#' @param window numeric, in Mbp of window to either side of the peak to look. Default is
#' 4
#' 

#Define function==========================================================================
#-----------------------------------------------------------------------------------------
  
Emfinger_mediation <- function(
    wr_dir = "choose",
    K = "choose",
    genoprobs = "choose",
    map = "choose",
    pheno = "choose",
    covar = "choose",
    test_these = "choose",
    mediator = "choose",
    cc_colors = "new",
    QTL_list = "choose",
    mediator_list ="choose",
    plot_mediation = TRUE,
    save_mediation = TRUE,
    show_fits = FALSE,
    lod_drop_thr = 0.4,
    ln_prior_c = "complete",
    diagrams = "choose",
    save_rescan=FALSE,
    probability_effect_thr = 0.4,
    save_plot_obj =FALSE,
    covar_pheno = NULL,
    covar_mediator = NULL,
    weights_all = NULL,
    weights_pheno = NULL,
    weights_mediator = NULL,
    pval = TRUE,
    pval_model_type = "mediation",
    LOD_thr = 6,
    mediator_type = "RNA/Protein",
    markers = "choose",
    chr_breaks = "choose",
    DO_set = "first",
    pdf_name_prompt = FALSE,
    ncors = 20,
    RNA_type = "genes",
    usewindow = FALSE,
    window = 4,
    save_PDF_file = TRUE,
    perm_p = FALSE,
    combine_graphs = FALSE,
    rescan = TRUE,
    rankz_set = TRUE
){

#load libraries===========================================================================
  library("rstudioapi")
  start_time <- Sys.time()
  pckg_missing <- c("You are missing the following packages: ")
  numb_missing <- 0
  if (!require("BiocManager", quietly = TRUE)){
    pckg_missing <- c(pckg_missing, " BiocManager")
    numb_missing <- numb_missing + 1
  }
  if (!require("rstudioapi", quietly = TRUE)){
    pckg_missing <- c(pckg_missing, " rstudioapi")
    numb_missing <- numb_missing + 1
  }
  if (!require("bmediatR", quietly = TRUE)){
    pckg_missing <- c(pckg_missing, " bmediatR")
    numb_missing <- numb_missing + 1
  }
  if (!require("dplyr", quietly = TRUE)){
    pckg_missing <- c(pckg_missing, " dplyr")
    numb_missing <- numb_missing + 1
  }
  if (!require("tidyverse", quietly = TRUE)){
    pckg_missing <- c(pckg_missing, " tidyverse")
    numb_missing <- numb_missing + 1
  }
  if (!require("intermediate", quietly = TRUE)){
    pckg_missing <- c(pckg_missing, " intermediate")
    numb_missing <- numb_missing + 1
  }
  if (!require("ggplot2", quietly = TRUE)){
    pckg_missing <- c(pckg_missing, " ggplot2")
    numb_missing <- numb_missing + 1
  }
  if (!require("stringr", quietly = TRUE)){
    pckg_missing <- c(pckg_missing, " stringr")
    numb_missing <- numb_missing + 1
  }
  if (!require("cowplot", quietly = TRUE)){
    pckg_missing <- c(pckg_missing, " cowplot")
    numb_missing <- numb_missing + 1
  }
  if (!require("grid", quietly = TRUE)){
    pckg_missing <- c(pckg_missing, " grid")
    numb_missing <- numb_missing + 1
  }
  if (!require("qtl2", quietly = TRUE)){
    pckg_missing <- c(pckg_missing, " qtl2")
    numb_missing <- numb_missing + 1
  }
  if (!require("ggrepel", quietly = TRUE)){
    pckg_missing <- c(pckg_missing, " ggrepel")
    numb_missing <- numb_missing + 1
  }
  if (!require("gridGraphics", quietly = TRUE)){
    pckg_missing <- c(pckg_missing, " gridGraphics")
    numb_missing <- numb_missing + 1
  }
  if (!require("ggpubr", quietly = TRUE)){
    pckg_missing <- c(pckg_missing, " ggpubr")
    numb_missing <- numb_missing + 1
  }
  if (!require("YandellIntermediate", quietly = TRUE)){
    pckg_missing <- c(pckg_missing, " YandellIntermediate")
    numb_missing <- numb_missing + 1
  }
  #if (!require("magick", quietly = TRUE)){
  #  pckg_missing <- c(pckg_missing, " magick")
  #  numb_missing <- numb_missing + 1
  #}
  if (numb_missing == 0){
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
    require("YandellIntermediate")
  }
  if (numb_missing >= 1){
    print(paste0("ERROR: You need ", numb_missing, " packages."), quote = FALSE)
    print(pckg_missing, quote = FALSE)
    break
  }
  
  
#set working directory==================================================================
  if (wr_dir == "choose"){
    cat("Set your working directory...\n")
    setwd(selectDirectory())
    cat(paste0("...you're saving to ", getwd()))
  }
  
  if (wr_dir != "choose"){
    cat("...setting the working directory...\n")
    setwd(wr_dir)
    cat(paste0("...you're saving to ", getwd(), "\n"))
  } 

#define minifunctions==================================================================

#note some of these related to the manhattan-type plots use code modified from
  #https://r-graph-gallery.com/101_Manhattan_plot.html
  #others come directly from the bmediatR github
  
#define QTL plot function--------------------------------------------------------------
   QTL_plot <- function(qtl.temp){
     
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
    
    #derive the positional information from the rownames
    #for (i in 1:nrow(qtl.temp)){
    #  qtl.temp$chr[i] <- str_split(rownames(qtl.temp)[i], pattern="_")[[1]][1]
    #  qtl.temp$pos[i] <- as.numeric(str_split(rownames(qtl.temp)[i], pattern="_")[[1]][2])/(10^6)
    #  qtl.temp$order[i] <- as.numeric(str_split(rownames(qtl.temp)[i], pattern="_")[[1]][1])
    #}
    
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
    grob <- grobTree(textGrob(paste0("Tested QTL "), x=0.1,  y=.95, hjust=0,
                              gp=gpar(col="red", fontsize=20, fontface="italic")))
    grob2 <- grobTree(textGrob(paste0("Original position"), x=0.1,  y=.9, hjust=0,
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
      
      #add title
      labs(title = paste0(selected_QTL$lodcolumn[1],
                          " rescanned for common mice"))+
      
      #add significant LOD
      geom_hline(aes(yintercept=6,
                     linetype="Significant LOD"), color="black")+
      
      #add annotations
      annotation_custom(grob)+
      annotation_custom(grob2)
      
      #highlight the point where the original QTL was
      if (("marker.id" %in% colnames(selected_QTL))==TRUE){
        plot_QTL <- plot_QTL +
          geom_point(data=qtl_plot_obj[which(qtl_plot_obj$markers==selected_QTL$marker.id[cntr2]),], 
                                 color="orange", size=8) 
      }
      
      if (("marker.id" %in% colnames(selected_QTL))!=TRUE){
        plot_QTL <- plot_QTL +
      geom_point(data=qtl_plot_obj[which(qtl_plot_obj$markers==as.character(selected_QTL[cntr2,grep(colnames(selected_QTL), pattern="marker"), drop=FALSE])),], 
                              color="orange", size=8) 
      }
      #add significant LOD
    plot_QTL <- plot_QTL +
      geom_hline(aes(yintercept=selected_QTL$lod[cntr2],
                     linetype="Original LOD"), color="red")
    
    
    
    #show(plot_QTL)
    return(plot_QTL)
  }

#define coef plot functions------------------------------------------------------------
  effects_w_ci_plot <- function(qtl_effects, cc_colors, title_info){
    qtl_effects$strain <- c("A/J","B6/J","129","NOD","NZO","CAST","PWK","WSB")
    names(cc_colors)<-c("A/J","B6/J","129","NOD","NZO","CAST","PWK","WSB")
    plt_effect <- ggplot(qtl_effects, aes(x = target_effect, y = mediator_effect, color = strain))
    plt_effect <- plt_effect + scale_color_manual(values = cc_colors)
    plt_effect <- plt_effect + geom_point()
    # removes background
    plt_effect <- plt_effect + ggplot2::theme(panel.background = element_blank())
    # removes grid
    plt_effect <- plt_effect + ggplot2::theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    #add axis lines
    #add title
    plt_effect <- plt_effect + labs(title = title_info)
    #set the legend
    plt_effect <- plt_effect + labs(color = "Strain")
    plt_effect <- plt_effect + ggplot2::theme(axis.line = element_line(colour = "black"))
    #QTL effect 1 error
    plt_effect <- plt_effect + geom_errorbar(aes(xmin=target_effect-target_SE, xmax=target_effect+target_SE), width=.01)
    #QTL effect 2 error
    plt_effect <- plt_effect + geom_errorbar(aes(ymin=mediator_effect-mediator_SE, ymax=mediator_effect+mediator_SE), width=.01)
    plt_effect <- plt_effect +   xlab("QTL effects on Target")+
      ylab("QTL effects on Mediator")
    
    regrsn <- lm(formula = target_effect ~ mediator_effect,
                 data= qtl_effects)
    coeff<-coefficients(regrsn)          
    intercept<-coeff[1]
    slope<- coeff[2]
    
    plt_effect <- plt_effect + geom_abline(intercept = intercept, slope = slope, color="black", 
                                           linetype="dashed", size=1.0)
    grob <- grobTree(textGrob(paste0("Adj. R-squared: ", summary.lm(regrsn)[9]), x=0.1,  y=.95, hjust=0,
                              gp=gpar(col="black", fontsize=9, fontface="italic")))
    plt_effect <- plt_effect +  annotation_custom(grob)
    
    
    return(plt_effect)
  }
  
#define mediation LOD drop plotting function-------------------------------------------
  LOD_drop <- function(med, thr){
    #med is the intermediate mediation scan object with modifications
    #thr is the lod_drop_thr
    thr <- as.numeric(thr)
    drop_obj <- med
    
    if (mediator_type == "RNA/Protein"){
    drop_obj$order <- as.numeric(drop_obj$chr)
    
    drop_obj[which(drop_obj$chr=="X"),"order"]<-20
    drop_obj[which(drop_obj$chr=="Y"),"order"]<-21
    drop_obj[which(drop_obj$chr=="M"),"order"]<-22
    drop_obj$chr <- drop_obj$order
    }
    
    drop_obj <- drop_obj %>% 
      
      # Compute chromosome size
      group_by(chr) %>% 
      summarise(chr_len=max(pos)) %>% 
      
      # Calculate cumulative position of each chromosome
      mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
      select(-chr_len) %>%
      
      # Add this info to the initial dataset
      left_join(drop_obj, ., by=c("chr"="chr")) %>%
      
      # Add highlight and annotation information
      mutate( is_highlight_mediator=ifelse(median_drop > thr, "yes", "no")) %>%
      mutate( is_annotate_mediator=ifelse(median_drop > thr, "yes", "no")) %>%
      mutate( is_highlight_antimediator=ifelse(median_drop < (thr*(-1)), "yes", "no")) %>%
      mutate( is_annotate_antimediator=ifelse(median_drop < (thr*(-1)), "yes", "no")) %>%
      
      # Add a cumulative position of each SNP
      arrange(chr, pos) %>%
      mutate( BPcum=pos+tot)
    
    #x axis mods: we do not want to display the cumulative 
    #position of SNP in bp, but just show the chromosome name instead.
    if (mediator_type == "RNA/Protein"){
    axisdf = drop_obj %>% group_by(order) %>% 
      summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
    axisdf$order[which(axisdf$order==20)]<-"X"
    axisdf$order[which(axisdf$order==21)]<-"Y"
    axisdf$order[which(axisdf$order==22)]<-"M"
    axisdf$order <- factor(axisdf$order, levels=c(as.character(1:19),"X","Y","M"))
    
    drop_obj$chr[which(drop_obj$chr==20)]<-"X"
    drop_obj$chr[which(drop_obj$chr==21)]<-"Y"
    drop_obj$chr[which(drop_obj$chr==22)]<-"M"
    }
    if (mediator_type == "Phenotypes"){
      axisdf = drop_obj %>% group_by(chr) %>% 
        summarize(center=( max(pos) + min(pos) ) / 2 )
    }
    
    if (DO_set == "second"){
      colnames(drop_obj)[which(colnames(drop_obj)=="Gene.name")]<-"symbol"
    }

    #define Grobs
    grob <- grobTree(textGrob(paste0("Fractional drop greater than ",thr), x=0.1,  y=.95, hjust=0,
                              gp=gpar(col="orange", fontsize=20, fontface="italic")))
    
    grob2 <- grobTree(textGrob(paste0("Antimediator effect greater than ",thr), x=0.1,  y=.90, hjust=0,
                              gp=gpar(col="purple", fontsize=20, fontface="italic")))
    
    grob3 <- grobTree(textGrob(paste0("Original LOD"), x=0.7,  y=.95, hjust=0,
                              gp=gpar(col="black", fontsize=20, fontface="italic")))
    
    grob4 <- grobTree(textGrob(paste0("New median LOD"), x=0.7,  y=.90, hjust=0,
                               gp=gpar(col="red", fontsize=20, fontface="italic")))
    
    #create the ggplot object
    if (mediator_type=="RNA/Protein"){
    plot_LOD <- ggplot(drop_obj, aes(x=BPcum, y=LOD))
    }
    if (mediator_type=="Phenotypes"){
      plot_LOD <- ggplot(drop_obj, aes(x=pos, y=LOD))  
    }
    
    plot_LOD <- plot_LOD+
      # Show all points
      geom_point( aes(color=as.factor(chr)), alpha=0.8, size=1.3) +
      scale_color_manual(values = rep(c("lightgrey", "slategray4"), length(unique(drop_obj$chr)) )) +
      
      #show the selected points with greater change from median
      geom_point(data=subset(drop_obj, is_highlight_mediator=="yes"), color="orange", size=4) 
      
      
      # Add label using ggrepel to avoid overlapping
      if (mediator_type=="RNA/Protein"){
      plot_LOD <- plot_LOD +
      geom_label_repel(data = drop_obj[which(drop_obj$is_annotate_mediator=="yes"),], mapping = aes(label=drop_obj[which(drop_obj$is_annotate_mediator=="yes"),"symbol"]), size=8) 
      }
      if (mediator_type=="Phenotypes"){
      plot_LOD <- plot_LOD +
        geom_label_repel(data = drop_obj[which(drop_obj$is_annotate_mediator=="yes"),], mapping = aes(label=drop_obj[which(drop_obj$is_annotate_mediator=="yes"),"ID"]), size=8) 
      }
    
      #show the selected points with greater change from median
      plot_LOD <- plot_LOD +
      geom_point(data=subset(drop_obj, is_highlight_antimediator=="yes"), color="purple", size=4) 
      
      # Add label using ggrepel to avoid overlapping
      if (mediator_type=="RNA/Protein"){
        plot_LOD <- plot_LOD +
      geom_label_repel(data = drop_obj[which(drop_obj$is_annotate_antimediator=="yes"),], mapping = aes(label=drop_obj[which(drop_obj$is_annotate_antimediator=="yes"),"symbol"]), size=8) 
      }
      if (mediator_type=="Phenotypes"){
        plot_LOD <- plot_LOD +
          geom_label_repel(data = drop_obj[which(drop_obj$is_annotate_antimediator=="yes"),], mapping = aes(label=drop_obj[which(drop_obj$is_annotate_antimediator=="yes"),"ID"]), size=8) 
      }
      
      # custom X axis:
    if(mediator_type == "RNA/Protein"){
      plot_LOD <- plot_LOD+
      scale_x_continuous( label = axisdf$order, breaks= axisdf$center ) +
      #scale_y_continuous(expand = c(0, 0) ) +     
      # remove space between plot area and x axis
      ylim(0,c(max(drop_obj$LOD)+.25*max(drop_obj$LOD)))
    }
    if(mediator_type == "Phenotypes"){
      plot_LOD <- plot_LOD+
        scale_x_continuous( label = axisdf$chr, breaks= axisdf$center ) +
        #scale_y_continuous(expand = c(0, 0) ) +     
        # remove space between plot area and x axis
        ylim(0,c(max(drop_obj$LOD)+.25*max(drop_obj$LOD)))
    }
    
      # Custom the theme:
    plot_LOD <- plot_LOD +
      theme_bw() +
      theme( 
        legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
      )+
      
      #change the axis lines
      ggplot2::theme(axis.line = element_line(colour = "black"))
    
    if(mediator_type=="RNA/Protein"){
      plot_LOD <- plot_LOD + xlab("Chromosomal position")
    }
    if(mediator_type=="Phenotypes"){
      if (selected_mediator$Object.names[1]=="Clinical"){
      plot_LOD <- plot_LOD + xlab("Clinical phenotype")
      }
      if (selected_mediator$Object.names[1]=="Ex vivo"){
        plot_LOD <- plot_LOD + xlab("Islet ex vivo phenotype")
      }
      if (selected_mediator$Object.names[1]=="Lipids"){
        plot_LOD <- plot_LOD + xlab("Lipid class")
      }
      if (selected_mediator$Object.names[1]=="Metabolites"){
        plot_LOD <- plot_LOD + xlab("Metabolite class")
      }
      if (selected_mediator$Object.names[1]=="Modules"){
        plot_LOD <- plot_LOD + xlab("Module")
      }
      plot_LOD <- plot_LOD + theme (axis.text.x = element_text (angle = 45, hjust = 1))+
        theme(axis.text.x = element_text(size=10))
    }
      
      #change the axis labels
      plot_LOD <- plot_LOD +
      ylab("LOD")+
      theme(axis.text=element_text(size=18),
            axis.title=element_text(size=20))+
      
      #add title
      labs(title = paste0(selected_QTL$lodcolumn[cntr2], " QTL chr",
                          selected_chr,":",selected_pos,
                          " mediated by ",
                          selected_mediator$Tissue.compartment[1],
                          " ", selected_mediator$Type[1]))+
      
      #add original LOD
      geom_hline(aes(yintercept=selected_QTL$lod[cntr2],
                     linetype="Oringal LOD"), color="black")+
      
      #add median of the LOD drop values
      geom_hline(aes(yintercept=drop_obj$new_median[1],
                     linetype="Median new LOD"), color="red")+
      
      #add grob annotation
      annotation_custom(grob)+
      annotation_custom(grob2)+
      annotation_custom(grob3)+
      annotation_custom(grob4)+
      
      #add the legend
      labs(linetype ="LOD type")
    
    return(plot_LOD)
    
  }
  
#define the plotting of the bmediatR probabilities
  plot_bmediatR_odds <- function(post_odds, thr){
    post_odds_object <- post_odds
    #add order for the plotting:
    if(mediator_type == "RNA/Protein"){
    post_odds_object$order <- as.numeric(post_odds_object$chr)
    
    post_odds_object[which(post_odds_object$chr=="X"),"order"]<-20
    post_odds_object[which(post_odds_object$chr=="Y"),"order"]<-21
    post_odds_object[which(post_odds_object$chr=="M"),"order"]<-22
    post_odds_object$chr <- post_odds_object$order
    }
    
    post_odds_object <- post_odds_object %>% 
      
      # Compute chromosome size
      group_by(chr) %>% 
      summarise(chr_len=max(pos)) %>% 
      
      # Calculate cumulative position of each chromosome
      mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
      select(-chr_len) %>%
      
      # Add this info to the initial dataset
      left_join(post_odds_object, ., by=c("chr"="chr")) %>%
      
      # Add highlight and annotation information
      mutate( complete_prob = (exp(post_odds$complete))/(1+exp(post_odds$complete))) %>%
      mutate( is_highlight_complete=ifelse(complete_prob > thr, "yes", "no")) %>%
      mutate( is_annotate_complete=ifelse(complete_prob > thr, "yes", "no")) %>%
      
      # Add a cumulative position of each SNP
      arrange(chr, pos) %>%
      mutate( BPcum=pos+tot)
    
    #x axis mods: we do not want to display the cumulative 
    #position of SNP in bp, but just show the chromosome name instead.
    if (mediator_type=="RNA/Protein"){
    axisdf = post_odds_object %>% group_by(order) %>% 
      summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
    axisdf$order[which(axisdf$order==20)]<-"X"
    axisdf$order[which(axisdf$order==21)]<-"Y"
    axisdf$order[which(axisdf$order==22)]<-"M"
    axisdf$order <- factor(axisdf$order, levels=c(as.character(1:19),"X","Y","M"))
    
    post_odds_object$chr[which(post_odds_object$chr==20)]<-"X"
    post_odds_object$chr[which(post_odds_object$chr==21)]<-"Y"
    post_odds_object$chr[which(post_odds_object$chr==22)]<-"M"
    }
    
    if (mediator_type=="Phenotypes"){
      axisdf = post_odds_object %>% group_by(chr) %>% 
        summarize(center=( max(pos) + min(pos) ) / 2 )
    }
    
    #make sure the symbol is used:
     if (DO_set == "second"){
          colnames(post_odds_object)[which(colnames(post_odds_object)=="Gene.name")] <- "symbol"
        }

    #define Grobs
    grob <- grobTree(textGrob(paste0("Probability of complete mediation above ",thr), x=0.1,  y=.95, hjust=0,
                              gp=gpar(col="darkgreen", fontsize=20, fontface="italic")))
    
    #create the ggplot object
    if(mediator_type == "RNA/Protein"){
      plot_complete <- ggplot(post_odds_object, aes(x=BPcum, y=complete_prob))
    }
    if(mediator_type == "Phenotypes"){
      plot_complete <- ggplot(post_odds_object, aes(x=pos, y=complete_prob))
    }
    plot_complete <- plot_complete+
      # Show all points
      geom_point( aes(color=as.factor(chr)), alpha=0.8, size=1.3) +
      scale_color_manual(values = rep(c("lightgrey", "slategray4"), length(unique(post_odds_object$chr)) )) +
      
      #add threshold
      geom_hline(aes(yintercept=thr,
                     linetype="Threshold"), color="lightgreen")+
      
      #show the selected points with greater change from median
      geom_point(data=subset(post_odds_object, is_highlight_complete=="yes"), color="darkgreen", size=8) +
      
      #add highest probability
      geom_hline(aes(yintercept=1,
                     linetype="complete"), color="black")
      
      # Add label using ggrepel to avoid overlapping
      if (mediator_type=="RNA/Protein"){
        plot_complete <- plot_complete +
      geom_label_repel(data = post_odds_object[which(post_odds_object$is_annotate_complete=="yes"),], mapping = aes(label=post_odds_object[which(post_odds_object$is_annotate_complete=="yes"),"symbol"]), size=8) 
      }
      if (mediator_type=="Phenotypes"){
        plot_complete <- plot_complete +
      geom_label_repel(data = post_odds_object[which(post_odds_object$is_annotate_complete=="yes"),], mapping = aes(label=post_odds_object[which(post_odds_object$is_annotate_complete=="yes"),"ID"]), size=8) 
      }
        
      # custom X axis:
      if (mediator_type == "RNA/Protein"){
      plot_complete <- plot_complete +
      scale_x_continuous( label = axisdf$order, breaks= axisdf$center )
      }
    if (mediator_type == "Phenotypes"){
      plot_complete <- plot_complete +
        scale_x_continuous( label = axisdf$chr, breaks= axisdf$center )
    }
      #scale_y_continuous(expand = c(0, 0) ) +     
      # remove space between plot area and x axis
    plot_complete <- plot_complete +  
    ylim(0,1.25)+
      
      # Custom the theme:
      theme_bw() +
      theme( 
        legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
      )+
      
      #change the axis lines
      ggplot2::theme(axis.line = element_line(colour = "black"))
      
      #change the axis labels
    if (mediator_type == "RNA/Protein"){
      plot_complete <- plot_complete +
        xlab("Chromosomal position")
    }
    if (mediator_type == "Phenotypes"){
      if (selected_mediator$Object.names[1]=="Clinical"){
        plot_complete <- plot_complete + xlab("Clinical phenotype")
      }
      if (selected_mediator$Object.names[1]=="Ex vivo"){
        plot_complete <- plot_complete + xlab("Islet ex vivo phenotype")
      }
      if (selected_mediator$Object.names[1]=="Lipids"){
        plot_complete <- plot_complete + xlab("Lipid class")
      }
      if (selected_mediator$Object.names[1]=="Metabolites"){
        plot_complete <- plot_complete + xlab("Metabolite class")
      }
      if (selected_mediator$Object.names[1]=="Modules"){
        plot_complete <- plot_complete + xlab("Module type")
      }
      plot_complete <- plot_complete + theme (axis.text.x = element_text (angle = 45, hjust = 1))+
        theme(axis.text.x = element_text(size=10))
    }
      
    plot_complete <- plot_complete +
      ylab("Probability of complete mediation")+
      theme(axis.text=element_text(size=18),
            axis.title=element_text(size=20))+
      
      #add title
      labs(title = paste0(selected_QTL$lodcolumn[cntr2], " QTL chr",
                          selected_chr,":",selected_pos,
                          " mediated by ",
                          selected_mediator$Tissue.compartment[1],
                          " ", selected_mediator$Type[1]))+
      
      
      #add grob annotation
      annotation_custom(grob)+
      
      #add the legend
      labs(linetype =c("Threshold"))
    
    #show(plot_complete)
    
    #start the prep for the plots of the partial mediation
    
    post_odds_object <- post_odds_object %>%
      
      # Add highlight and annotation information
      mutate( partial_prob = (exp(post_odds$partial))/(1+exp(post_odds$partial))) %>%
      mutate( is_highlight_partial=ifelse(partial_prob > thr, "yes", "no")) %>%
      mutate( is_annotate_partial=ifelse(partial_prob > thr, "yes", "no")) 
    
    #define Grobs
    grob <- grobTree(textGrob(paste0("Probability of partial mediation above ",thr), x=0.1,  y=.95, hjust=0,
                              gp=gpar(col="aquamarine4", fontsize=20, fontface="italic")))
    
    #create the ggplot object
    if (mediator_type == "RNA/Protein"){
    plot_partial <- ggplot(post_odds_object, aes(x=BPcum, y=partial_prob))
    }
    if (mediator_type == "Phenotypes"){
      plot_partial <- ggplot(post_odds_object, aes(x=pos, y=partial_prob))
    }
    
      # Show all points
      plot_partial <- plot_partial + 
      geom_point( aes(color=as.factor(chr)), alpha=0.8, size=1.3) +
      scale_color_manual(values = rep(c("lightgrey", "slategray4"), length(unique(post_odds_object$chr)) )) +
      
      
      #add threshold
      geom_hline(aes(yintercept=thr,
                     linetype="Threshold"), color="aquamarine3")+
      
      #show the selected points with greater change from median
      geom_point(data=subset(post_odds_object, is_highlight_partial=="yes"), color="aquamarine4", size=2) +
      
      #add highest probability
      geom_hline(aes(yintercept=1,
                     linetype="complete"), color="black")
      
      # Add label using ggrepel to avoid overlapping
      if (mediator_type=="RNA/Protein"){  
      plot_partial <- plot_partial+
        geom_label_repel(data = post_odds_object[which(post_odds_object$is_annotate_partial=="yes"),], mapping = aes(label=post_odds_object[which(post_odds_object$is_annotate_partial=="yes"),"symbol"]), size=2) 
      }
      if (mediator_type=="Phenotypes"){  
        plot_partial <- plot_partial+
          geom_label_repel(data = post_odds_object[which(post_odds_object$is_annotate_partial=="yes"),], mapping = aes(label=post_odds_object[which(post_odds_object$is_annotate_partial=="yes"),"ID"]), size=2) 
      }
      
          
      # custom X axis:
        if (mediator_type == "RNA/Protein"){
          plot_partial <- plot_partial +
            scale_x_continuous( label = axisdf$order, breaks= axisdf$center )
        }
      if (mediator_type == "Phenotypes"){
        plot_partial <- plot_partial +
          scale_x_continuous( label = axisdf$chr, breaks= axisdf$center )
      }
     
      # remove space between plot area and x axis
      plot_partial <- plot_partial + ylim(0,1.25)+
      
      # Custom the theme:
      theme_bw() +
      theme( 
        legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
      )+
      
      #change the axis lines
      ggplot2::theme(axis.line = element_line(colour = "black"))
      
      #change the axis labels
        if (mediator_type == "RNA/Protein"){
          plot_partial <- plot_partial +
            xlab("Chromosomal position")
        }
      if (mediator_type == "Phenotypes"){
        if (selected_mediator$Object.names[1]=="Clinical"){
          plot_partial <- plot_partial + xlab("Clinical phenotype")
        }
        if (selected_mediator$Object.names[1]=="Ex vivo"){
          plot_partial <- plot_partial + xlab("Islet ex vivo phenotype")
        }
        if (selected_mediator$Object.names[1]=="Lipids"){
          plot_partial <- plot_partial + xlab("Lipid class")
        }
        if (selected_mediator$Object.names[1]=="Metabolites"){
          plot_partial <- plot_partial + xlab("Metabolite class")
        }
        if (selected_mediator$Object.names[1]=="Modules"){
          plot_partial <- plot_partial + xlab("Module type")
        }
        plot_partial <- plot_partial + theme (axis.text.x = element_text (angle = 45, hjust = 1))+
          theme(axis.text.x = element_text(size=10))
      }
        
      plot_partial <- plot_partial+
      ylab("Probability of partial mediation")+
      theme(axis.text=element_text(size=18),
            axis.title=element_text(size=20))+
      
      #add title
      labs(title = paste0(selected_QTL$lodcolumn[cntr2], " QTL chr",
                          selected_chr,":",selected_pos,
                          " mediated by ",
                          selected_mediator$Tissue.compartment[1],
                          " ", selected_mediator$Type[1]))+
      
      #add grob annotation
      annotation_custom(grob)+
      
      #add the legend
      labs(linetype =c("Threshold"))
    
    #show(plot_partial)
    
    #start the prep for the plots of the colocal effects
    
    post_odds_object <- post_odds_object %>%
      
      # Add highlight and annotation information
      mutate( colocal_prob = (exp(post_odds$colocal))/(1+exp(post_odds$colocal))) %>%
      mutate( is_highlight_colocal=ifelse(colocal_prob > thr, "yes", "no")) %>%
      mutate( is_annotate_colocal=ifelse(colocal_prob > thr, "yes", "no")) 
    
    #define Grobs
    grob <- grobTree(textGrob(paste0("Probability of colocal effects above ",thr), x=0.1,  y=.95, hjust=0,
                              gp=gpar(col="blue", fontsize=20, fontface="italic")))
    
    #create the ggplot object
    if(mediator_type == "RNA/Protein"){
    plot_colocal <- ggplot(post_odds_object, aes(x=BPcum, y=colocal_prob))
    }
    if(mediator_type == "Phenotypes"){
      plot_colocal <- ggplot(post_odds_object, aes(x=pos, y=colocal_prob))
    }
    
      # Show all points
    plot_colocal <- plot_colocal+
      geom_point( aes(color=as.factor(chr)), alpha=0.8, size=1.3) +
      scale_color_manual(values = rep(c("lightgrey", "slategray4"), length(unique(post_odds_object$chr)) )) +
      
      #add threshold
      geom_hline(aes(yintercept=thr,
                     linetype="Threshold"), color="lightblue")+
      
      #show the selected points with greater change from median
      geom_point(data=subset(post_odds_object, is_highlight_colocal=="yes"), color="blue", size=2) +
      
      #add highest probability
      geom_hline(aes(yintercept=1,
                     linetype="complete"), color="black")
      
      # Add label using ggrepel to avoid overlapping
      if (mediator_type == "RNA/Protein"){
        plot_colocal <- plot_colocal +
      geom_label_repel(data = post_odds_object[which(post_odds_object$is_annotate_colocal=="yes"),], mapping = aes(label=post_odds_object[which(post_odds_object$is_annotate_colocal=="yes"),"symbol"]), size=2) 
      }
    if (mediator_type == "Phenotypes"){
      plot_colocal <- plot_colocal +
        geom_label_repel(data = post_odds_object[which(post_odds_object$is_annotate_colocal=="yes"),], mapping = aes(label=post_odds_object[which(post_odds_object$is_annotate_colocal=="yes"),"ID"]), size=2) 
    }
        
      # custom X axis:
      if (mediator_type == "RNA/Protein"){
        plot_colocal <- plot_colocal +
          scale_x_continuous( label = axisdf$order, breaks= axisdf$center )
      }
    if (mediator_type == "Phenotypes"){
      plot_colocal <- plot_colocal +
        scale_x_continuous( label = axisdf$chr, breaks= axisdf$center )
    }
    plot_colocal <- plot_colocal +
      #scale_y_continuous(expand = c(0, 0) ) +     
      # remove space between plot area and x axis
      ylim(0,1.25)+
      
      # Custom the theme:
      theme_bw() +
      theme( 
        legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
      )+
      
      #change the axis lines
      ggplot2::theme(axis.line = element_line(colour = "black"))
      
      #change the axis labels
      if (mediator_type == "RNA/Protein"){
        plot_colocal <- plot_colocal +
          xlab("Chromosomal position")
      }
    if (mediator_type == "Phenotypes"){
      if (selected_mediator$Object.names[1]=="Clinical"){
        plot_colocal <- plot_colocal + xlab("Clinical phenotype")
      }
      if (selected_mediator$Object.names[1]=="Ex vivo"){
        plot_colocal <- plot_colocal + xlab("Islet ex vivo phenotype")
      }
      if (selected_mediator$Object.names[1]=="Lipids"){
        plot_colocal <- plot_colocal + xlab("Lipid class")
      }
      if (selected_mediator$Object.names[1]=="Metabolites"){
        plot_colocal <- plot_colocal + xlab("Metabolite class")
      }
      if (selected_mediator$Object.names[1]=="Modules"){
        plot_colocal <- plot_colocal + xlab("Module type")
      }
      plot_colocal <- plot_colocal + theme (axis.text.x = element_text (angle = 45, hjust = 1))+
        theme(axis.text.x = element_text(size=10))
    }
    
    plot_colocal <- plot_colocal+
      ylab("Probability of colocal effects")+
      theme(axis.text=element_text(size=18),
            axis.title=element_text(size=20))+
      
      #add title
      labs(title = paste0(selected_QTL$lodcolumn[cntr2], " QTL chr",
                          selected_chr,":",selected_pos,
                          " mediated by ",
                          selected_mediator$Tissue.compartment[1],
                          " ", selected_mediator$Type[1]))+
      
       
      #add grob annotation
      annotation_custom(grob)+
      
      #add the legend
      labs(linetype =c("Threshold"))
    
    #show(plot_colocal)  
    
    if(ln_prior_c == "reactive"){
    #plot the probability of reactivity
    
    #start the prep for the plots of the reactive effects
    
    post_odds_object <- post_odds_object %>%
      
      # Add highlight and annotation information
      mutate( reactive_prob = (exp(post_odds$reactive))/(1+exp(post_odds$reactive))) %>%
      mutate( is_highlight_reactive=ifelse(reactive_prob > thr, "yes", "no")) %>%
      mutate( is_annotate_reactive=ifelse(reactive_prob > thr, "yes", "no")) 
    
    #define Grobs
    grob <- grobTree(textGrob(paste0("Probability of any reactive effects above ",thr), x=0.1,  y=.95, hjust=0,
                              gp=gpar(col="goldenrod", fontsize=20, fontface="italic")))
    
    #create the ggplot object
    if(mediator_type == "RNA/Protein"){
      plot_reactive <- ggplot(post_odds_object, aes(x=BPcum, y=reactive_prob))
    }
    if(mediator_type == "Phenotypes"){
      plot_reactive <- ggplot(post_odds_object, aes(x=pos, y=reactive_prob))
    }
    
    
    plot_reactive <- plot_reactive+
      # Show all points
      geom_point( aes(color=as.factor(chr)), alpha=0.8, size=1.3) +
      scale_color_manual(values = rep(c("lightgrey", "slategray4"), length(unique(post_odds_object$chr)) )) +
      
      #add threshold
      geom_hline(aes(yintercept=thr,
                     linetype="Threshold"), color="gold")+
      
      #show the selected points with greater change from median
      geom_point(data=subset(post_odds_object, is_highlight_reactive=="yes"), color="goldenrod", size=2) +
      
      #add highest probability
      geom_hline(aes(yintercept=1,
                     linetype="complete"), color="black")
      
      # Add label using ggrepel to avoid overlapping
    if(mediator_type == "RNA/Protein"){
      plot_reactive <- plot_reactive +
      geom_label_repel(data = post_odds_object[which(post_odds_object$is_annotate_reactive=="yes"),], mapping = aes(label=post_odds_object[which(post_odds_object$is_annotate_reactive=="yes"),"symbol"]), size=2) 
    }
    if(mediator_type == "Phenotypes"){
      plot_reactive <- plot_reactive +
        geom_label_repel(data = post_odds_object[which(post_odds_object$is_annotate_reactive=="yes"),], mapping = aes(label=post_odds_object[which(post_odds_object$is_annotate_reactive=="yes"),"ID"]), size=2) 
    }
      
      # custom X axis:
      if (mediator_type == "RNA/Protein"){
        plot_reactive <- plot_reactive +
          scale_x_continuous( label = axisdf$order, breaks= axisdf$center )
      }
      if (mediator_type == "Phenotypes"){
        plot_reactive <- plot_reactive +
        scale_x_continuous( label = axisdf$chr, breaks= axisdf$center )
      }
      
    plot_reactive <-plot_reactive +
      #scale_y_continuous(expand = c(0, 0) ) +     
      # remove space between plot area and x axis
      ylim(0,1.25)+
      
      # Custom the theme:
      theme_bw() +
      theme( 
        legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
      )+
      
      #change the axis lines
      ggplot2::theme(axis.line = element_line(colour = "black"))
      
      #change the axis labels
      if (mediator_type == "RNA/Protein"){
        plot_reactive <- plot_reactive +
          xlab("Chromosomal position")
      }
      if (mediator_type == "Phenotypes"){
      if (selected_mediator$Object.names[1]=="Clinical"){
        plot_reactive <- plot_reactive + xlab("Clinical phenotype")
      }
      if (selected_mediator$Object.names[1]=="Ex vivo"){
        plot_reactive <- plot_reactive + xlab("Islet ex vivo phenotype")
      }
      if (selected_mediator$Object.names[1]=="Lipids"){
        plot_reactive <- plot_reactive + xlab("Lipid class")
      }
      if (selected_mediator$Object.names[1]=="Metabolites"){
        plot_reactive <- plot_reactive + xlab("Metabolite class")
      }
        if (selected_mediator$Object.names[1]=="Modules"){
          plot_reactive <- plot_reactive + xlab("Module type")
        }
        plot_reactive <- plot_reactive + theme (axis.text.x = element_text (angle = 45, hjust = 1))+
        theme(axis.text.x = element_text(size=10))
      }
    plot_reactive <- plot_reactive+
      ylab("Probability of reactive effects")+
      theme(axis.text=element_text(size=18),
            axis.title=element_text(size=20))+
      
      #add title
      labs(title = paste0(selected_QTL$lodcolumn[cntr2], " QTL chr",
                          selected_chr,":",selected_pos,
                          " mediated by ",
                          selected_mediator$Tissue.compartment[1],
                          " ", selected_mediator$Type[1]))+
      
      
      #add grob annotation
      annotation_custom(grob)+
      
      #add the legend
      labs(linetype =c("Threshold"))
    
    #show(plot_reactive)
    
    
    #load the diagrams for the relevant elements. 
    plot_complete <- cowplot::plot_grid(diagrams$complete,
                                        plot_complete,
                                        rel_widths = c(1,4)) 
    #save these plots as 
    plot_partial <- cowplot::plot_grid(diagrams$partial,
                                       plot_partial,
                                       rel_widths = c(1,4)) 
    plot_colocal <- cowplot::plot_grid(diagrams$colocal,
                                       plot_colocal,
                                       rel_widths = c(1,4)) 
    plot_reactive <- cowplot::plot_grid(diagrams$reactive,
                                        plot_reactive,
                                        rel_widths = c(1,4)) 
    plot_all <- list(plot_complete, plot_partial, 
                     plot_colocal, plot_reactive)
    names(plot_all)<-c("complete", "partial", "colocal","reactive")
    }
    
    if (ln_prior_c == "complete"){
      plot_all <- list(plot_complete, plot_partial, 
                       plot_colocal, diagrams$summary)
      names(plot_all)<-c("complete", "partial", "colocal","diagrams")
    }
    
    return(plot_all)
  }

#load the data=========================================================================
#K------------------------------------------------------------------------
  cat("...loading the kinship matrix...\n")
    #if the K is chosen by user  
  if (is.character(K)==TRUE){
  if (K == "choose"){
    K <- readRDS(file.choose())
  }
  
    #if the K is not chosen in the run
  if (K!="choose"){
     K <- readRDS(K)
  }
  }
#genoprobs----------------------------------------------------------------
  cat("...loading the genome probabilities (genoprobs) file...\n")
  if (is.character(genoprobs)==TRUE){
  if (genoprobs == "choose")  {
    cat("...choose the genome probabilities (genoprobs) file: \n")
    genoprobs <- readRDS(file.choose())
  }
  if (genoprobs != "choose")  {
      genoprobs <- readRDS(genoprobs)
  }
  }
  
#map----------------------------------------------------------------------
  cat("...loading the chromosomal map\n")
  if (is.character(map)==TRUE){
  if (map == "choose"){
    cat("...choose the map file\n")
    maps <- readRDS(file.choose())
  }
  if (map != "choose"){
    maps <- readRDS(map)
  }
  }
  if (is.character(map)!=TRUE){
    maps <- map
  }
  
#load the phenotypes file-------------------------------------------------
  cat("...loading the phenotypes file...\n")
  if (is.character(pheno)==TRUE){
  if (pheno == "choose"){
    pheno <- readr::read_csv(file.choose())
    clnms <- colnames(pheno)
    pheno <- data.frame(pheno)
    colnames(pheno)<-clnms
    rownames(pheno)<-pheno$Mouse
  }
  if (pheno != "choose"){
    pheno <- readr::read_csv(pheno)
    clnms <- colnames(pheno)
    pheno <- data.frame(pheno)
    colnames(pheno)<-clnms
    rownames(pheno)<-pheno$Mouse
  }
  }
  
#load the mediator file---------------------------------------------------
  cat("...load the mediators...\n")
  if (is.character(mediator)==TRUE){
  if (mediator == "choose"){
    mediators <- readRDS(file.choose())
  }
  if (mediator != "choose"){
    mediators <- readRDS(mediator)
  }
  }
  if (is.character(mediator)!=TRUE){
    mediators <- mediator
  }
  #print(head(mediators))
  
#load the covariates------------------------------------------------------
  cat("...loading the covariates for the outcome variable...\n")
  if (is.character(covar)==TRUE){
  if (covar == "choose"){
    covar <- data.frame(read.csv(file.choose()), stringsAsFactors = FALSE)
    rownames(covar)<-covar$X
    covar <- covar[-1]
  }
  if (covar != "choose"){
    covar <- data.frame(read.csv(covar), stringsAsFactors = FALSE)
    rownames(covar)<-covar$X
    covar <- covar[-1]
  }
  }
  print(head(covar))
  
  if (is.null(covar_pheno)!=TRUE){
    cat("...loading the covariates for the mediator variable...\n")
    covar_pheno <- covar
    if (is.character(covar_pheno)==TRUE){
    if (covar_mediator == "choose"){
      covar_mediator <- data.frame(read.csv(file.choose()), stringsAsFactors = FALSE)
      rownames(covar_mediator)<-covar$X
      covar_mediator <- covar_mediator[-1]
    }
    if (covar_mediator != "choose"){
      covar_mediator <- data.frame(read.csv(covar_mediator), stringsAsFactors = FALSE)
      rownames(covar_mediator)<-covar_mediator$X
      covar_mediator <- covar_mediator[-1]
    }
    }
  }
  
  #load the markers---------------------------------------------------------
  cat("...loading the marker file\n")
  
  if (is.character(markers)==TRUE){
  #if the markers are chosen by user
  if (markers == "choose"){
    cat("choose the marker file: \n")
    if (DO_set == "first"){
      mrkrs <- data.frame(read.csv(file.choose()), stringsAsFactors = FALSE)
    }
    if (DO_set == "second"){
      mrkrs <- readRDS(file.choose())
    }
  }
  
  #if the path for the markers file is provided
  if (markers != "choose"){
    if (is.character(markers)==TRUE){
      if (DO_set == "first"){
        mrkrs <- data.frame(read.csv(markers), stringsAsFactors = FALSE)
      }
      if (DO_set == "second"){
        mrkrs <- readRDS(markers)
      }
    }
      mrkrs <- markers
  }
  }
  
  if (is.character(markers)!=TRUE){
    mrkrs <- markers
  }
  
  
  #load the chromosomal breaks----------------------------------------------
  cat("...loading the chromosomal breaks file\n")
  if (is.character(chr_breaks)==TRUE){
  #if the breaks are chosen by user
  if (chr_breaks == "choose"){
    cat("choose the chromosomal breaks file: \n")
    chr_breaks <- data.frame(read.csv(file.choose()), stringsAsFactors = FALSE)
  }
  
  #if the path for the chr_breaks file is provided
  if (chr_breaks != "choose"){
    chr_breaks <- data.frame(read.csv(chr_breaks), stringsAsFactors = FALSE)
  }  
  }

  #load the covariates------------------------------------------------------
  if (is.null(weights_all)!=TRUE){
    cat("...loading the weights for the outcome variable...\n")
    if (weights_all == "choose"){
      weights_all <- data.frame(read.csv(file.choose()), stringsAsFactors = FALSE)
      rownames(weights_all)<-weights_all$X
      weights_all <- weights_all[-1]
    }
    if (weights_all != "choose"){
      weights_all <- data.frame(read.csv(weights_all), stringsAsFactors = FALSE)
      rownames(weights_all)<-weights_all$X
      weights_all <- weights_all[-1]
    }
    if(weights_pheno == TRUE){
      cat("...loading the weights for the mediator variable...\n")
      weights_pheno <- weights_all
      if (weights_mediator == "choose"){
        weights_mediator <- data.frame(read.csv(file.choose()), stringsAsFactors = FALSE)
        rownames(weights_mediator)<-weights_mediator$X
        weights_mediator <- weights_mediator[-1]
      }
      if (weights_mediator != "choose"){
        weights_mediator <- data.frame(read.csv(weights_mediator), stringsAsFactors = FALSE)
        rownames(weights_mediator)<-weights_mediator$X
        weights_mediator <- weights_mediator[-1]
      }
    }
    
  }
    
  
#load the overall list of variables to test-------------------------------------
  cat("...loading the QTL list...\n")
  if (is.character(QTL_list)==TRUE){
  if (QTL_list == "choose"){
    QTL_lists <- data.frame(read.csv(file.choose()), stringsAsFactors = FALSE)
  }
  if (QTL_list != "choose"){
      QTL_lists <- data.frame(read.csv(QTL_list), stringsAsFactors = FALSE)
  }
  } 
  if (is.character(QTL_list)!=TRUE){
    QTL_lists <- QTL_list
  }
  QTL_lists <- unique(QTL_lists)
  print(head(QTL_lists))
  
#load the subset of variables to test-------------------------------------
  cat("...loading the file with the variables to test...\n")
  if (is.character(test_these)==TRUE){
  if (test_these == "choose"){
    test_these <- data.frame(read.csv(file.choose()), stringsAsFactors = FALSE)
    colnames(test_these)<-"Phenotype"
  }
  if (test_these != "choose"){
    test_these <- data.frame(read.csv(test_these), stringsAsFactors = FALSE)
    colnames(test_these)<-"Phenotype"
  }
  }
  
  #make sure the phenotype you want to mediate has a QTL in your QTL list
  print(unique(test_these$Phenotype %in% QTL_lists$lodcolumn))
  if (length(unique(test_these$Phenotype %in% QTL_lists$lodcolumn))>1){
    return("ERROR: Your desired phenotype has no QTL in your QTL list.")
    break
  }
  if (length(unique(test_these$Phenotype %in% QTL_lists$lodcolumn))==1){
    if (unique(test_these$Phenotype %in% QTL_lists$lodcolumn)=="FALSE"){
    return("ERROR: Your desired phenotype has no QTL in your QTL list.")
    break
    }
  }
  
#load the set of objects in the mediator object to look at----------------
  cat("...loading the selection of mediator objects to test...\n")
  if (is.character(mediator_list)==TRUE){
  if (mediator_list == "choose"){
    mediator_lists <- data.frame(read.csv(file.choose()), stringsAsFactors = FALSE)
  }
  if (mediator_list != "choose"){
       mediator_lists <- data.frame(read.csv(mediator_list), stringsAsFactors = FALSE)
  }
  }
  if (is.character(mediator_list)!=TRUE){
    mediator_lists <- mediator_list
  }
  
#load the model diagrams------------------------------------------------------------
  cat("...loading the diagrams...\n")
  if (is.character(diagrams)==TRUE){
  if (diagrams == "choose"){
    diagrams <- readRDS(file.choose())
    
  }
  if (diagrams != "choose"){
      diagrams <- readRDS(diagrams)   
  }
  }
#cc colors----------------------------------------------------------------
  if (cc_colors == "old"){
    cc_colors <- qtl2::CCcolors
  }
  if (cc_colors == "new"){
    #set colors
    newClrs <- qtl2::CCcolors
    newClrs[1]<-"#000000"
    newClrs[2]<-"#96989A"
    newClrs[3]<-"#E69F00"
    newClrs[4]<-"#0072B2"
    newClrs[5]<-"#619BFF"
    newClrs[6]<-"#009E73"
    newClrs[7]<-"#D55E00"
    newClrs[8]<-"#CC79A7"
    cc_colors <- newClrs
  }

#process the data========================================================================
 
#set repeat loops---------------------------------------------------------
  
  #repeat loop for determining traits of interest-------------------------
  #set counters
  cntr <- 0
  
  #begin repeat loop------------------------------------------------------
  repeat{
    #increment counter
    cntr <- cntr+1
    #grab the next phenotype
    selected_QTL <- unique(QTL_lists[which(QTL_lists$lodcolumn == test_these$Phenotype[cntr]),])
    #set the next QTL in the list-----------------------------------------
    #set counter
    cntr2 <-0
    #begin second repeat loop
    repeat{
      #increment second counter
      cntr2 <- cntr2 + 1
      #get the next peak chromosome & position
      selected_chr <- selected_QTL$chr[cntr2]
      selected_pos <- selected_QTL$pos[cntr2]
      print(paste0("Now on ", 
                   selected_QTL$lodcolumn[cntr2],
                   ", chromosome ",
                   selected_chr, " at ",
                   selected_pos, "Mbp"),
            quote = FALSE)
      
      #set 3rd counter
      cntr3 <- 0
      
      #set the repeat loop for the different mediators
      repeat{
        # increment 3rd counter
        cntr3 <- cntr3+1
        # find the next object
        selected_mediator <- mediator_lists[cntr3,]
        print(paste0("now using ", selected_mediator$Tissue.compartment[1],
              " ", selected_mediator$Type[1]), quote = FALSE)
        mediator_selected <- mediators[[selected_mediator$Major_object[1]]][[selected_mediator$Object.names[1]]]
        annotation_selected <- mediators[[selected_mediator$Major_info_object[1]]][[selected_mediator$Info_object[1]]]
        
        
        # find the common mice
        cat("...finding common mice...\n")
        common_mice <- Reduce(intersect, list(rownames(mediator_selected),
                                              rownames(pheno),
                                              rownames(genoprobs[[1]])),
                                              rownames(covar)
                              )
        print(paste0("there are ",length(common_mice)," common mice"), quote = FALSE)
        
        #make temporary genoprobs
        temporary_genoprobs <- pull_genoprobpos(genoprobs, maps, 
                                                chr = selected_chr, 
                                                pos = selected_pos)
                                                
        # set for common mice
        temporary_genoprobs <- temporary_genoprobs[common_mice,]
        pheno_subset <- pheno[common_mice,]
        mediator_selected <- mediator_selected[common_mice,]
        #print(head(mediator_selected))
        # rank z the mediators & phenotype:
        if (rankz_set == TRUE){
        pheno_subset[,colnames(pheno[selected_QTL$lodcolumn[cntr2]])]<-rankz(pheno_subset[,colnames(pheno[selected_QTL$lodcolumn[cntr2]])])
        }
        
        print(head(pheno_subset[,colnames(pheno[selected_QTL$lodcolumn[cntr2]])]))
        print(head(annotation_selected))
        
        covar_subset <- covar[common_mice,]
        if (mediator_type == "RNA/Protein"){
            if (DO_set == "first"){
        colnames(annotation_selected)[5]<-"pos"
        colnames(annotation_selected)[4]<-"chr"
            }
            if (DO_set == "second"){
              if (RNA_type == "genes"){
        annotation_selected$pos <- annotation_selected$start
        annotation_selected$pos <- annotation_selected$pos/(10^6)
        print(annotation_selected)
              }
            if (RNA_type == "isoforms"){
              annotation_selected$pos <- annotation_selected$Gene.start..bp.
              annotation_selected$pos <- annotation_selected$pos/(10^6)
            }
        #colnames(annotation_selected)[8]<-"chr" 
            }
            if (usewindow == TRUE){
              annotation_selected <- subset(annotation_selected, chr == selected_chr)
              #print(head(annotation_selected))
              annotation_selected <- subset(annotation_selected, pos >= (selected_pos-window))
              #print(head(annotation_selected))
              annotation_selected <- subset(annotation_selected, pos <= (selected_pos+window))
              #print(head(annotation_selected))
              clnms_mediator <- data.frame(colnames(mediator_selected))
              #print(head(clnms_mediator))
              colnames(clnms_mediator)<-"ID"
              #for (i in 1:nrow(clnms_mediator)){
              #  clnms_mediator$gene.id[i] <- str_split(clnms_mediator$ID[i], pattern="[.]")[[1]][1] 
              #}
              #print(head(clnms_mediator))
              if (selected_mediator$Type[1] == "Proteins") {
              clnms_mediator <- subset(clnms_mediator, ID %in% annotation_selected$protein.id)
              }
              if (selected_mediator$Type[1] == "RNA"){
                if (RNA_type == "genes"){
                  clnms_mediator <- subset(clnms_mediator, ID %in% annotation_selected$gene.id)
                }
                if (RNA_type == "isoforms"){
                  clnms_mediator <- subset(clnms_mediator, ID %in% annotation_selected$transcript.id)
                }
              }
              #print(head(mediator_selected))
              print(head(clnms_mediator))
              mediator_selected <- mediator_selected[,clnms_mediator$ID]
              #print(head(mediator_selected))
              #print(head(annotation_selected))
            }
        }
        if (mediator_type == "Phenotypes"){
          colnames(annotation_selected)[which(colnames(annotation_selected)=="Class")]<-"chr"
          colnames(annotation_selected)[which(colnames(annotation_selected)=="Position")]<-"pos"
        }
        if (rankz_set == TRUE){
          for (i in 1:ncol(mediator_selected)){
            mediator_selected[,i]<-rankz(mediator_selected[,i])
          }
        }
        
        if (rescan == TRUE){
        #revised QTL plot for the subset
        cat("...performing revised QTL scan...\n")
        summary(K)
        qtl.temp <- scan1(genoprobs = genoprobs, 
                          pheno = pheno_subset[,colnames(pheno[selected_QTL$lodcolumn[cntr2]]), drop = FALSE], 
                          kinship = K, addcovar = covar_subset, cores= ncors)
      
        if (save_rescan == TRUE){
          cat("...saving revised QTL scan...\n")
          saveRDS(qtl.temp, 
                    file=paste0(selected_QTL$lodcolumn[cntr2], "_QTL_chr",
                                selected_chr,"_",selected_pos,
                                "_rescanned_using_",
                                selected_mediator$Tissue.compartment[1],
                                "_", selected_mediator$Type[1], "common_mice.rds"))
        }
        if (save_rescan =="choose"){
          flname <- base::readline(prompt="What's the file name for the qtl rescan object (including the .rds)? ")
          saveRDS(qtl.temp, file=flname)
        }
        
        #determine if the peak is still there
        cat("...determining if peak is still there...\n")
        if (!("marker.id" %in% colnames(selected_QTL))==TRUE){
          rowname_check <- paste0(as.numeric(selected_chr),"_",(as.numeric(selected_pos)*(10^6)))
        }
        if(("marker.id" %in% colnames(selected_QTL))==TRUE){
          rowname_check <- selected_QTL$marker.id[cntr2]
        }
        
        cat(paste0("...LOD threshold is set to ",LOD_thr,"...\n"))
        print(paste0("the LOD at your selected position is ",qtl.temp[rowname_check,1]), quote=FALSE)
        if (as.numeric(qtl.temp[rowname_check,1])<LOD_thr){
          LOD_pass <- FALSE
          cat(paste0("...your peak is now less than ", LOD_thr, " so further analysis for this subset will be discontinued...\n"))
          post_odds <- data.frame(matrix(nrow=1, ncol=1))
          colnames(post_odds)<-"Peak below threshold"
          post_odds[1,1] <- paste0("...your peak is now less than ", LOD_thr, " so further analysis for this subset will be discontinued...\n")
          plot_QTL <- QTL_plot(qtl.temp)
          plot_LOD <- NULL
          plot_probs <- list(NULL,NULL, NULL)
          names(plot_probs) <- c("complete","partial", "colocal")
          output <- list(post_odds, plot_QTL, plot_LOD, plot_probs)
          names(output) <- c("tabular_results", "QTL_plot","LOD_drop_plot","bmediatR_plots")
          return(output)
        }
        if (as.numeric(qtl.temp[rowname_check,1])>=LOD_thr){
          cat("...your peak passes check, so this subset will be mediated...\n")
          LOD_pass <- TRUE
        }
        }
#=======================================================================================================  
        if (rescan == FALSE){
          LOD_pass <- TRUE 
          cat("... you chose not to rescan your dataset...")
        }
        if (LOD_pass == TRUE){
          
        
#look at the LOD drop=====================================================================================        
        if (DO_set == "second"){
            if (RNA_type == "genes"){
                annotation_selected$ID <- annotation_selected$gene.id
                #annotation_selected <- unique(annotation_selected[c(1,5:8,11)])
                #gene_ids <- data.frame(colnames(mediator_selected))
                #colnames(gene_ids)<-"ID"
                #gene_ids$gene.id <- gene_ids$ID
                #for (i in 1:nrow(gene_ids)){
                #    gene_ids$gene.id[i]<-str_split(gene_ids$ID[i], pattern="[.]")[[1]][1]
                #}
                #annotation_selected <- gene_ids[c("ID","gene.id")] %>% inner_join(., annotation_selected, by=c("gene.id"="gene.id"))
                #mediator_selected <- mediator_selected[,annotation_selected$ID]
            }
             if (RNA_type == "isoforms"){
                annotation_selected$ID <- annotation_selected$transcript.id
                #gene_ids <- data.frame(colnames(mediator_selected))
                #colnames(gene_ids)<-"ID"
                #gene_ids$transcript.id <- gene_ids$ID
                #for (i in 1:nrow(gene_ids)){
                #    gene_ids$transcript.id[i]<-str_split(gene_ids$ID[i], pattern="_")[[1]][1]
                #}
                #annotation_selected <- gene_ids %>% inner_join(., annotation_selected, by=c("transcript.id"="transcript.id"))
                #mediator_selected <- mediator_selected[,annotation_selected$ID]
            }
        }
#print(head(annotation_selected))
#print(head(mediator_selected))
        
          
        #run the scan
          cat("...calculating the LOD drop...\n")
          assign("annotation_selected_obj", annotation_selected, envir = .GlobalEnv)
          assign("mediator_selected_obj", mediator_selected, envir = .GlobalEnv)
          inter_mediation <- mediation.scan(target = as.matrix(pheno_subset[test_these$Phenotype[cntr]]),
                                            mediator = as.matrix(mediator_selected),
                                            annotation = annotation_selected,
                                            covar = as.matrix(covar_subset),
                                            method = "double-lod-diff",
                                            qtl.geno = temporary_genoprobs,
                                            verbose = TRUE
          )
          #reset the mito to "M" instead of "MT". This is necessary because
          #the plotting functions will give errors if it is left as "MT"
          inter_mediation$chr[which(inter_mediation$chr=="MT")]<-"M"
          
          #determine fraction drop from starting QTL
          inter_mediation$new_median <- median(inter_mediation$LOD)
          for (i in 1:nrow(inter_mediation)){
            inter_mediation$median_drop[i] <- (inter_mediation$new_median[i]-
                                                 inter_mediation$LOD[i])/inter_mediation$new_median[i]
            inter_mediation$original_drop[i] <- (selected_QTL$lod[cntr2]-
                                                   inter_mediation$LOD[i])/selected_QTL$lod[cntr2]
          }
          
          print(head(inter_mediation))
        
#potentially look at different bmediatR effects for specific loci          
          
          #set 4th counter
          cntr4 <-0
          
          if (show_fits == TRUE){
            cat("...calculating fits for individual components...\n")
          #start repeat for the different components of the selected mediator
            #increment 4th counter

            repeat{
            cntr4 <- cntr4+1
            if (cntr4 %% 1000 == 0){print(cntr4)}
            #begin mediation analysis
            
            target_effect <- fit1(genoprobs = temporary_genoprobs,
                                   pheno = pheno_subset[test_these$Phenotype[cntr]],
                                   K = K[[selected_chr]],
                                   blup = T)
            
            mediator_effect <- fit1(genoprobs = temporary_genoprobs,
                                    pheno = mediator_selected[cntr4],
                                    K = K[[selected_chr]],
                                    blup = F)
            
            title_info <- paste0("QTL effects for ", selected_QTL$lodcolumn[cntr2]," ",
                                 selected_chr,
                            ":",selected_pos," mediated by ", 
                            selected_mediator$Tissue.compartment[cntr3],
                            " ",selected_mediator$Type[cntr3],": ",
                            colnames(mediator_selected)[cntr4]
                            )
            
            qtl_effects <- cbind(target_effect$coef, target_effect$SE,
                                 mediator_effect$coef, mediator_effect$SE) 
            qtl_effects <- data.frame(qtl_effects)
            colnames(qtl_effects)<-c("target_effect", 
                                     "target_SE", 
                                     "mediator_effect", 
                                     "mediator_SE")
            qtl_effects <- qtl_effects[1:8,] 
            #plot the coefficients---------
            
            #use the plot function
            plt_effect <- effects_w_ci_plot(qtl_effects, cc_colors, title_info)
            show(plt_effect)
            
            done_withit <- base::readline(prompt="Done looking at mediation effects? Y/N: ")
            if (done_withit=="Y"){break}
            if (cntr4 >= ncol(mediator_selected)){break}
            }    
            }
            
#determine & graph the combined effects=============================================================            
            #bmediatR set---------------------------
              #from the help of bmediatR
              #y is the matrix designating the single outcome (phenotype)
                #that is measured by the QTL you're mediating
              #M is the matrix of mediators (e.g. RNA)
              #X is the genome probabilities for the QTL position
              #Z is the covariate matrix for both QTL and mediator
                #you can supply independent covariates for outcome (the QTL)
                #and mediators as Z_y and Z_M respectively
              #see notes below about the model options
              
          #run the bmediatR scan for all elements in the selected mediator
          cat("...now running Bayesian mediation analysis (bmediatR)...\n")
            bmediatR_mediation <- bmediatR(
              y = as.matrix(pheno_subset[test_these$Phenotype[cntr]]), 
              M = as.matrix(mediator_selected), 
              X = temporary_genoprobs,
              Z = as.matrix(covar_subset),
              Z_y = covar_pheno,
              Z_M = covar_mediator,
              w = weights_all,
              w_y = weights_pheno,
              w_M = weights_mediator,
              ln_prior_c = ln_prior_c
            )
            
##==========#also from the vignette on using the bmediatR and the paper describing it=================================##
            #The causal models can be described as directed acyclic graphs (DAGs)
            #There are three possible edges (a, b, c) and we define each to be either present
            #or absent using an indicator vector V = (Va, Vb, Vc), where for example V = (1, 0, 0) denotes
            #presence of a only:
            ##
            ##
            ##                        :--edge a---mediator--edge b---:
            ##                        |                              |          
            ##                        |                              |
            ##                      QTL-----------edge c--------phenotype
            ##
            ##
            # for our purposes, the ones of interest are 
            #
            # V=(1,0,1) "co-local" effects where the QTL alters both the phenotype and 
            # mediator, but for which the mediator and phenotype are independent
            # 
            # V=(1,1,1) "partial mediation" where the QTL alters the phenotype
            # and the mediator, and there is a causal link where the mediator also
            # alters the phenotype
            # 
            # V=(1,1,0) "complete mediation" where the QTL alters only the 
            # mediator, and the mediator then alters the phenotype
            # 
            # in some special cases, and this has a lot of caveats suggested in
            # the paper, there are what they term "reactive" effects where the phenotype
            # alters the theoretical mediator. Two specific instances of this are
            # potentially of interest
            #
            # V=(0,*,1) "complete reactive mediation" where the QTL alters the
            # phenotype and the phenotype then alters the theoretical mediator
            #
            # V=(1,*,1) "partial reactive mediation" where the QTL alters the
            # phenotype and both the QTL and phenotype alter the theoretical mediator
            #
            # Interpreting the "reactive" mediation effects should take into context
            # what would make these biologically plausible. 
            #
            # All other effects suggest either no phenotype-QTL link (Vc = 0),
            # no effect of phenotype on the theoretical mediator (Va=0),
            # or no effect of the theoretical mediator on the QTL (Vb = 0)
            # or some combination thereof. This is referred to collectively as
            # "no mediation effect". 
            #
            # These vector designations are important to keep clear in mind because they 
            # are used as column headings for the bmediatR object column headings
            # 
##==================================================================================================================##

#join the bmediatR obj's posterior odds data to the info about the mediators
post_odds <- data.frame(bmediatR_mediation$ln_post_odds)
print(head(post_odds))
if (mediator_type=="RNA/Protein"){
  if (selected_mediator$Type[1] == "RNA"){
    if (RNA_type == "genes"){
      if (DO_set == "first"){
post_odds$gene.id <- rownames(post_odds)
post_odds <- post_odds %>%
  left_join(annotation_selected, ., by=c("gene.id"="gene.id"))
      }
      if (DO_set == "second"){
        post_odds$ID <- rownames(post_odds)
        post_odds <- post_odds %>%
  left_join(annotation_selected, ., by=c("ID"="ID"))
        }
      }
  }
  if (RNA_type == "isoforms"){
   post_odds$ID <- rownames(post_odds)
post_odds <- post_odds %>%
  left_join(annotation_selected, ., by=c("ID"="ID"))
  
  }
  }
  if (selected_mediator$Type[1] == "Proteins"){
    post_odds$protein.id <- rownames(post_odds)
    post_odds <- post_odds %>%
      left_join(annotation_selected, ., by=c("protein.id"="protein.id"))
  }
  if (any(post_odds$chr == "MT")){
post_odds$chr[which(post_odds$chr == "MT")]<-"M"
}

if (mediator_type=="Phenotypes"){
  post_odds$ID <- rownames(post_odds)
  post_odds <- post_odds %>%
    left_join(annotation_selected, ., by=c("ID"="ID"))
  
}
#if determining mediator pvalues and not using reactive analysis
if (pval == TRUE){
  if (ln_prior_c != "reactive"){
    
    if (mediator_type=="RNA/Protein"){
    #determine approximate p-values
    approx_pvals <- get_approx_pval(
      bmediatR_object = bmediatR_mediation,
      model_type = pval_model_type,
      med_annot = annotation_selected,
      med_var = "gene.id"
    )
    
    colnames(approx_pvals)[which(colnames(approx_pvals)=="ln_post_odds")]<-pval_model_type
    }
    
    if (mediator_type=="Phenotypes"){
      #determine approximate p-values
      approx_pvals <- get_approx_pval(
        bmediatR_object = bmediatR_mediation,
        model_type = pval_model_type,
        med_annot = annotation_selected,
        med_var = "ID"
      )
      
      colnames(approx_pvals)[which(colnames(approx_pvals)=="ln_post_odds")]<-pval_model_type
    }
    
    #determine permuted p-values
    #this is impractical given the limits of my computer. It takes many minutes to run one permutation and
    #it will be impossibly long to accomplish for 23k transcripts in multiple tissues even for
    #one of the phenotypes/QTLs. This has to be done on the CHTC
    #note I originally went with the built-in function but it re-does the bmediatR scan and gives errors
    #permuted p-values----------------------------------------------------------------
    if (perm_p == TRUE){
    actual_po <- bmediatR_mediation
    model_type <- pval_model_type
    within_cov <- as.matrix(covar) %>% unique
    Z <- as.matrix(covar)
    original_name <- NULL
    perm_mat <- NULL
    num_perm <- 1000
    for (i in 1:nrow(within_cov)) {
      sub_Z <- Z[apply(Z, 1, function(x) paste(x, collapse = "")) == paste(within_cov[i,], collapse = ""),, drop = FALSE]
      sub_original_name <- rownames(sub_Z)
      sub_perm_mat <- matrix(NA, nrow = nrow(sub_Z), ncol = num_perm)
      for (j in 1:num_perm) {
        sub_perm_mat[,j] <- sample(sub_original_name)
      }
      original_name <- c(original_name, sub_original_name)
      perm_mat <- rbind(perm_mat, sub_perm_mat)
    }
    M <- as.matrix(mediator_selected)
    if (is.null(Z)) { Z <- matrix(1, nrow = nrow(as.matrix(mediator_selected))); rownames(Z) <- rownames(as.matrix(mediator_selected)) }
    rownames(perm_mat) <- original_name
    perm_mat <- perm_mat[rownames(M),]
    
    perm_pval <- rep(NA, ncol(M))
    perm_po_mat <- matrix(NA, nrow = num_perm, ncol = ncol(M))
    colnames(perm_po_mat) <- colnames(M)
    for(i in 1:num_perm) {
      ## Rename things
      perm_M <- M
      rownames(perm_M) <- as.character(perm_mat[,i])
      for (j in 1:ncol(M)) {
        perm_med <- bmediatR(y = as.matrix(pheno_subset[test_these$Phenotype[cntr]]), 
                             M = perm_M[,j,drop = FALSE], 
                             X = temporary_genoprobs, 
                             Z = Z, 
                             w = w,
                             ln_prior_c = ln_prior_c,
                             verbose = FALSE
                             )
        perm_po_mat[i, j] <- perm_med$ln_post_odds[,model_type]
      } 
      print(paste("Perm", i, "done out of", num_perm))
    }
    }
  }
}

if (plot_mediation == TRUE){
  #plot the QTL plot for the subset
  cat("...plotting revised QTL scan...\n")
  if (rescan == TRUE){
  plot_QTL <- QTL_plot(qtl.temp)
  }
  if (rescan == FALSE){
    plot_QTL <- NULL
  }
  #plot the LOD drop
  cat("...plotting LOD drop...\n")
  plot_LOD<-LOD_drop(inter_mediation, lod_drop_thr)
  #plot the probabilities
  cat("...making probability graphs...\n")
  plot_probs <- plot_bmediatR_odds(post_odds = post_odds, thr = probability_effect_thr)
  #combine the plots into a single R file
  if (combine_graphs == TRUE){
  cat("...combining graphs...\n")
  if (ln_prior_c == "reactive"){
    plot_all <- plot_grid(plot_QTL, plot_LOD, plot_probs$complete, plot_probs$partial, plot_probs$colocal, plot_probs$reactive, ncol = 2, nrow = 3)
  }
  if (ln_prior_c == "complete"){
    plot_all <- plot_grid(plot_QTL, plot_LOD, plot_probs$complete, plot_probs$partial, plot_probs$colocal, plot_probs$diagrams, ncol = 2, nrow = 3)
  }
  }
  if (save_plot_obj == TRUE){
  cat("...saving graph object...\n")
  saveRDS(plot_all, file=paste0(selected_QTL$lodcolumn[cntr2], "_QTL_chr",
                                selected_chr,"_",selected_pos,
                                "_mediated_vs_",
                                selected_mediator$Tissue.compartment[1],
                                "_", selected_mediator$Type[1], "plots.rds"))
  }
  
  #save the pdf file
  if (save_PDF_file == TRUE){
    cat("...saving graph pdf...\n")
    if (pdf_name_prompt == FALSE){
      pdf_file <- paste0(selected_QTL$lodcolumn[cntr2], "_QTL_chr",
                         selected_chr,"_",selected_pos,
                         "_mediated_vs_",
                         selected_mediator$Tissue.compartment[1],
                         "_", selected_mediator$Type[1], "_allplots.pdf")
    }
    if (pdf_name_prompt == TRUE){
      pdf_file <- base::readline(prompt="choose the pdf file name (inc the .pdf!): ")
    }
    pdf(file= pdf_file,
        width = 52, height = 18)
    show(plot_all)
    dev.off()}
  }

#combine the LOD drop & bmediatR objects===============================================
cat("...joining LOD drop and probabilities...\n")
if (mediator_type == "RNA/Protein"){
  if (RNA_type == "genes"){
  #colnames(inter_mediation)[1]<-"gene.id"
  #colnames(post_odds)[1]<-"gene.id"
  
post_odds <- post_odds %>%
  left_join(inter_mediation, by=c("gene.id"="gene.id"))
  }

  if (RNA_type == "isoforms"){
  post_odds <- post_odds %>%
  left_join(inter_mediation, by=c("ID"="ID"))
  }
}

if (mediator_type != "RNA/Protein"){
  post_odds <- post_odds %>%
  left_join(inter_mediation, by=c("ID"="ID"))
  colnames(post_odds)[1]<-"ID"
}
if (save_mediation == TRUE){
  cat("...saving mediations...\n")
  write.csv(post_odds, 
            file=paste0(selected_QTL$lodcolumn[cntr2], "_QTL_chr",
                        selected_chr,"_",selected_pos,
                        "_mediated_vs_",
                        selected_mediator$Tissue.compartment[1],
                        "_", selected_mediator$Type[1], ".csv"))
}
if (save_mediation =="choose"){
  flname <- base::readline(prompt="What's the file name for the mediation (including the .csv)? ")
  write.csv(post_odds, file=flname)
}
output <- list(post_odds, plot_QTL, plot_LOD, plot_probs, mediator_selected, covar_subset, annotation_selected, temporary_genoprobs, pheno_subset[test_these$Phenotype[cntr]])
names(output) <- c("tabular_results", "QTL_plot","LOD_drop_plot","bmediatR_plots", "selected_mediator_matrix", "subset_covariates","annotation_selected", "temporary_genoprobs", "pheno_subset")
return(output)
}            
#break the inner loop selecting the mediators
if (cntr3 >= nrow(mediator_lists)){
  break
}
}
#break the inner loop selecting the QTL
      if (cntr2 >= nrow(selected_QTL)){break}
    }
    
#break the outer loop selecting the phenotypes    
if (cntr >= nrow(test_these)){
  cat("...Finished mediating...\n")
  end_time <- Sys.time()
  print(paste0("Start time was ", start_time), quote=FALSE)
  print(paste0("End time was ", end_time), quote=FALSE)
      break
    }
  }

  
#end function----------------------------------------------------------------------------
  #======================================================================================
}