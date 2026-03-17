#' Attie_RNA_WGCNA_metabolites
#' 
#' @usage 
#' 
#' This is the  modified version for the metabolites data
#' it takes 3 main arguments which are the filenames for the 
#' relevant files and the number of cores to use
#' 
#' runs the default settings of the program, without showing the 
#' intermediate graphs and prompting for continuance
#' 
#' typing this into the command line:
#' Attie_WGCNA_metabolites(display_part = TRUE)
#' will run the function with default settings but will also display
#' intermediate steps and prompt the user for continuing
#' 
#' typing this into the command line:
#' Attie_WGCNA_metabolites(useblockwise= TRUE) uses the blockwiseModules 
#' command from the WGCNA package instead of the step-by-step
#' of this script. That will give slightly different results
#' for module assignments. I'm not sure yet what internal settings
#' cause this
#' 
#' 
#' @description 
#' The SCRIPT to end ALL SCRIPTS
#' Written for n00bs by n00bs
#' God help us all
#' 
#' notes regarding the installation:
#' This implements the basic WGCNA package
#' For use in analyzing mutliple source datas (e.g., RNAs from mutliple tissues)
#' 
#' Dependencies::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#' this was created in Windows 10 x86-64 build 19044
#' with Microsoft's R Open based on R 4.0.2
#' 
#' see https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/
#' and https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/
#' for information and tutorials on the WGCNA package
#' 
#' requires: rstudioapi, WGCNA, dplyr, tidyverse
#' rstudioapi is installed by Rstudio
#' WGCNA and its dependencies are installed using Bioconductor
#' dplyr and tidyverse are installed by the installed.packages command
#' in the past, I've had issues with the version control for some of the 
#' dependencies, depending on the version of R
#' The problems' sources aren't clear but I'd go through the missing
#' dependencies and see if you can install them before the needed package
#' 
#' Installing these packages requires installing Rtools, which for
#' R Open (2020-06-22) uses the version at https://cran.r-project.org/bin/windows/Rtools/rtools40.html
#' 
#' Notes regarding the linear algebra::::::::::::::::::::::::::::::::::::::::::::::::::
#' 
#' This packages uses a lot of linear algebra. The base version of R will default
#' to the system's basic linear algebra subcommands (BLAS)
#' The default BLAS for most computers is unreasonable slow. In some cases for 
#' large datasets particularly may take long enough to cause timeouts and
#' other errors
#' 
#' it is HIGHLY RECOMMENDED to upgrade the BLAs to an alternate, faster 
#' version. 
#' Guides on how to do this can be found online. An example which
#' benchmarks the various versions of the BLAS alternates available
#' at the time this function was made can be found at:
#' https://csantill.github.io/RPerformanceWBLAS/
#' 
#' Based on this, I chose to use R Open since I created this script
#' on Windows and R Open's Intel MKL benchmarks performed extremely
#' well
#' 
#' In this case the upgraded BLAS improved the speed for the computations
#' from hours to < 10 minutes for datasets of size 20,000+ probes/tissue type across
#' 371 tissue samples in multiple tissue types
#' 
#' Note 2: if switching from a newer version of R to R Open, which is 
#' optimized for R 4.0.2, on some but not all computers we experienced
#' R being unable to see the packages of the previous install of R
#' I did not have this occur on my laptop, but it did occur on one 
#' of our desktops, after which I had to reinstall the packages
#' this can create some unusual issues since the installation for
#' the newer package versions will create dependency errors during the
#' installation process as these newer versions look for dependencies 
#' not natively installed in 4.0.2 and for which you'll have to look
#' for legacy versions
#' 
#' As a consequence, I had a difficult time looking for the correct
#' dependencies and getting them installed. In particular I had issues
#' with rlang, Rtools, and GO.db for the WGCNA package's dependencies
#' 
#' #the choice of network type is important and one should read the WGCNA 
#' tutorials
#' 
#' #notes on defaults:::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#' 
#' The default analysis type is not to use the blockwisemodules function
#' The default network type is 'signed'
#' The default beta is 12
#' The default module size is 10
#' the default TOM type is 'signed'
#' the default clustering method is 'average'
#' the default distance cutoff is 0.25
#' the default kME threshold for pruning is 0.365
#' 
#' @author 
#' Chris Emfinger, PhD. Programatically challenged
#' 
#' @details
#' This assumes a specific type of format for the loaded data
#' it must be in CSV format
#' samples must columns
#' probes must be rows
#' 1st 2 columns must be the probe ID (metabolite) and
#' the second must be some other variable. Some metabolites appear
#' more than once so there must be something else to assign names
#' to them
#' 
#' @param default is a TRUE/FALSE boolean
#' Its default is TRUE
#' if TRUE, the beta (soft power), and other options are 
#' set to the Attie lab defaults. 
#' 
#' if FALSE, the script will ask for some specific command
#' parameters
#' 
#' @param display_part is a TRUE/FALSE boolean
#' its default is FALSE
#' if TRUE, it will display certain aspects of the
#' process, the clustering tree, etc
#' 
#' @param useblockwise is a TRUE/FALSE boolean
#' this designates whether the script will use its step-by-step
#' analysis method or the blockwiseModules command
#' default is FALSE
#' 
#' @param beta sets the softpower. default is 12. Alternately, "soft_power_check"
#' will allow the user to examine the scale-free topology and determine the 
#' beta value themselves
#' 
#' @param probe_path designates the path for the file containing
#' the probes (the metabolites)
#' if NULL, it prompts the user to load the files
#' 
#' @param ncores tells the program how many cores to use. default is 1
#' 
#' @param netType tells the program the kind of adjacency matrix type
#' default is "signed". Can be 'unsigned', 'signed', 'signed hybrid',
#' see the adjacency function information for details
#' 
#' @param modSize tells the program the minimum number of transcripts
#' per module. Default is 20 transcripts
#' 
#' @param TOMtype sets the type of TOM matrix calculation
#' options are:
#' 'unsigned', 'none', 'unsigned 2', 'signed Nowick',
#' 'signed 2', 'signed Nowick 2'. See the TOMsimilarityFromExpr
#' function for details
#' 
#' @param thrsh_kME sets the threshold of kME for 'pruning' poorly
#' correlated transcripts; Default is 0.365
#' 
#' @param fst_clust_method sets the method for clustering
#' default is "average"
#' options are 'single', 'complete', 'average', 'mcquitty',
#' 'ward.D', 'ward.D2', 'centroid', 'median',
#' see the fastcluster::hclust function information for details.
#' 
#' @param MEDissThres sets the dissimilarity threshold for
#' merging modules. default is 0.25
#' 
#' @param cr_use sets the cor() function use option for how it handles
#' missing data. Options are "all.obs" and "pairwise.complete.obs"
#' 
#' @param wd_path sets the writing path. if FALSE, it will prompt the user for 
#' 
#' 

#start==============================================================
#start the function:
Attie_RNA_WGCNA_Metabolites <- function(probe_path = NULL,
                        wd_path = TRUE,
                        ncores = 4,
                        beta =  "soft_power_check",
                        netType = "signed",
                        modSize = 20,
                        TOMtype = "signed",
                        fst_clust_method = "average",
                        thrsh_kME = 0.365,
                        MEDissThres = 0.25,
                        eigenfile = paste0("Module_eigengenes_md",modSize,".csv"),
                        allmodules = paste0("Module_membership_md",modSize,".csv"),
                        default = TRUE, 
                        display_part = FALSE,
                        cr_use = "pairwise.complete.obs",
                        useblockwise = FALSE){

# basic checks======================================================
#===================================================================

# check & load  required libraries==================================================
library("rstudioapi")
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
if (!require("WGCNA", quietly = TRUE)){
  pckg_missing <- c(pckg_missing, " WGCNA")
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
if (numb_missing == 0){
  require("rstudioapi")
  require("WGCNA")
  require("dplyr")
  require("tidyverse")
}
if (numb_missing >= 1){
  print(paste0("ERROR: You need ", numb_missing, " packages."), quote = FALSE)
  print(pckg_missing, quote = FALSE)
  break
}

#set the wd===============================================================
if (wd_path == FALSE){
print("choose the directory where you'd like to save the files", quote = FALSE)
setwd(selectDirectory())
}
if (wd_path == TRUE){
  print("choose the directory where you'd like to save the files", quote = FALSE)
  setwd(getwd())
}

#enable multithreading & set directory path=======================================================
enableWGCNAThreads()


#==========================================================================
#load & process data for analysis==========================================
#set options; note the tutorial says this cannot be omitted
options(stringsAsFactors = FALSE)

#load data
print("load the data", quote = FALSE)
if (is.null(probe_path) == TRUE){
  probe_path <- file.choose()
}
combined_probes <- data.frame(read.csv(probe_path), 
                              stringsAsFactors = FALSE)

#convert data to a matrix
print("converting uploaded data for processing", quote = FALSE)
clnms <- colnames(combined_probes)
clnms[1]<-"Accession"
clnms[2]<-"Compound"
colnames(combined_probes)<-clnms
combined_probes$Compound<-combined_probes$Accession
combined_probes$Accession <- paste0(combined_probes$Compound,"_", rownames(combined_probes))
rownames(combined_probes) <- combined_probes$Accession
combined_probes <- combined_probes[-c(1:2)]
combined_probes <- t(combined_probes)

#load the gene ID info
#print("load the file with the gene info", quote = FALSE)
#gene_info <- data.frame(read.csv(info_path), 
#                        stringsAsFactors = FALSE)


#adjust the soft power threshold===================================================
if (beta == "soft_power_check"){
  cat("now you'll determine the soft power threshold...\n")
  
  # Choose a set of soft-thresholding powers
  powers = c(seq(4,10,by=1), seq(12,20, by=2));
  
  #reformat into the set object
  
  # Initialize a list to hold the results of scale-free analysis
  nSets = 1
  powerTables = vector(mode = "list", length = nSets);
  
  # Call the network topology analysis function for each set in turn
  for (set in 1:nSets){
    powerTables[[set]] = list(data = pickSoftThreshold(combined_probes, 
                                                       powerVector=powers,
                                                       blockSize = NULL,
                                                       networkType = netType,
                                                       verbose = 2)[[2]])
  }
  collectGarbage()
  
  # Plot the results:
  colors = c("black", "red")
  # Will plot these columns of the returned scale free analysis tables
  plotCols = c(2,5,6,7)
  
  colNames = c("Scale Free Topology Model Fit", "Mean connectivity", "Median connectivity",
               "Max connectivity")
  
  # Get the minima and maxima of the plotted points
  ylim = matrix(NA, nrow = 2, ncol = 4)
  for (set in 1:nSets){
    for (col in 1:length(plotCols)){
      ylim[1, col] = min(ylim[1, col], 
                         powerTables[[set]]$data[, plotCols[col]], 
                         na.rm = TRUE);
      ylim[2, col] = max(ylim[2, col], 
                         powerTables[[set]]$data[, plotCols[col]], 
                         na.rm = TRUE);
    }
  }
  
  # Plot the quantities in the chosen columns vs. the soft thresholding power
  sizeGrWindow(8, 6)
  par(mfcol = c(2,2))
  par(mar = c(4.2, 4.2 , 2.2, 0.5))
  cex1 = 0.7
  
  for (col in 1:length(plotCols)) for (set in 1:nSets){
    if (set==1){
      plot(powerTables[[set]]$data[,1],
           -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
           xlab="Soft Threshold (power)",
           ylab=colNames[col],type="n",
           ylim = ylim[, col],
           main = colNames[col])
      addGrid()
    }
    if (col==1){
      text(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
           labels=powers,cex=cex1,col=colors[set]);
    }
    else
      text(powerTables[[set]]$data[,1], powerTables[[set]]$data[,plotCols[col]],
           labels=powers,cex=cex1,col=colors[set])
  }
  
  #select the power you'd like
  cat("now you'll select the soft power (beta)\n")
  beta <- base::readline(prompt="What softpower would you like? ")
  beta <- as.numeric(beta)
  
}




#cluster & extract MEs & ME IDs==============================================
#============================================================================

#if using the blockwiseModules command=========================================================
if (useblockwise == TRUE){
  print("You've chosen to use the single-step blockwiseModules command.", quote = FALSE)
  maxblock <- base::readline(prompt="enter the max block size: ")
  maxblock <- as.numeric(maxblock)
  
  bwnet <- blockwiseModules(combined_probes, maxBlockSize = maxblock,
                            power = beta, networkType = netType, 
                            TOMType = TOMtype, minModuleSize = modSize,
                            deepsplit=2,
                            reassignThreshold = 0, mergeCutHeight = MEDissThres,
                            numericLabels = FALSE,
                            loadTOM = FALSE,
                            saveTOMs = FALSE,
                            saveTOMFileBase = "combo_RNAs_TOMs_RZ_blockwise",
                            verbose = 3
                            )
  
  #extract the ME assigments
  print("getting the ME assignments", quote = FALSE)
  MEs<-data.frame(bwnet$colors, stringsAsFactors = FALSE)
  colnames(MEs)<-"module"
  tissue_origin <- MEs
  tissue_origin$module <- rownames(tissue_origin)
  colnames(tissue_origin)<-"tissue_origin"
  for (i in 1:nrow(tissue_origin)){
    tissue_origin$tissue_origin[i]<-str_split(tissue_origin$tissue_origin[i], pattern="_")[[1]][1]
  }
  colnames(tissue_origin)<-"Compound"
  MEs <- cbind(MEs, tissue_origin)
  
  #add the kME values for the pruned modules to the file
  print("adding the kME values", quote = FALSE)
  kMEs <- signedKME(combined_probes, bwnet$MEs)
  kval <- cbind(MEs[1],MEs[1])
  colnames(kval)<-c("kME_module", "module")
  for (i in 1:nrow(kval)){
    kval$kME_module[i] <- kMEs[rownames(kval)[i], paste0("kME",kval$module[i])]
  }
  kval<-kval[-2]
  ME_data <- cbind(MEs, kval)
  MEs3 <- subset(ME_data, module != "grey")
  
  #write all of the ME assignments
  print("saving all ME assignments", quote = FALSE)
  write.csv(MEs3, file = allmodules, row.names = FALSE)
  
  #extract the eigengene values
  print("extracting eigengene values", quote = FALSE)
  eigen <- data.frame(bwnet$MEs, stringsAsFactors = FALSE)
  print("saving eigengene values")
  write.csv(eigen, file = eigenfile)
  
  print("saving WGCNA objects", quote = FALSE)
  save(list=c("bwnet", "MEs", "kMEs", "eigen"), file="blockwise_WGCNA_objects_ms10.RData")
  
  
  
  # Plot the dendrogram and the module colors underneath for block 1
  print("plotting the dendrogram and modules", quote = FALSE)
  
  pdf(file = "Dendrogram_with_colors_assigned.pdf", width = 12, height = 4)
  plotDendroAndColors(bwnet$dendrograms[[1]], bwnet$colors,
                      "Module colors", main = "Gene dendrogram and module colors in block 1",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = FALSE, guideHang = 0.05,
                      rowText = bwnet$colors,
                      rowTextAlignment = "center",addTextGuide =TRUE,
                      rowTextIgnore = "grey"
                      )
  #base::readline(prompt="Press enter to end the function.")
  dev.off()
}

#if not using the blockwiseModules command===================================================
if (useblockwise == FALSE){
  print("you're using the step-by-step method", quote = FALSE)
  
  #construct the adjacency matrix
  print("constructing the signed adjacency matrix", quote = FALSE)
  adjacency <- adjacency(combined_probes, 
                        type = netType,
                        power = beta
                        )
  
  #calculate topological overlap
  print("calculating topological overlap", quote = FALSE)
  probeTOM <- TOMsimilarity(adjacency, 
                            TOMType = TOMtype,
                            verbose = 3)
  print("removing adjacency for memory", quote = FALSE)
  rm(adjacency) 
  print("clearing unused R memory", quote = FALSE)
  gc()
  
  #calculate the dissimilarity
  print("calculating topological dissimilarity", quote = FALSE)
  dissTOM <- 1-probeTOM
  print("clearing the TOM for memory", quote = FALSE)
  rm(probeTOM)
  print("clearing unused R memory", quote = FALSE)
  gc()
  
  #cluster
  print(paste0("clustering by ",fst_clust_method), quote = FALSE)
  geneTree <- fastcluster::hclust(as.dist(dissTOM), 
                                  method = fst_clust_method)
  
  if (display_part == TRUE){
    print("displaying the clustering tree")
    plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
       labels = FALSE, hang = 0.04)
  base::readline(prompt="Press ENTER to continue or hit Esc to interrupt")
  }
  
  #cut the tree
  print("cutting the tree and defining modules", quote = FALSE)
  dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                              deepSplit = 2,
                              pamRespectsDendro = TRUE,
                              minClusterSize = modSize,
                              verbose=3)
 
  dynamicColors <- labels2colors(dynamicMods)
  print(paste0("There are ", length(unique(dynamicColors)), " modules"))
  table(dynamicColors)
  
  if (display_part == TRUE){
    print("displaying the initial modules", quote = FALSE)
    plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                        dendroLabels = FALSE, hang = 0.03,
                        addGuide = FALSE, guideHang = 0.05,
                        main = "Gene dendrogram and module colors")
    base::readline(prompt="Press ENTER to continue or hit Esc to interrupt")
  }
  
  ##Calculate eigengenes=================================================
  print("calculating initial eigengenes", quote = FALSE)
  MEList = moduleEigengenes(combined_probes, 
                   colors= dynamicColors, 
                   impute = FALSE, 
                   nPC = 1, 
                   align = "along average", 
                   excludeGrey = FALSE, 
                   grey = if (is.numeric(colors)) 0 else "grey",
                   subHubs = TRUE,
                   trapErrors = FALSE, 
                   returnValidOnly = FALSE, 
                   softPower = beta,
                   scale = TRUE,
                   verbose = 3, indent = 0)
  
    
  MEs = MEList$eigengenes
  
  # Calculate dissimilarity of module eigengenes
  print("calculating dissimilarity of initial eigengenes", quote = FALSE)
  MEDiss = 1-cor(MEs, use = cr_use);
  # Cluster module eigengenes
  METree = fastcluster::hclust(as.dist(MEDiss), method = fst_clust_method)
  
  if (display_part == TRUE){
    print("displaying the initial module eigengene clustering", quote = FALSE)
    plot(METree, main = "Clustering of module eigengenes",
         xlab = "", sub = "")
    base::readline(prompt="Press ENTER to continue or hit Esc to interrupt")
  }
  
  if (display_part == TRUE){
    print("the red line is the threshold", quote = FALSE)
    # Plot the cut line into the dendrogram
    abline(h=MEDissThres, col = "red")
    base::readline(prompt="Press ENTER to continue or hit Esc to interrupt")
  }
  
  
  # Call an automatic merging function
  print("merging the modules", quote = FALSE)
  merge = mergeCloseModules(combined_probes, 
                            dynamicColors,
                            unassdColor="grey",
                            useAbs = FALSE,
                            relabel = TRUE,
                            getNewMEs = TRUE,
                            getNewUnassdME = TRUE,
                            iterate = TRUE,
                            cutHeight = MEDissThres, verbose = 3)
  # The merged module colors
  mergedColors = merge$colors;
  # Eigengenes of the new merged modules:
  mergedMEs = merge$newMEs
  
  if (display_part == TRUE){
    print("displaying the new module eigengene clustering", quote = FALSE)
    plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                        c("Dynamic Tree Cut", "Merged dynamic"),
                        dendroLabels = FALSE, hang = 0.03,
                        addGuide = FALSE, guideHang = 0.05)
    base::readline(prompt="Press ENTER to continue or hit Esc to interrupt")
  }
  
  # determine kMEs
  print("Determining the kME values", quote = FALSE)
  kMEs <- signedKME(combined_probes, mergedMEs)
  
  #prune the modules
  print(paste0("pruning genes to those with any abs kME >  ", thrsh_kME), quote = FALSE)
  gn_list <- colnames(combined_probes)
  MEcolors <- merge$colors
  MEcolors <- data.frame(gn_list, MEcolors)
  for (i in 1:nrow(MEcolors)){
    if (any(abs(kMEs[MEcolors$gn_list[i],]) > thrsh_kME) != TRUE){
      MEcolors$MEcolors[i] <- "grey"
    }
  }
  dynamicColors_2 <- MEcolors$MEcolors
  table(dynamicColors_2)
  #MEList$validColors <- MEcolors$MEcolors
  
  
print("saving the R objects", quote = FALSE)
    save(list=c("merge","dynamicColors_2", "geneTree", "kMEs"), file="WGCNA_objects_ms10.Rdata")

#extract metadata==================================================================
#extract tissue of origin==========================================================  
print("now extracting relevant compound data", quote = FALSE)  
merge$colors <- dynamicColors_2
ME_data <- data.frame(MEcolors$gn_list, merge$colors)
colnames(ME_data)<-c("Accession", "module")
rownames(ME_data)<-ME_data$Accession

#write all modules
print("now saving all modules' eigenvalues", quote = FALSE)
write.csv(merge$newMEs, file = eigenfile)

tissue_origin <- ME_data[1]
tissue_origin$Accession <- rownames(tissue_origin)
colnames(tissue_origin)<-"tissue_origin"
for (i in 1:nrow(tissue_origin)){
  tissue_origin$tissue_origin[i]<-str_split(tissue_origin$tissue_origin[i], pattern="_")[[1]][1]
}
colnames(tissue_origin)<-"Compound"
ME_data <- cbind(ME_data, tissue_origin)

#add the kME values for the pruned modules to the file=============================================
print("adding the kME values for the assigned modules", quote = FALSE)
print("note- the kME assignments may not be maximum for", quote = FALSE)
print("the assigned module. This is because the dendrogram", quote = FALSE)
print("is factored into the module assignments.", quote = FALSE)
kval <- ME_data[1:2]
colnames(kval)<-c("kME_module", "module")
for (i in 1:nrow(kval)){
  kval$kME_module[i] <- kMEs[kval$kME_module[i], paste0("kME",kval$module[i])]
}
kval<-kval[-2]
ME_data <- cbind(ME_data, kval)
ME_data <- subset(ME_data, module != "grey")

print("saving the modules", quote = FALSE)  
write.csv(ME_data, file = allmodules, row.names = FALSE)

#show the dendrogram with text================================================================
print("displaying the new module eigengene clustering", quote = FALSE)

pdf(file = "Dendrogram_with_colors_assigned.pdf", width = 12, height = 4)
pltree <- plotDendroAndColors(geneTree, colors=dynamicColors_2,
                              "Merged & pruned", rowText = dynamicColors_2,
                              rowTextAlignment = "center",addTextGuide =TRUE,
                              dendroLabels = FALSE, hang = 0.03,
                              rowTextIgnore = "grey",
                              addGuide = FALSE, guideHang = 0.05)
dev.off()

#base::readline(prompt="Press Enter to continue: ")
}








#end the function====================================================
#
}
