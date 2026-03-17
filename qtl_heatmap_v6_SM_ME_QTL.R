#' New QTL heatmap function for mouse qtl plotting
#' 
#' @author Chris Emfinger, PhD. Desperately needs a map
#' 
#' @param qtl.temp is the combined qtl object. It needs to have markers
#' that match the map's marker names
#' 
#' @param mrkrs is the marker object. rownames of the qtl object must
#' match the marker id in the qtl object
#' 
#' @param low.thr is the low threshold for color display. Every value
#' below this is set to 0 for display
#' 
#' @param high.thr is the high threshold for color display. Every value
#' above this is set to the high threshold.
#' 
#' @param mid.thr is the middle threshold for the color display
#' this is mostly 
#' 
#' @param chr_breaks is the csv file with the chromosomal breaks for the genome of inter
#' 
#' @param set2 is the value you want everything below your threshold set to
#' 
#' @param lowcolor is the color for the gradient for the plotting fills
#' 
#' @param midcolor is the color for the middle of the plotting fill gradient
#' 
#' @param highcolor is the color for the high value of the plotting fill gradient
#' 
#' this assumes an input from the AutoQTL script. 
#' 

qtl_heatmap <- function(qtl.temp, 
                        mrkrs, 
                        low.thr = 5.5,
                        mid.thr = 6,
                        high.thr = 7, 
                        chr_breaks,
                        set2 = 0,
                        lowcolor = "white",
                        midcolor = 'orange',
                        highcolor = 'darkred'){
  
#load libraries =======================================================
  library(rstudioapi)
  library(tidyverse)
  library(ggplot2)
  library(RColorBrewer)
  library(stringr)
  library(scales)

# Note, clustering method = "ward.D", "ward.D2", "complete", or "average"  
  
#data processing=======================================================
  #get the rowname/colname info
  clnms <- colnames(qtl.temp)
  markers_ID <- rownames(qtl.temp)
  qtl.temp <- data.frame(qtl.temp)
  
  #define clustering order
  if (ncol(qtl.temp) > 3){
    qtl.cluster = hclust(as.dist(1.0 - cor(qtl.temp)), method = "ward.D2")
    qtl.temp = qtl.temp[,qtl.cluster$order]
  }
  if (ncol(qtl.temp) <= 3){
    cat("you input has too few qtl to cluster \n"
    )
  }
  
  #set to thresholds
  clnms <- colnames(qtl.temp)
  qtl.hi <- which(qtl.temp > high.thr)
  qtl.lo <- which(qtl.temp < low.thr)
  qtl.temp <- as.matrix(qtl.temp)
  qtl.temp[qtl.hi] <- high.thr
  qtl.temp[qtl.lo]<- set2
  qtl.temp <- data.frame(qtl.temp)
  #colnames(qtl.temp)<-clnms
  
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
  
  qtl.temp$order <- qtl.temp$chr
  
  qtl.temp[which(qtl.temp$chr=="X"),"order"]<-20
  qtl.temp[which(qtl.temp$chr=="Y"),"order"]<-21
  qtl.temp[which(qtl.temp$chr=="M"),"order"]<-22
  qtl.temp$chr <- qtl.temp$order
  qtl.temp$order <- as.numeric(qtl.temp$order)
  
  #add the max positions for each chromosome
  qtl.pos <- qtl.temp[1:nrow(chr_breaks),]
  for (i in 1:ncol(qtl.pos)){
    qtl.pos[,i]<-0
  }
  qtl.pos$markers <- chr_breaks$chr
  qtl.pos$chr <- chr_breaks$chr
  qtl.pos$position <- chr_breaks$initial.end
  qtl.pos$order <- chr_breaks$chr
  qtl.pos$position <- qtl.pos$position/(10^6)
  qtl.pos[which(qtl.pos$chr=="X"),"order"]<-20
  qtl.pos[which(qtl.pos$chr=="Y"),"order"]<-21
  qtl.pos[which(qtl.pos$chr=="M"),"order"]<-22
  qtl.pos$order <- as.numeric(qtl.pos$order)
  
  qtl.temp <- rbind(qtl.temp, qtl.pos)
  
  #create the plotting object
  qtl_plot_obj <- qtl.temp %>% 
    
    # Compute chromosome size
    group_by(order) %>% 
    summarise(chr_len=max(position)) %>% 
    
    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(chr_len)-chr_len) %>%
    select(-chr_len) %>%
    
    
    # Add this info to the initial dataset
    left_join(qtl.temp, ., by=c("order"="order")) %>%
    
    # Add a cumulative position of each SNP
    arrange(order, position) %>%
    mutate( BPcum=position+tot)
  
  #define the initial start points in x
  qtl_plot_obj$x <- qtl_plot_obj$BPcum
  qtl_plot_obj <- qtl_plot_obj %>%
    mutate(x = lag(BPcum, default = first(BPcum)))
  qtl_plot_obj$x[1]<-0
  
  #x axis mods: we do not want to display the cumulative 
  #position of SNP in bp, but just show the chromosome name instead.
  axisdf = qtl_plot_obj %>% group_by(order) %>% 
    summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
  axisdf$order[which(axisdf$order==20)]<-"X"
  axisdf$order[which(axisdf$order==21)]<-"Y"
  axisdf$order[which(axisdf$order==22)]<-"M"
  axisdf$order <- factor(axisdf$order, levels=c(as.character(1:19),"X","Y","M"))
  
  qtl_plot_obj2 <- qtl_plot_obj[-which(colnames(qtl_plot_obj) %in% c("position","tot"))]
  
  #note- the melt function will not use numeric elements as IDs
  #so you must convert the BPcum column to characters and back
  qtl_plot_obj2$BPcum <- as.character(qtl_plot_obj2$BPcum)
  qtl_plot_obj2$x <- as.character(qtl_plot_obj2$x)
  qtl_plot_obj2$order <- as.character(qtl_plot_obj2$order)
  melted_QTL_obj <- reshape2::melt(qtl_plot_obj2, ids=c("markers","chr", 
                                                   "order",
                                                   "BPcum","x"))
  melted_QTL_obj$BPcum <- as.numeric(melted_QTL_obj$BPcum)
  melted_QTL_obj$x <- as.numeric(melted_QTL_obj$x)
  chrs <- c(1:19,"X","Y","M")
  melted_QTL_obj <-data.frame(melted_QTL_obj, stringsAsFactors = FALSE)
  melted_QTL_obj$ymin <- 0
  melted_QTL_obj$ymax <- 1
  melted_subset <- subset(melted_QTL_obj, !(markers %in% chrs))
  
  plot_heatmap_qtl <- ggplot(data = melted_subset) + 
    
    #plot the marker segments as boxes
    geom_rect(data= melted_subset, 
              mapping=aes(xmin = melted_subset$x, 
              xmax = melted_subset$BPcum, 
              ymin = melted_subset$ymin, 
              ymax = melted_subset$ymax, 
              fill = melted_subset$value),
              alpha=1) +
    
    #remove extraneous grids/labels
    theme_bw() +
    
    theme( 
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.background = element_rect(colour = "black", 
                                      linewidth = 0.25)
    )+
    
    #add back y axis
    theme(
      axis.line.y = element_line(colour = "black", linewidth = 0.25),
    ) +
    
    #change color scheme
    scale_fill_gradient2(low = lowcolor, 
                         mid = midcolor, 
                         high = highcolor, 
                         midpoint = mid.thr,
                         guide="colorbar")+
    
    #add chromosomal labels
    scale_x_continuous( label = axisdf$order, 
                        breaks= axisdf$center,
                        expand = c(0,2))+
    
    #add the vertical lines for each chromosome
    geom_vline(xintercept = subset(melted_QTL_obj, markers %in% chrs)$BPcum,
               color = "black", linewidth = 0.25)+
  
    #facet by scan  
    facet_grid(rows=vars(variable))+
    
    #set spacing between panels to 0
    theme(panel.spacing.y=unit(0, "lines"))+
    
    #rotate the faceting labels
    theme(strip.text.y = element_text(angle = 0))+
    
    #change the legend label
    labs(fill="LOD")
  
  return(plot_heatmap_qtl)
  
}






