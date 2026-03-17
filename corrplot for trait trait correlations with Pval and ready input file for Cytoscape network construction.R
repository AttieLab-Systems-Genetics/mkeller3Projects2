###############################################
# Correlation + P-value + Cytoscape Edge List
# Fixed: consistent column names & style hints
###############################################

# Install packages if missing
packages <- c("tidyverse", "corrplot", "gdata", "gplots", "limma")
installed <- rownames(installed.packages())
for (p in packages) {
  if (!(p %in% installed)) install.packages(p, dependencies = TRUE)
}

library(tidyverse)
library(corrplot)
library(gdata)
library(gplots)
library(limma)
library(rstudioapi)

###############################################
# Load data
###############################################

setwd(selectDirectory())
data <- data.frame(read.csv(file.choose()), stringsAsFactors = FALSE)
rownames(data) <- data$ID
data2 <- data[, 4:ncol(data)]

###############################################
# Compute P-values
###############################################

cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

M_trait <- cor(data2, use = "pairwise.complete.obs")
p.mat <- cor.mtest(data2)

###############################################
# Save correlation and P-value matrices
###############################################

flnm <- readline(prompt = "Enter output filename for correlation matrix: ")
write.csv(M_trait, file = flnm, row.names = TRUE, quote = FALSE)

flnm <- readline(prompt = "Enter output filename for P-value matrix: ")
write.csv(p.mat, file = flnm, row.names = TRUE, quote = FALSE)

###############################################
# P-value cutoff
###############################################

get_p_cutoff <- function() {
  cat("\nSelect a P-value cutoff:\n")
  cat("  1 = 0.05\n")
  cat("  2 = 0.01\n")
  cat("  3 = 0.001\n")
  cat("  4 = custom numeric value\n\n")
  
  choice <- readline(prompt = "Enter choice (1–4): ")
  choice <- suppressWarnings(as.numeric(choice))
  
  if (is.na(choice) || choice < 1 || choice > 4)
    stop("Invalid selection. Please enter 1, 2, 3, or 4.")
  
  if (choice == 1) return(0.05)
  if (choice == 2) return(0.01)
  if (choice == 3) return(0.001)
  
  custom <- readline(prompt = "Enter custom numeric P-value cutoff: ")
  custom <- suppressWarnings(as.numeric(custom))
  if (is.na(custom)) stop("Invalid numeric cutoff.")
  
  return(custom)
}

p_cutoff <- get_p_cutoff()
cat("Using P-value cutoff:", p_cutoff, "\n")

###############################################
# Create Cytoscape edge list
###############################################

create_cytoscape_edges <- function(cor_mat, p_mat, p_thresh = 0.05) {
  
  cor_df <- cor_mat %>%
    as.data.frame() %>%
    rownames_to_column("traitA") %>%
    pivot_longer(-traitA, names_to = "traitB", values_to = "correlation")
  
  p_df <- p_mat %>%
    as.data.frame() %>%
    rownames_to_column("traitA") %>%
    pivot_longer(-traitA, names_to = "traitB", values_to = "pvalue")
  
  merged <- left_join(cor_df, p_df, by = c("traitA", "traitB")) %>%
    filter(traitA < traitB) %>%              # upper triangle only
    filter(pvalue <= p_thresh) %>%           # P-value filter
    mutate(
      absCorrelation = abs(correlation),     # width
      sign = ifelse(correlation >= 0, "positive", "negative"),
      sign_numeric = ifelse(correlation >= 0, 1, -1), # color
      minusLog10P = -log10(pvalue)           # transparency
    ) %>%
    rename(source = traitA, target = traitB) # consistent names
  
  return(merged)
}

edge_list <- create_cytoscape_edges(M_trait, p.mat, p_cutoff)
cat("Number of significant edges:", nrow(edge_list), "\n")

###############################################
# Save Cytoscape edge list (quote-free)
###############################################

flnm <- readline(prompt = "Enter filename for Cytoscape edge list: ")
write.csv(edge_list, file = flnm, row.names = FALSE, quote = FALSE)

###############################################
# Save style hints CSV
###############################################

style_file <- readline(prompt = "Enter filename for Cytoscape style hints CSV: ")

style_hints <- edge_list %>%
  select(source, target, sign_numeric, absCorrelation, minusLog10P) %>%
  rename(
    color = sign_numeric,
    width = absCorrelation,
    transparency = minusLog10P
  )

write.csv(style_hints, file = style_file, row.names = FALSE, quote = FALSE)

cat("\nCytoscape edge list saved as:", flnm, "\n")
cat("Cytoscape style hints file saved as:", style_file, "\n")
cat("Done.\n")




#Corr plot function to display relationships between MEs and clinical traits

col <- colorRampPalette(c("#4477AA", "#77AADD","#FFFFFF", "#EE9988","#BB4444"))

corrplot(M_trait, 
         method = "square", 
         #col=COL2('BrBG',200), 
         col = col(200),
         order = "hclust", 
         type = "lower", 
         tl.col="black",
         p.mat = p.mat, 
         sig.level = .01, 
         insig = "blank",
         diag = FALSE)













#this function can also be used to generate a different kind of heatmap
#draw heatmap of correlations
hm <- heatmap.2(M_trait,
                symm=TRUE,
                dendrogram="both",
                col='bluered',
                #distfun=function(c) as.dist(1-c),
                hclustfun=hclust,
                density.info='density',
                key.title=NA,
                keysize = 1,
                margins=c(12,12),
                tracecol="black",
                trace='none')


#output files for clustered matrices
flnm <- base::readline(prompt = "Enter filename for data: ")
write.csv(M_trait[hm$rowInd, hm$colInd],
          file = flnm,
          row.names = TRUE)

flnm <- base::readline(prompt = "Enter filename for data: ")
write.csv(p.mat[hm$rowInd, hm$colInd],
          file = flnm,
          row.names = TRUE)





########################################################################
#These functions can be used to regress one or more variables, e.g., sex
########################################################################


#create function to regress out Sex
regress_out_sex <- function(mat, Sex) {
  apply(mat, 2, function(y) {
    out <- rep(NA, length(Sex))
    isna <- !is.na(y)
    out[isna] <- resid(lm(y ~ Sex))
    out
  })
}


#create a new data matrix with sex regressed out

Sex <- as.factor(data$Sex)

data3 <- regress_out_sex(data2, Sex)
data3 <- as.data.frame(data3)
data4 <- cbind(data[,1:2], data3)

#output adjusted matrix
flnm <- base::readline(prompt = "Enter filename for data: ")
write.csv(data4, file = flnm, row.names = TRUE)

#recompute correlation and p-value matrices from adjusted sex adj data

#compute correlation and p-values
M_trait_sex <- cor(data3, use = "pairwise.complete.obs")
p.mat_sex <- cor.mtest(data3)

#view data
head(M_trait_sex[, 1:5])
head(p.mat_sex[, 1:5])

#draw heatmap of correlations for sex adj data matrix
hm2 <- heatmap.2(M_trait_sex,
                symm=TRUE,
                dendrogram="both",
                col='bluered',
                #distfun=function(c) as.dist(1-c),
                hclustfun=hclust,
                density.info='density',
                key.title=NA,
                keysize = 1,
                margins=c(12,12),
                tracecol="black",
                trace='none')



par(mfrow=c(2,2))
boxplot(data$ins_14~Sex, col=c("pink", "royalblue"), main="Unadjusted")
boxplot(data3$ins_14~Sex, col=c("pink", "royalblue"), main="Adjusted for sex")
boxplot(data$MEturquoise~Sex, col=c("pink", "royalblue"), main="Unadjusted")
boxplot(data3$MEturquoise~Sex, col=c("pink", "royalblue"), main="Adjusted for sex")



#output files for clustered matrices after sex adjustment
flnm <- base::readline(prompt = "Enter filename for data: ")
write.csv(M_trait_sex[hm2$rowInd, hm2$colInd],
          file = flnm,
          row.names = TRUE)

flnm <- base::readline(prompt = "Enter filename for data: ")
write.csv(p.mat_sex[hm2$rowInd, hm2$colInd],
          file = flnm,
          row.names = TRUE)


#Compute ME vs ME correlation matrix for directed network in Cytoscape
ME <- data.frame(read.csv(file.choose()), stringsAsFactors = FALSE)
rownames(ME) <- data$mouse
ME2 <- ME[,4:ncol(ME)]

#compute correlation and p-values
ME2_corr <- cor(ME2, use = "pairwise.complete.obs")
p.mat_ME2 <- cor.mtest(ME2)

#draw heatmap of ME:ME correlations
hm2 <- heatmap.2(ME2_corr,
                 symm=TRUE,
                 dendrogram="both",
                 col='bluered',
                 #distfun=function(c) as.dist(1-c),
                 hclustfun=hclust,
                 density.info='density',
                 key.title=NA,
                 keysize = 1,
                 margins=c(12,12),
                 tracecol="black",
                 trace='none')

#output Correlation and P-value matrices for MEs
flnm <- base::readline(prompt = "Enter filename for data: ")
write.csv(ME2_corr, file = flnm, row.names = TRUE)

flnm <- base::readline(prompt = "Enter filename for data: ")
write.csv(p.mat_ME2, file = flnm, row.names = TRUE)

#Regress out sex differences among MEs and recompute correlation matrix

MESex <- as.factor(ME$Sex)

ME2_sex <- regress_out_sex(ME2, Sex)
ME2_sex <- as.data.frame(ME2_sex)
ME2_sex2 <- cbind(ME[,1:3], ME2_sex)

#recompute correlation and p-value matrices from sex adj data

#compute correlation and p-values
ME2_sex_corr <- cor(ME2_sex, use = "pairwise.complete.obs")
p.mat_ME2_sex <- cor.mtest(ME2_sex)


#draw heatmap of ME:ME correlations
hm3 <- heatmap.2(ME2_sex_corr,
                 symm=TRUE,
                 dendrogram="both",
                 col='bluered',
                 #distfun=function(c) as.dist(1-c),
                 hclustfun=hclust,
                 density.info='density',
                 key.title=NA,
                 keysize = 1,
                 margins=c(12,12),
                 tracecol="black",
                 trace='none')

#output files for clustered matrices after sex adjustment
flnm <- base::readline(prompt = "Enter filename for data: ")
write.csv(ME2_sex_corr[hm3$rowInd, hm3$colInd],
          file = flnm,
          row.names = TRUE)

flnm <- base::readline(prompt = "Enter filename for data: ")
write.csv(p.mat_ME2_sex[hm3$rowInd, hm3$colInd],
          file = flnm,
          row.names = TRUE)

