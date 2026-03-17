#install packages and librares
install.packages("corrplot")
library(rstudioapi)
library(tidyverse)
library(corrplot)
library(gdata)
library(gplots)
library(limma)

#set directory
setwd(selectDirectory())

#load data
data <- data.frame(read.csv(file.choose()), stringsAsFactors = FALSE)
rownames(data) <- data$ID
data2 <- data[,6:ncol(data)]

# create function with R cor.test
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
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

#compute correlation and p-values
M_trait <- cor(data2, use = "pairwise.complete.obs")
p.mat <- cor.mtest(data2)

#view data
head(M_trait[, 1:5])
head(p.mat[, 1:5])


#Corr plot function to display relationships between MEs and clinical traits

col <- colorRampPalette(c("#4477AA", "#77AADD","#FFFFFF", "#EE9988","#BB4444"))

corrplot(M_trait, 
         method = "circle", 
         #col=COL2('BrBG',200), 
         col = col(200),
         order = "hclust", 
         type = "lower", 
         tl.col="black",
         p.mat = p.mat, 
         sig.level = 1E-2, 
         insig = "blank",
         diag = FALSE)


#output files for clustered matrices
flnm <- base::readline(prompt = "Enter filename for data: ")
write.csv(M_trait,
          file = flnm,
          row.names = TRUE)

flnm <- base::readline(prompt = "Enter filename for data: ")
write.csv(p.mat,
          file = flnm,
          row.names = TRUE)









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

