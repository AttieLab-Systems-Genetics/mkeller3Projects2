library(rstudioapi)
library(tidyverse)

#set directory
setwd(selectDirectory())

#load data
data <- data.frame(read.csv(file.choose()), stringsAsFactors = FALSE)
rownames(data) <- data$Ensembl.Symbol
data2 <- data[,-1]
data3 <- t(data2)

#create rankz fundtion
rankz = function(x) { 
  x = rank(x, na.last = "keep", ties.method = "average") / (sum(!is.na(x)) + 1)
  return(qnorm(x))
}

#setup plotting window
par(mfrow = c(2, 2))

#view histograms of non-transformed data
hist(data3[,1])
hist(data3[,100])
hist(data3[,200])
hist(data3[,300])
hist(data3[,400])
hist(data3[,500])
hist(data3[,600])
hist(data3[,700])
hist(data3[,800])
hist(data3[,900])
hist(data3[,1000])
hist(data3[,2000])
hist(data3[,3000])
hist(data3[,4000])
hist(data3[,5000])
hist(data3[,6000])
hist(data3[,7000])
hist(data3[,8000])
hist(data3[,9000])
hist(data3[,10000])

#perform rankz transformation on all data
for (i in 1:ncol(data3)){
  data3[,i] <- rankz(data3[,i])
}

#view histograms of transformed data to verify rankz transformed correctly
hist(data3[,1])
hist(data3[,100])
hist(data3[,200])
hist(data3[,300])
hist(data3[,400])
hist(data3[,500])
hist(data3[,600])
hist(data3[,700])
hist(data3[,800])
hist(data3[,900])
hist(data3[,1000])
hist(data3[,2000])
hist(data3[,3000])
hist(data3[,4000])
hist(data3[,5000])
hist(data3[,6000])
hist(data3[,7000])
hist(data3[,8000])
hist(data3[,9000])
hist(data3[,10000])

#transposed data so genes are rows and mice as columns
data4<-t(data3)




#output a csv file with transformed data
flnm <- base::readline(prompt = "Enter filename for data: ")
write.csv(data4, file = flnm)
