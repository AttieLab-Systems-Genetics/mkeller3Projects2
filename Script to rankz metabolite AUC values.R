library(rstudioapi)
library(tidyverse)

#set directory
setwd(selectDirectory())

#load data
data <- data.frame(read.csv(file.choose()), stringsAsFactors = FALSE)
rownames(data) <- data$X
data2<-data[-1]


#Define rankz fundtion
rankz = function(x) { 
  x = rank(x, na.last = "keep", ties.method = "average") / (sum(!is.na(x)) + 1)
  return(qnorm(x))
}


par(mfrow = c(2, 2))

#Check distributions of untransformed data
hist(data2[,10])
hist(data2[,20])
hist(data2[,30])
hist(data2[,40])
hist(data2[,50])
hist(data2[,60])
hist(data2[,70])
hist(data2[,80])
hist(data2[,90])
hist(data2[,100])
hist(data2[,150])
hist(data2[,200])
hist(data2[,250])
hist(data2[,300])
hist(data2[,350])
hist(data2[,355])




#Perform rankz transformation of all data
for (i in 1:ncol(data2)){
  data2[,i] <- rankz(data2[,i])
}

#Check distributions of transformed data
hist(data2[,10])
hist(data2[,20])
hist(data2[,30])
hist(data2[,40])
hist(data2[,50])
hist(data2[,60])
hist(data2[,70])
hist(data2[,80])
hist(data2[,90])
hist(data2[,100])
hist(data2[,150])
hist(data2[,200])
hist(data2[,250])
hist(data2[,300])
hist(data2[,350])
hist(data2[,355])


flnm <- base::readline(prompt = "Enter filename for data: ")

write.csv(data2, file = flnm, row.names=TRUE)
