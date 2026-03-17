library(rstudioapi)
library(tidyverse)

#set directory
setwd(selectDirectory())

#load data
data <- data.frame(read.csv(file.choose()), stringsAsFactors = FALSE)
rownames(data)<-data$mouse
data2<-data[-1]


#Define rankz fundtion
rankz = function(x) { 
  x = rank(x, na.last = "keep", ties.method = "average") / (sum(!is.na(x)) + 1)
  return(qnorm(x))
}


#Check distributions of untransformed data
hist(data2[,10])
hist(data2[,20])
hist(data2[,30])
hist(data2[,40])
hist(data2[,50])
hist(data2[,60])
hist(data2[,70])
hist(data2[,80])
hist(data2[,1])


#Perform rankz transformation of all data
for (i in 1:ncol(data2)){
  data2[,i] <- rankz(data2[,i])
}



flnm <- base::readline(prompt = "Enter filename for data: ")

write.csv(data2, file = flnm, row.names=TRUE)
