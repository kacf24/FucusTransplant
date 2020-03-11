#Load in required packages
library(readr)
library(ggplot2)
library(tidyverse)

#Set correct working directory and read in files
setwd("~/Desktop/ForGithub")
tax <- read_csv("taxonomy.csv")
ASVs <- read_csv("ASVs.csv")
meta <- read_csv("metadata.csv")

#Coverts data into a long format for easy subsetting
data_long <- gather(ASVs, Site, measurement, colnames(ASVs[,2:ncol(ASVs)]), factor_key=TRUE)
data_long_met<-merge(data_long,meta,by.y="samples",by.x="Site")
data_long_met_tax<-merge(data_long_met,tax,by.x="samples",by.y="ASVid")

#Subset to correct data, in this case removes water column samples, holdfast samples, and days 1/2
s_data_long_met_tax<-subset(data_long_met_tax, treatment != "none")
s_data_long_met_tax<-subset(s_data_long_met_tax,tissue!="holdfast"&day=="4")

#Creates a dataframe and transforms data into presence/absence 
pa_data <- s_data_long_met_tax
pa_data[which(pa_data$measurement<1),3] <- 0
pa_data[which(pa_data$measurement>=1),3] <- 1

#calculates the percentage of samples an ASV in present in 
agg_padata <- do.call(data.frame, aggregate(pa_data$measurement, list(pa_data$samples,pa_data$treatment), FUN = function(x) c(mn = mean(x))))

#Creates a cutoff between 0 and 1
agg_padata<-subset(agg_padata,agg_padata$x >= 0.5)

#Creates long and wide dataframes, which are both useful for different types of analyses
agg_padata$x3 <- 1
agg_padata<-agg_padata[,-c(3)]
x<-spread(agg_padata,Group.2,x3)
names(agg_padata) <- c("ASVid","Treatment","Present")
vc<-merge(x,tax,by.x="Group.1",by.y="ASVid")

#Save your results
write.csv(vc,"TESTwideDataRV.csv",row.names = FALSE)
write.csv(agg_padata,"TESTlongDataRV.csv")

