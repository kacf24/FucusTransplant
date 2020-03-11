#Load in required packages
library(readr)
library(ggplot2)
library(tidyverse)
library(dplyr)

#Set correct working directory and read in files
setwd("~/Desktop/ForGithub")
tax <- read_csv("taxonomy.csv")
ASVs <- read_csv("ASVs.csv")
meta <- read_csv("metadata.csv")

#A custom function
`%notin%` <- Negate(`%in%`)

#Converts data into proprotional
for(i in 2:ncol(ASVs)){
 ASVs[,i]<-ASVs[,i]/colSums(ASVs)[i]}

#Coverts data into a long format for easy subsetting
data_long <- gather(ASVs, Site, measurement, colnames(ASVs[,2:227]), factor_key=TRUE)
data_long_met<-merge(data_long,meta,by.y="samples",by.x="Site")
data_long_met_tax<-merge(data_long_met,tax,by.x="samples",by.y="ASVid")

#Subset to correct data
s_data_long_met_tax<-subset(data_long_met_tax,day=="4")

#Finds overall summed proportions of each
test<-do.call(data.frame, aggregate(s_data_long_met_tax$measurement, list(s_data_long_met_tax$Class), FUN = function(x) c(mn = sum(x))))

#Gets rid of things that are less than 1 and lumps them into an "Other" categroy
non_zero_data_long_met_tax<-subset(s_data_long_met_tax,s_data_long_met_tax$Class %notin% test[which(test$x<1),1])
zero_data_long_met_tax<-subset(s_data_long_met_tax,s_data_long_met_tax$Class %in% test[which(test$x<1),1])
zero_data_long_met_tax$Class<-"Other"
combined<-rbind(non_zero_data_long_met_tax,zero_data_long_met_tax)

#Sums each class for tissues and treatments
test<-do.call(data.frame, aggregate(combined$measurement, list(combined$Class,combined$tissue,combined$treatment), FUN = function(x) c(mn = sum(x))))
colnames(test)<-c("class","tissue","treatment","x")

#Finds the (n=) for each tissue and treatment and calculates mean
day_4 <- subset(meta,day =="4")
day_4_data <- do.call(data.frame, aggregate(day_4$sample_ID, list(day_4$treatment,day_4$tissue), FUN = function(x) c(mn = length(x))))
day_4_data$combo <- paste0(day_4_data$Group.1,day_4_data$Group.2)
test$combo <- paste0(test$treatment,test$tissue)
test_day_4_data<-merge(test,day_4_data,by="combo")
test$mn <- test_day_4_data$x.x/test_day_4_data$x.y


#Renames Tissues and Treatments
test$treatment <- as.character(test$treatment)
test$tissue <- as.character(test$tissue)
test[which(test$treatment=="none" & test$tissue=="na"),2]<-"holdfast"
test[which(test$treatment=="none"),3]<-"WC"
test[which(test$treatment=="dry"),3]<-"D"
test[which(test$treatment=="watered"),3]<-"W"
test[which(test$treatment=="control"),3]<-"PC"
test[which(test$treatment=="none_Fv"),3]<-"Fv"
test[which(test$treatment=="none_Fs"),3]<-"Fs"
test[which(test$treatment=="none_Fd"),3]<-"Fd"

test[which(test$tissue=="holdfast"),2]<-"H"
test[which(test$tissue=="reproductive"),2]<-"R"
test[which(test$tissue=="vegetative"),2]<-"V"


test2<-do.call(data.frame, aggregate(combined$measurement, list(combined$Class), FUN = function(x) c(mn = sum(x))))

#Reorders and renames for legend
test$class<-as.character(test$class)
test$sort <- 1
test[test$class == "Acidimicrobiia",7] <- 1
test[test$class == "Other",7] <- 2
test[test$class == "Sphingobacteriia",7] <- 3
test[test$class == "Verrucomicrobiae",7] <- 4
test[test$class == "Planctomycetacia",7] <- 5
test[test$class == "Flavobacteriia",7] <- 6
test[test$class == "Betaproteobacteria",7] <- 7
test[test$class == "Alphaproteobacteria",7] <- 8
test[test$class == "Gammaproteobacteria",7] <- 9

test[test$class == "Acidimicrobiia",1] <- "Acidimicrobiia (A)"
test[test$class == "Sphingobacteriia",1] <- "Sphingobacteriia (B)"
test[test$class == "Verrucomicrobiae",1] <- "Verrucomicrobiae (V)"
test[test$class == "Planctomycetacia",1] <- "Planctomycetacia (Pl)"
test[test$class == "Flavobacteriia",1] <- "Flavobacteriia (B)"
test[test$class == "Betaproteobacteria",1] <- "Betaproteobacteria (Pr)"
test[test$class == "Alphaproteobacteria",1] <- "Alphaproteobacteria (Pr)"
test[test$class == "Gammaproteobacteria",1] <- "Gammaproteobacteria (Pr)"

test2<-test[order(test[,7]),]
test2$class<-as_factor(test2$class)
test <- test2


#Creates plot
colz <- c("#3800DE","#CCCCCC", "#C5CF08", "#D67240", "#000000", "#E1FF21", "#179BB5", "#8C1ED6", "#25CC4F")
tiff("Barplot.tiff", width = 174, height = 174, units = 'mm', res = 72)
ggplot(data=test,aes(x=test$tissue,y=test$mn,fill=test$class))+geom_bar(position="fill",stat = "identity")+ facet_grid(cols = vars(test$treatment))+scale_fill_manual(values = colz,name="Class")  + theme_classic()+ theme(strip.text.x = element_text(size = 16,face="bold"),axis.text.y = element_text(size=10,face="bold"),axis.text.x = element_text(size=12,face="bold"),axis.title.x = element_text(size=16, face="bold"),axis.title.y = element_text(size=16, face="bold"),legend.text=element_text(size=10,face="bold"),legend.title = element_text(size = 12,face="bold"))+ xlab("Tissue")+ylab("Mean Proportion")
dev.off()
