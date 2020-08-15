library(readr)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(tidyr)
#setwd("~/Desktop/barplots")
tax <- read_csv("newTaxonomy.csv")
ASVs <- read_csv("ASVs.csv")
meta <- read_csv("metadata.csv")


`%notin%` <- Negate(`%in%`)

for(i in 2:ncol(ASVs)){
 ASVs[,i]<-ASVs[,i]/colSums(ASVs)[i]}

min(colSums(ASVs))
max(colSums(ASVs[,2:227]))

data_long <- gather(ASVs, Site, measurement, colnames(ASVs[,2:227]), factor_key=TRUE)
data_long_met<-merge(data_long,meta,by.y="samples",by.x="Site")
data_long_met_tax<-merge(data_long_met,tax,by.x="samples",by.y="ASVid")

s_data_long_met_tax<-subset(data_long_met_tax,day=="4")


test<-do.call(data.frame, aggregate(s_data_long_met_tax$measurement, list(s_data_long_met_tax$Class), FUN = function(x) c(mn = sum(x))))


non_zero_data_long_met_tax<-subset(s_data_long_met_tax,s_data_long_met_tax$Class %notin% test[which(test$x<1),1])

zero_data_long_met_tax<-subset(s_data_long_met_tax,s_data_long_met_tax$Class %in% test[which(test$x<1),1])

zero_data_long_met_tax$Class<-"Other"

combined<-rbind(non_zero_data_long_met_tax,zero_data_long_met_tax)

yeez <- subset(meta,day =="4")
yeezy <- do.call(data.frame, aggregate(yeez$sample_ID, list(yeez$treatment,yeez$tissue), FUN = function(x) c(mn = length(x))))
yeezy$combo <- paste0(yeezy$Group.1,yeezy$Group.2)

test<-do.call(data.frame, aggregate(combined$measurement, list(combined$Class,combined$treatment,combined$tissue), FUN = function(x) c(mn = sum(x))))
test$combo <- paste0(test$Group.2,test$Group.3)
test_yeezy<-merge(test,yeezy,by="combo")
test_yeezy<-test_yeezy[,-c(1,6,7)]
names(test_yeezy) <- c("Class","Treatment","Tissue","SumOfProportions","NumberofSamples")
test_yeezy$MeanProportion <- test_yeezy$SumOfProportions/test_yeezy$NumberofSamples
test_yeezy$combo<-paste0(test_yeezy$Class,test_yeezy$Treatment,test_yeezy$Tissue)


teezy2<-do.call(data.frame, aggregate(combined$measurement, list(combined$sample_ID,combined$Class,combined$tissue,combined$treatment), FUN = function(x) c(mn = sum(x))))
teezy2$combo <- paste0(teezy2$Group.2,teezy2$Group.4,teezy2$Group.3)
teezy3<-merge(teezy2,test_yeezy,by="combo")
teezy3$dev <- teezy3$MeanProportion - teezy3$x
teezy3$sqdev <- teezy3$dev^2
teezy4<-do.call(data.frame, aggregate(teezy3$sqdev, list(teezy3$Class,teezy3$Treatment,teezy3$Tissue), FUN = function(x) c(mn = sum(x))))
names(teezy4) <- c("Class","Treatment","Tissue","SumOfSquareDeviances")
teezy4$combo <- paste0(teezy4$Treatment,teezy4$Tissue)
finalish<-merge(teezy4,yeezy,by="combo")

finalish$stdev <- sqrt(finalish$SumOfSquareDeviances/(finalish$x-1))
finalish$sem <- finalish$stdev/sqrt(finalish$x)
finalish$combo <-paste0(finalish$Class,finalish$combo)
final<-merge(test_yeezy,finalish,"combo")


final<-subset(final,Tissue.x!="na")
final$Class.x <- as.character(final$Class.x)
final$Tissue.x <- as.character(final$Tissue.x)
final[which(final$Class.x=="Alphaproteobacteria"),2] <- "Alphaproteobac"
final[which(final$Class.x=="Bacteroidia"),2] <- "B"
final[which(final$Class.x=="Gammaproteobacteria"),2] <- "Gammaproteobac"
final[which(final$Class.x=="Other"),2] <- "O"
final[which(final$Class.x=="Acidimicrobiia"),2] <- "Ac"
final[which(final$Class.x=="Planctomycetacia"),2] <- "P"
final[which(final$Class.x=="Verrucomicrobiae"),2] <- "V"
final[which(final$Tissue.x=="holdfast"),4] <- "Holdfast"
final[which(final$Tissue.x=="reproductive"),4] <- "Receptacle"
final[which(final$Tissue.x=="vegetative"),4] <- "Vegetative"

final$Treatment.x <- as.character(final$Treatment.x)

final$Treatment.x<-factor(final$Treatment.x,levels=c("none_Fs","dry","watered","control","none_Fv","none_Fd"))
final$Class.x <-factor(final$Class.x,levels=c("Acidimicrobiia","Bacteroidia","Alphaproteobac","Gammaproteobac","Other","Planctomycetacia","Verrucomicrobiae"))

tiff("Barplot.tiff", width = 180, height = 180, units = 'mm', res = 300)
ggplot(final,aes(x=Class.x,y=MeanProportion,fill=Treatment.x))+geom_bar(stat="identity",position=position_dodge())+ facet_grid(rows=vars(final$Tissue.x)) + theme_bw() + theme(text = element_text(size=8)) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + scale_fill_manual(name = "Treatment", breaks=c("none_Fs","dry","watered","control","none_Fv","none_Fd") , labels = c("Fs","D","W","PC","Fv","Fd"), values = c("#F0E442","#E69F00","#009E73","black","#0072B2","#D55E00","gray")) +  geom_errorbar(aes(ymin=MeanProportion-sem, ymax=MeanProportion+sem), position=position_dodge(.9))+ xlab("Class") + ylab("Mean Proportion")
dev.off()
+geom_boxplot()+xlab("Class") + facet_grid(cols = vars(final$Treatment.x),rows=vars(final$Tissue.x)) 

ggplot(df2, aes(x=dose, y=len, fill=supp)) + 
  geom_bar(stat="identity", color="black", position=position_dodge()) +
  geom_errorbar(aes(ymin=len, ymax=len+sd), width=.2, position=position_dodge(.9)) 


ggplot(data=test,aes(x=test$tissue,y=test$mn,fill=test$class))+geom_bar(position="fill",stat = "identity")+ facet_grid(cols = vars(test$treatment))+scale_fill_manual(values = colz,name="Class")  + theme_classic()+ theme(strip.text.x = element_text(size = 16,face="bold"),axis.text.y = element_text(size=10,face="bold"),axis.text.x = element_text(size=12,face="bold"),axis.title.x = element_text(size=16, face="bold"),axis.title.y = element_text(size=16, face="bold"),legend.text=element_text(size=10,face="bold"),legend.title = element_text(size = 12,face="bold"))+ xlab("Tissue")+ylab("Mean Proportion")
dev.off()
