#Load in required packages
library(mvabund)
library(readr)

#Set correct working directory and read in files
setwd("~/Desktop/ForGithub")
tASVs <- read_csv("ASVs.csv")
meta <- read_csv("metadata.csv")

#Creates a pseudovariable for interaction term.
meta$interaction<-paste0(meta$day,meta$treatment)

#Merges together files
s_tASVs<-merge(meta,tASVs)

#Subset to the correct test, in this case it is the vegetative tissue from the transplants with three sample days
s_tASVs2<-subset(s_tASVs,tissue=="vegetative" & treatment !="none_Fd" & tissue != "na" & treatment != "none_Fv" & treatment != "none_Fs")

#Takes out extra metadata
x<-s_tASVs2[,-c(1:11)]

#Converts the data to presence/absence
for(j in 1:ncol(x)){
  for(k in 1:nrow(x)){
    if(x[k,j] > 0){x[k,j] <- 1}
  }
}


#Sets up the mvabund object. It only takes numeric data, not factors
Herb_spp <- mvabund(x)

#Create your MGLM model, formula is essentially the same as normal glm. Family is the distribtuion that is being used. I use binomial for P/A data, and negative binomial (default) for abundance data 
mod3 <- manyglm(Herb_spp ~s_tASVs2$transect+s_tASVs2$day+s_tASVs2$treatment+s_tASVs2$day:s_tASVs2$treatment)

#Your diagnostic plot, you're looking for something without a very strong pattern (i.e., not linear, exponential, quadratic). Can be a bit slow to display.
#plot(mod3)



#This is what gives significance, specify your model "mod3", how you want it do to pairwise computations, and if you want adjusted pvalues for univariate testing. Show time can be nice, but often breaks.
x<-anova(mod3,p.uni="adjusted",pairwise.comp=s_tASVs2$interaction,show.time = "all")

