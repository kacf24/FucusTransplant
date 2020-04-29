#Load the used libraries
library(vegan)
library(svglite)
library(readr)
library(ggplot2)
library(foreach)
library(doMC)

#Set the number of cores you'd like to use for parallel computing
registerDoMC(4)

#Upload the files you will be using, and merge
setwd("~/Desktop/ForGithub")
tASVs <- read_csv("tASVs.csv")
meta <- read_csv("metadata.csv")
metaASVs <- merge(tASVs, meta, by = "samples")

#Perform a bit of cleanup
rm(tASVs)
rm(meta)

#####
#####All tissues from natural Fucus spp.
#####

#Subset your data to correct set
nat <- subset(metaASVs, treatment == "none_Fv" |
  treatment == "none_Fs" |
  treatment == "none_Fd")

#Remove abundance data from metadata
nat_abund <- as.matrix(nat[, 2:6321])

#Sets a seed for reproducability
set.seed(1124)

#Runs NMDS 1000 times based on randomly selected seeds
simResults <- foreach(i = 1:1000) %dopar% {
  seed <- sample(1:99999999, 1)
  set.seed(seed)
  nmds_nat = metaMDS(nat_abund, trymax = 500, distance = "bray")
  print(c(nmds_nat$stress, seed))
}

#Stores the results in a dataframe and finds the optimum seed (lowest stress)
Results <- data.frame(matrix(nrow = 1000, ncol = 2))
names(Results) <- c("Stress", "Seed")
for (i in 1:nrow(Results)) {
  Results[i, 1] <- simResults[[i]][1]
  Results[i, 2] <- simResults[[i]][2]
}

optSeed <- Results[which.min(Results$Stress), 2]
#19534983 was my optimum seed

#Sets seed to optimum, and recreates NMDS
set.seed(optSeed)
nmds_nat = metaMDS(nat_abund, trymax = 500, distance = "bray")

#Pulls out data so it can be used in a ggplot figure
data.scores_nat = as.data.frame(scores(nmds_nat))
data.scores_nat$Species <- nat$species
data.scores_nat$Tissue <- nat$tissue

#Creates your ggplot figure
plot_nat = ggplot(data.scores_nat, aes(x = NMDS1, y = NMDS2)) +
          geom_point(stroke = 3, size = 4, aes(shape = Tissue, colour = Species)) +
          theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"),
          axis.text.x = element_text(colour = "black", face = "bold", size = 12),
          legend.text = element_text(size = 12, face = "bold", colour = "black"),
          legend.position = "right", axis.title.y = element_text(face = "bold", size = 14),
          axis.title.x = element_text(face = "bold", size = 14, colour = "black"),
          legend.title = element_text(size = 14, colour = "black", face = "bold"),
          panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
          legend.key = element_blank()) +
          labs(x = "NMDS1", colour = "Time", y = "NMDS2", shape = "Type") +
          scale_colour_manual(name = "Species", values = c("green", "blue", "black")) +
          scale_shape_manual(name = "Tissue", values = c(19, 17, 3))

#Plots and saves it
plot_nat
ggsave(width = 6, height = 6, units = "in", "nat.svg")


#####
#####Transplant treatment
#####
tp <- subset(metaASVs, treatment != "none" & treatment != "none_Fd")
tp <- subset(tp, day == "4")


tp_abund <- as.matrix(tp[, 2:6321])

set.seed(1124)
simResults <- foreach(i = 1:1000) %dopar% {
  seed <- sample(1:99999999, 1)
  set.seed(seed)
  nmds_tp = metaMDS(tp_abund, trymax = 500, distance = "bray")
  print(c(nmds_tp$stress, seed)) }

Results <- data.frame(matrix(nrow = 1000, ncol = 2))
names(Results) <- c("Stress", "Seed")
for (i in 1:nrow(Results)) {
  Results[i, 1] <- simResults[[i]][1]
  Results[i, 2] <- simResults[[i]][2]

}
optSeed <- Results[which.min(Results$Stress), 2]
#23627768


set.seed(optSeed)
nmds_tp = metaMDS(tp_abund, trymax = 500, distance = "bray")
#plot(nmds)

data.scores_tp = as.data.frame(scores(nmds_tp))
data.scores_tp$Species <- tp$treatment
data.scores_tp$Tissue <- tp$tissue

plot_tp = ggplot(data.scores_tp, aes(x = NMDS1, y = NMDS2)) +

  geom_point(stroke = 3, size = 4, aes(shape = Tissue, colour = Species)) +
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"),
        axis.text.x = element_text(colour = "black", face = "bold", size = 12),
        legend.text = element_text(size = 12, face = "bold", colour = "black"),
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14),
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"),
        legend.title = element_text(size = 14, colour = "black", face = "bold"),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key = element_blank()) +
  labs(x = "NMDS1", colour = "Time", y = "NMDS2", shape = "Type")

plot_tp
ggsave(width = 6, height = 6, units = "in", "tp_hf.svg")
