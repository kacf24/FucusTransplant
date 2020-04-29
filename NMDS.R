#Load the used libraries.
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
optSeed <- 19534983
#Sets seed to optimum, and recreates NMDS
set.seed(optSeed)
nmds_nat = metaMDS(nat_abund, trymax = 500, distance = "bray")

#Pulls out data so it can be used in a ggplot figure
data.scores_nat = as.data.frame(scores(nmds_nat))
data.scores_nat$Species <- nat$species
data.scores_nat$Tissue <- nat$tissue

#Creates your ggplot figure
plot_nat = ggplot(data.scores_nat, aes(x = NMDS1, y = NMDS2)) +
  geom_point(stroke = 1, size = 2, aes(shape = Tissue, colour = Species)) +
  theme(axis.text.y = element_text(colour = "black", size = 10, face = "bold"),
        axis.text.x = element_text(colour = "black", face = "bold", size = 10),
        legend.text = element_text(size = 10, face = "bold", colour = "black"),
        legend.position = "bottom", axis.title.y = element_text(face = "bold", size = 10),
        axis.title.x = element_text(face = "bold", size = 10, colour = "black"),
        legend.title = element_text(size = 10, colour = "black", face = "bold"),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key = element_blank()) +
  labs(x = "NMDS1", colour = "Time", y = "NMDS2", shape = "Type") +
  scale_colour_manual(name = "Species", breaks = c("Fs", "Fv", "Fd"), values = c("#F0E442", "#0072B2", "#D55E00")) +
  scale_shape_manual(name = "Tissue", labels = c("H", "R", "V"), values = c(15, 17, 19)) +
  guides(col = guide_legend(title.position = "top", ncol = 1), shape = guide_legend(title.position = "top", ncol = 1))

#Plots and saves it
plot_nat
ggsave(width = 89, height = 120, units = "mm", "nat.svg")


#####
#####All tissues from transplants, Fs, and Fv
#####

#Subset your data to correct set
tp <- subset(metaASVs, treatment != "none" & treatment != "none_Fd")
tp <- subset(tp, day == "4")

#Remove abundance data from metadata
tp_abund <- as.matrix(tp[, 2:6321])

#Sets a seed for reproducability
set.seed(1124)

#Runs NMDS 1000 times based on randomly selected seeds
simResults <- foreach(i = 1:1000) %dopar% {
  seed <- sample(1:99999999, 1)
  set.seed(seed)
  nmds_tp = metaMDS(tp_abund, trymax = 500, distance = "bray")
  print(c(nmds_tp$stress, seed))
}

#Stores the results in a dataframe and finds the optimum seed (lowest stress)
Results <- data.frame(matrix(nrow = 1000, ncol = 2))
names(Results) <- c("Stress", "Seed")
for (i in 1:nrow(Results)) {
  Results[i, 1] <- simResults[[i]][1]
  Results[i, 2] <- simResults[[i]][2]

}
optSeed <- Results[which.min(Results$Stress), 2]
#23627768 was my optimum seed
optSeed <- 23627768
#Sets seed to optimum, and recreates NMDS
set.seed(optSeed)
nmds_tp = metaMDS(tp_abund, trymax = 500, distance = "bray")

#Pulls out data so it can be used in a ggplot figure
data.scores_tp = as.data.frame(scores(nmds_tp))
data.scores_tp$Species <- tp$treatment
data.scores_tp$Tissue <- tp$tissue

#Creates your ggplot figure
plot_tp = ggplot(data.scores_tp, aes(x = NMDS1, y = NMDS2)) +

  geom_point(stroke = 1, size = 2, aes(shape = Tissue, colour = Species)) +
  theme(axis.text.y = element_text(colour = "black", size = 10, face = "bold"),
        axis.text.x = element_text(colour = "black", face = "bold", size = 10),
        legend.text = element_text(size = 10, face = "bold", colour = "black"),
        legend.position = "bottom", axis.title.y = element_text(face = "bold", size = 10),
        axis.title.x = element_text(face = "bold", size = 10, colour = "black"),
        legend.title = element_text(size = 10, colour = "black", face = "bold"),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key = element_blank()) +
  labs(x = "NMDS1", colour = "Time", y = "NMDS2", shape = "Type") +
  scale_shape_manual(name = "Tissue", labels = c("H", "R", "V"), values = c(15, 17, 19)) +
  scale_colour_manual(name = "Sample Type", breaks = c("dry", "watered", "none_Fs", "control", "none_Fv"), labels = c("D", "W", "Fs", "PC", "Fv"), values = c("#000000", "#E69F00", "#F0E442", "#0072B2", "#009E73", "#0072B2")) +
  guides(col = guide_legend(title.position = "top", ncol = 1), shape = guide_legend(title.position = "top", ncol = 1))


#Plots and saves it
plot_tp
ggsave(width = 89, height = 120, units = "mm", "tp.svg")

#####
#####All samples
#####

#Subset your data to correct set
all <- metaASVs

#Remove abundance data from metadata
all_abund <- as.matrix(all[, 2:6321])

#Sets a seed for reproducability
set.seed(1124)

#Runs NMDS 1000 times based on randomly selected seeds
simResults <- foreach(i = 1:1000) %dopar% {
  seed <- sample(1:99999999, 1)
  set.seed(seed)
  nmds_all = metaMDS(all_abund, trymax = 500, distance = "bray")
  print(c(nmds_all$stress, seed))
}

#Stores the results in a dataframe and finds the optimum seed (lowest stress)
Results <- data.frame(matrix(nrow = 1000, ncol = 2))
names(Results) <- c("Stress", "Seed")
for (i in 1:nrow(Results)) {
  Results[i, 1] <- simResults[[i]][1]
  Results[i, 2] <- simResults[[i]][2]

}
optSeed <- Results[which.min(Results$Stress), 2]
#80450461 was my optimum seed

optSeed <- 80450461
#Sets seed to optimum, and recreates NMDS
set.seed(optSeed)
nmds_all = metaMDS(all_abund, trymax = 500, distance = "bray")

#Pulls out data so it can be used in a ggplot figure
data.scores_all = as.data.frame(scores(nmds_all))
data.scores_all$Species <- all$treatment
data.scores_all$Tissue <- all$tissue

#Creates your ggplot figure

plot_all = ggplot(data.scores_all, aes(x = NMDS1, y = NMDS2)) +
  geom_point(stroke = 1, size = 2, aes(colour = Species, shape = Tissue)) +
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"),
        axis.text.x = element_text(colour = "black", face = "bold", size = 12),
        legend.text = element_text(size = 12, face = "bold", colour = "black"),
        legend.position = "bottom", axis.title.y = element_text(face = "bold", size = 14),
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"),
        legend.title = element_text(size = 14, colour = "black", face = "bold"),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key = element_blank()) +
  xlim(c(-3, 2)) +
  ylim(c(-3, 3)) +
  labs(x = "NMDS1", colour = "Sample Type", y = "NMDS2", shape = "Tissue") +
  scale_shape_manual(name = "Tissue", breaks = c("holdfast", "reproductive", "vegetative", "na"), labels = c("H", "R", "V", "WC"), values = c(15, 17, 19, 5)) +
  scale_colour_manual(name = "Sample Type", breaks = c("dry", "watered", "none_Fs", "control", "none_Fv", "none_Fd", "none"), labels = c("D", "W", "Fs", "PC", "Fv", "Fd", "WC"), values = c("#000000", "#E69F00", "#56B4E9", "#D55E00", "#F0E442", "#0072B2", "#009E73", "#CC79A7")) +
  guides(col = guide_legend(title.position = "top", ncol = 1), shape = guide_legend(title.position = "top", ncol = 1))


#Plots and saves it
plot_all
ggsave(width = 89, units = "mm", "allSamples.svg")

#####
#####Vegetative Transplant Tissue
#####

#Subset your data to correct set
veg_tp <- subset(metaASVs, treatment != "none_Fv" &
  treatment != "none_Fs" &
  treatment != "none_Fd" &
  tissue == "vegetative")

#Remove abundance data from metadata
veg_tp_abund <- as.matrix(veg_tp[, 2:6321])

#Sets a seed for reproducability
set.seed(1124)

#Runs NMDS 1000 times based on randomly selected seeds
simResults <- foreach(i = 1:1000) %dopar% {
  seed <- sample(1:99999999, 1)
  set.seed(seed)
  nmds_veg_tp = metaMDS(veg_tp_abund, trymax = 500, distance = "bray")
  print(c(nmds_all$stress, seed))
}

#Stores the results in a dataframe and finds the optimum seed (lowest stress)
Results <- data.frame(matrix(nrow = 1000, ncol = 2))
names(Results) <- c("Stress", "Seed")
for (i in 1:nrow(Results)) {
  Results[i, 1] <- simResults[[i]][1]
  Results[i, 2] <- simResults[[i]][2]

}
optSeed <- Results[which.min(Results$Stress), 2]
#61562574 was my optimum seed

optSeed <- 61562574
#Sets seed to optimum, and recreates NMDS
set.seed(optSeed)
nmds_veg_tp = metaMDS(veg_tp_abund, trymax = 500, distance = "bray")

#Pulls out data so it can be used in a ggplot figure
data.scores_veg_tp = as.data.frame(scores(nmds_veg_tp))
data.scores_veg_tp$Species <- veg_tp$treatment
data.scores_veg_tp$Day <- as.character(veg_tp$day)

#Creates your ggplot figure
plot_veg_tp = ggplot(data.scores_veg_tp, aes(x = NMDS1, y = NMDS2)) +

  geom_point(stroke = 1, size = 2, aes(colour = Species, shape = Day)) +
  theme(axis.text.y = element_text(colour = "black", size = 10, face = "bold"),
        axis.text.x = element_text(colour = "black", face = "bold", size = 10),
        legend.text = element_text(size = 10, face = "bold", colour = "black"),
        legend.position = "bottom", axis.title.y = element_text(face = "bold", size = 10),
        axis.title.x = element_text(face = "bold", size = 10, colour = "black"),
        legend.title = element_text(size = 10, colour = "black", face = "bold"),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key = element_blank()) +
  labs(x = "NMDS1", colour = "Sample Type", y = "NMDS2", shape = "Tissue") +
  scale_shape_manual(name = "Day", labels = c("July 6", "July 11", "July 20"), values = c(3, 4, 6)) +
  scale_colour_manual(name = "Treatment", breaks = c("dry", "watered", "control"), labels = c("D", "W", "PC"), values = c("#E69F00", "#009E73", "#000000")) +
  guides(col = guide_legend(title.position = "top", ncol = 1), shape = guide_legend(title.position = "top", ncol = 1))


#Plots and saves it
plot_veg_tp
ggsave(width = 89, height = 120, units = "mm", "vegTissueTP.svg")

