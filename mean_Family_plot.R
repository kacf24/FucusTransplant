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

#Coverts data into a long format for easy subsetting
data_long <- gather(ASVs, Site, measurement, colnames(ASVs[, 2:227]), factor_key = TRUE)
data_long_met <- merge(data_long, meta, by.y = "samples", by.x = "Site")
data_long_met_tax <- merge(data_long_met, tax, by.x = "samples", by.y = "ASVid")

#Subset to correct data, in this case removes water column samples,
s_data_long_met_tax <- subset(data_long_met_tax, day == "4")

#Finds the number of reads each ASV accounts for and gets rid of ASVs with 0 reads
test <- do.call(data.frame, aggregate(s_data_long_met_tax$measurement, list(s_data_long_met_tax$samples), FUN = function(x) c(mn = sum(x))))
non_zero_data_long_met_tax <- subset(s_data_long_met_tax, s_data_long_met_tax$samples %notin% test[which(test$x == 0), 1])

#Finds the number of reads each ASV accounts for by tissue and treatment and removes records with 0 reads total
test <- do.call(data.frame, aggregate(non_zero_data_long_met_tax$measurement, list(non_zero_data_long_met_tax$samples, non_zero_data_long_met_tax$treatment, non_zero_data_long_met_tax$tissue), FUN = function(x) c(mn = sum(x))))
test$combo <- paste0(test$Group.1, test$Group.2, test$Group.3)
non_zero_data_long_met_tax$combo <- paste0(non_zero_data_long_met_tax$samples, non_zero_data_long_met_tax$treatment, non_zero_data_long_met_tax$tissue)
non_zero_data_long_met_tax2 <- subset(non_zero_data_long_met_tax, non_zero_data_long_met_tax$combo %notin% test[which(test$x == 0), 5])


#Finds the number of read each family accounts for by tissue and treatment
test <- do.call(data.frame, aggregate(non_zero_data_long_met_tax2$measurement, list(non_zero_data_long_met_tax2$Family, non_zero_data_long_met_tax2$treatment, non_zero_data_long_met_tax2$tissue), FUN = function(x) c(mn = sum(x))))

#Creates a metadata subset, finds the n of each treatment by tissue
day4_meta <- subset(meta, day == "4")
day4_samples <- do.call(data.frame, aggregate(day4_meta$sample_ID, list(day4_meta$treatment, day4_meta$tissue), FUN = function(x) c(mn = length(x))))
day4_samples$combo <- paste0(day4_samples$Group.1, day4_samples$Group.2)
test$combo <- paste0(test$Group.2, test$Group.3)

#merges together number of reads per family, n of samples, and calculates a mean number of reads of family/sample by treatment and tissue
test_day4_samples <- merge(test, day4_samples, by = "combo")
test_day4_samples <- test_day4_samples[, -c(1, 6, 7)]
colnames(test_day4_samples) <- c("Family", "Treatment", "Tissue", "Sum", "n")
test_day4_samples$avg <- test_day4_samples$Sum / test_day4_samples$n

#Creates a cutoff, and get rid of things below that
cutoff <- 170 / 2
test_day4_samples <- subset(test_day4_samples, test_day4_samples$avg > cutoff)


#Some housekeeping used to rename things
test_day4_samples$Tissue <- as.character(test_day4_samples$Tissue)
test_day4_samples[test_day4_samples$Tissue == "holdfast", 3] <- "H"
test_day4_samples[test_day4_samples$Tissue == "reproductive", 3] <- "R"
test_day4_samples[test_day4_samples$Tissue == "vegetative", 3] <- "V"
test_day4_samples[test_day4_samples$Tissue == "na", 3] <- "WC"

Phylum <- tax[, c(3, 6)]

test_day4_samples <- merge(test_day4_samples, Phylum, by = "Family", all.x = TRUE, all.y = FALSE)
test_day4_samples <- unique(test_day4_samples)
test_day4_samples$Family <- as.character(test_day4_samples$Family)
test_day4_samples$Phylum <- as.character(test_day4_samples$Phylum)
test_day4_samples2 <- test_day4_samples[order(test_day4_samples[, 7], test_day4_samples[, 1]),]
test_day4_samples2$Family <- as_factor(test_day4_samples2$Family)
test_day4_samples2$Family <- fct_rev(test_day4_samples2$Family)
test_day4_samples <- test_day4_samples2

test_day4_samples$Family <- as.character(test_day4_samples$Family)
test_day4_samples$Phylum <- as.character(test_day4_samples$Phylum)
test_day4_samples$Family <- gsub('_', ' ', test_day4_samples$Family)

test_day4_samples$Tissue <- as.character(test_day4_samples$Tissue)
test_day4_samples[test_day4_samples$Tissue == "holdfast", 3] <- "H"
test_day4_samples[test_day4_samples$Tissue == "reproductive", 3] <- "R"
test_day4_samples[test_day4_samples$Tissue == "vegetative", 3] <- "V"
test_day4_samples[test_day4_samples$Tissue == "na", 3] <- "WC"

test_day4_samples$Treatment <- as.character(test_day4_samples$Treatment)
test_day4_samples[test_day4_samples$Treatment == "none_Fv", 2] <- "Fv"
test_day4_samples[test_day4_samples$Treatment == "none_Fs", 2] <- "Fs"
test_day4_samples[test_day4_samples$Treatment == "none_Fd", 2] <- "Fd"
test_day4_samples[test_day4_samples$Treatment == "watered", 2] <- "W"
test_day4_samples[test_day4_samples$Treatment == "dry", 2] <- "D"
test_day4_samples[test_day4_samples$Treatment == "control", 2] <- "PC"
test_day4_samples[test_day4_samples$Treatment == "none", 2] <- "WC"
str(test_day4_samples)

#correctly reorders things
test_day4_samples$Family <- as_factor(test_day4_samples$Family)
test_day4_samples$Family <- fct_rev(test_day4_samples$Family)

#Makes and saves plot
tiff("Dotplot_updated_color.tiff", width = 174, height = 234, units = 'mm', res = 300)
ggplot(test_day4_samples, aes(x = test_day4_samples$Family, y = test_day4_samples$avg)) +
  geom_point(aes(shape = test_day4_samples$Treatment, color = test_day4_samples$Treatment), size = 3) +
  coord_flip() +
  facet_grid(col = vars(test_day4_samples$Tissue)) +
  xlab("Bacterial Family") +
  ylab("Mean Reads") +
  labs(shape = "Treatment") +
  scale_shape_manual(values = c(15, 16, 17, 3, 4, 16, 15)) +
  scale_color_manual(name = "Treatment", values = c("red", "green", "blue", "black", "darkorange4", "purple", "tan4")) +
  theme_bw() +
  theme(strip.text.x = element_text(size = 16, face = "bold"), axis.text.y = element_text(size = 8, face = "bold"), axis.text.x = element_text(size = 8, face = "bold", angle = -90), axis.title.x = element_text(size = 16, face = "bold"), axis.title.y = element_text(size = 16, face = "bold"), legend.text = element_text(size = 10, face = "bold"), legend.title = element_text(size = 12, face = "bold")) +
  scale_y_continuous(breaks = c(0, 2500, 5000, 7500, 10000), labels = c("0", "", "5000", "", "10000"))
dev.off()
