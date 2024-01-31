###########################################
#### LOADING PACKAGES AND SOURCE CODES ####
###########################################

library(data.table)
library(readxl)
source("http://zzlab.net/GAPIT/gapit_functions.txt")
library(tidyverse)
library(qqman)
library(stringr)

#########################################
#### LOADING FILES AND CLEANING DATA ####
#########################################

# Loading genotype map
myGM <- read.table("../data/kale_LCMS_ABTS_371_genotype_map.txt", header = TRUE)
# Creating a vector for myGD column renaming
SNP <- myGM$SNP

# Loading numerical data (myGD)
myGD <- data.frame(fread("../data/kale_LCMS_ABTS_371_numeric.txt", header = TRUE)) %>%
  rename("sample_ID" = "taxa") %>%
  mutate(sample_ID = str_extract(sample_ID, "(?<=BAM-sorted\\/)[^\\/]+(?=-sorted.bam)")) %>%
  mutate(sample_ID = gsub("\\-",":",toupper(sample_ID)))
colnames(myGD) <- c("sample_ID",SNP)

# Loading significant metabolite IDs identified upstream
met_ids_csv <- read.csv("../data/significant_metabolite_info.csv", header = TRUE)
# This vector is important for the loop later on
met_ids <- met_ids_csv$metabolite_ID
# Loading LCMS peak intensity data (phenotype - myY)
sig_metabolite_intensities <- read.csv("../data/kale_peaks_log2_noQC.csv", header = TRUE) %>%
  select(6, all_of(met_ids)) %>%
  # Rename column
  rename(sample_ID = OriginalID) %>%
  # Configuring the sample IDs
  mutate(sample_ID = gsub("\\.",":",toupper(sample_ID))) %>%
  # Filtering metabolite data for samples only in common with genotype data (n = 371 with some repeats)
  # Repeats removed in the next command
  filter(sample_ID %in% myGD$sample_ID)

# Removing of repeated samples
sig_metabolite_intensities <- sig_metabolite_intensities[which(!duplicated(sig_metabolite_intensities$sample_ID)),]
#
myY <- sig_metabolite_intensities %>% select(1,5)

#
faststructure <- data.frame(fread("../data/finalimputedout.7.meanQ")) 
sample_name <- data.frame(fread("../data/finalimputed_prune_0.5.fam")) %>% 
  select(1) %>%
  rename("sample_ID" = "V1") %>%
  mutate(sample_ID = str_extract(sample_ID, "(?<=BAM-sorted\\/)[^\\/]+(?=-sorted.bam)")) %>%
  mutate(sample_ID = gsub("\\-",":",toupper(sample_ID)))
#
myCV <- cbind(sample_name, faststructure) %>%
  filter(sample_ID %in% myY$sample_ID)

#
myKI <- read.table(paste0("../data/kale_LCMS_ABTS_371_kinship_centered.txt"), header = FALSE) %>%
  mutate(V1 = str_extract(V1, "(?<=BAM-sorted\\/)[^\\/]+(?=-sorted.bam)")) %>%
  mutate(V1 = gsub("\\-",":",toupper(V1)))

myGAPIT <- GAPIT(
  Y=myY,
  GD=myGD,
  GM=myGM,
  CV=myCV,
  KI=myKI,
  model="GLM")




GWAS_Bonn_corr_threshold <- -log10(0.05 / 274183) ; Suggested_threshold <- 5

filename <- "GAPIT.Association.GWAS_Results.GLM.align_772.csv"
GWAS_matrix <- fread(filename,header = F, sep = ",") %>%
  mutate(V1 = paste(V1,V2,sep=",")) %>%
  select(1,3:5) %>%
  rename("SNP" = "V1",
         "Chr" = "V3",
         "Pos" = "V4",
         "P_value" = "V5")

png("fast_structure_testing_MP.png", width = 1000, height = 600)
Mann_plot <- manhattan(
  GWAS_matrix,
  chr = "Chr",
  bp = "Pos",
  snp = "SNP",
  p = "P_value",
  col = c("blue3","goldenrod","blue3","goldenrod","blue3","goldenrod","blue3","goldenrod","blue3"),
  annotateTop = T,
  genomewideline = GWAS_Bonn_corr_threshold,
  suggestiveline = Suggested_threshold,
  main = "Testing for fast structure - align_772",
  cex.main = 1.6, cex.axis = 1.2, cex.lab = 1.3, las = 1)
dev.off()

# Calculate lambda value
GWAS_matrix <- GWAS_matrix %>% arrange(P_value)
P_values_vec <- GWAS_matrix$P_value
inflation <- function(ps) { 
  chisq <- qchisq(1 - ps, 1) 
  lambda <- median(chisq) / qchisq(0.5, 1)}
lambda <- round(inflation(P_values_vec),2)
# Plot QQ plot
png("fast_structure_testing_QQ.png", width = 600, height = 600)
QQ_plot <- qq(GWAS_matrix$P_value, main="", col = "blue", cex = 2, pch = 1,
              cex.main = 1.8, cex.axis = 1.2, cex.lab = 1.3, las = 1)
text(x = par("usr")[2] - 1.2, y = par("usr")[3] + 0.4, labels = paste0(" λ = ", lambda), 
     pos = 3, adj = 1, cex = 3)
dev.off()