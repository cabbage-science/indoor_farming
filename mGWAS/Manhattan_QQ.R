###########################################
#### LOADING PACKAGES AND SOURCE CODES ####
###########################################

library(qqman)
library(data.table)
library(readxl)
library(tidyverse)
library(stringr)

#######################
#### LOADING FILES ####
#######################

# Loading the significantly correlated metabolites - IDs and correlation coefficients
metabolite_corr <- read.csv("../data/significant_metabolite_correlation.csv")
# Extracting the metabolite IDs
metabolite_ids <- metabolite_corr$rownames.VIP_values_large.


# Defining parameters for manhattan plotting
GWAS_Bonn_corr_threshold <- -log10(0.05 / nrow(GWAS_matrix)) ; Suggested_threshold <- 5

# Change number everytime for each metabolite
metabo <- metabolite_ids[5]
#filename <- paste0("GAPIT_output_data/GLM_log2/GAPIT.Association.GWAS_Results.MLM.",metabolite_ids[1],".csv")
filename <- paste0("GAPIT_output_data/GLM_log2/GAPIT.Association.GWAS_Results.GLM.",metabo,".csv")
GWAS_matrix <- fread(filename,header = F, sep = ",") %>%
  mutate(V1 = paste(V1,V2,sep=",")) %>%
  select(1,3:5) %>%
  rename("SNP" = "V1",
         "Chr" = "V3",
         "Pos" = "V4",
         "P_value" = "V5")

######################################
#### MANHATTAN PLOT AND QQPLOTING ####
######################################

# Manhatton plot (900px x 600px)
# _GLM_log2_MP
# _MLM_log2_MP
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
  ylim=c(0,9),
  main = paste0("GLM for Kale metabolite ID ",metabo ," (Group 1)"),
  cex.main = 1.6, cex.axis = 1.2, cex.lab = 1.3, las = 1)

# QQplot (600px x 600px)
# _GLM_log2_QQ
# _MLM_log2_QQ
QQ_plot <- qq(GWAS_matrix$P_value, main="", col = "blue", cex = 2, pch = 1,
              cex.main = 1.8, cex.axis = 1.2, cex.lab = 1.3, las = 1)
