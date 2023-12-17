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
GWAS_Bonn_corr_threshold <- -log10(0.05 / 274183) ; Suggested_threshold <- 5

# Change number everytime for each metabolite
metabo <- metabolite_ids[1]
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
png("my_base_plot2.png", width = 900, height = 600)
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
  ylim=c(0,8),
  main = paste0("GLM for Kale metabolite ID ",metabo ," (Group 1)"),
  cex.main = 1.6, cex.axis = 1.2, cex.lab = 1.3, las = 1)
dev.off()

# QQplot (600px x 600px)
# _GLM_log2_QQ
# _MLM_log2_QQ
QQ_plot <- qq(GWAS_matrix$P_value, main="", col = "blue", cex = 2, pch = 1,
              cex.main = 1.8, cex.axis = 1.2, cex.lab = 1.3, las = 1)

############################
#### CALCULATING LAMBDA ####
############################

GWAS_matrix <- GWAS_matrix %>% arrange(P_value)
P_values_vec <- GWAS_matrix$P_value
inflation <- function(ps) { 
  chisq <- qchisq(1 - ps, 1) 
  lambda <- median(chisq) / qchisq(0.5, 1) 
  print(lambda)}
inflation(P_values_vec)


###################################################
#### FOR LOOP FOR SAVING MAHANTTAN AND QQPLOTS ####
###################################################

model <- "GLM"

for (k in 1:nrow(metabolite_corr)) {
  # Preparing data
  metabo <- metabolite_ids[k]
  filename <- paste0("GAPIT_output_data/",model,"_log2/GAPIT.Association.GWAS_Results.",model,".",metabo,".csv")
  GWAS_matrix <- fread(filename,header = F, sep = ",") %>%
    mutate(V1 = paste(V1,V2,sep=",")) %>%
    select(1,3:5) %>%
    rename("SNP" = "V1",
           "Chr" = "V3",
           "Pos" = "V4",
           "P_value" = "V5")
  # Plot manhattan plot
  png(paste0(metabo,"_",model,"_log2_MP.png"), width = 900, height = 600)
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
    main = paste0(model," for Kale metabolite ID ",metabo ," (Group 1)"),
    cex.main = 1.6, cex.axis = 1.2, cex.lab = 1.3, las = 1)
  dev.off()
  # Calculate lambda value
  GWAS_matrix <- GWAS_matrix %>% arrange(P_value)
  P_values_vec <- GWAS_matrix$P_value
  inflation <- function(ps) { 
    chisq <- qchisq(1 - ps, 1) 
    lambda <- median(chisq) / qchisq(0.5, 1) 
    print(lambda)}
  lambda <- round(inflation(P_values_vec),2)
  # Plot QQ plot
  png(paste0(metabo,"_",model,"_log2_QQ.png"), width = 600, height = 600)
  QQ_plot <- qq(GWAS_matrix$P_value, main="", col = "blue", cex = 2, pch = 1,
                cex.main = 1.8, cex.axis = 1.2, cex.lab = 1.3, las = 1)
  text(x = par("usr")[2] - 1.2, y = par("usr")[3] + 0.4, labels = paste0(" λ = ", lambda), 
       pos = 3, adj = 1, cex = 3)
  dev.off()
}
