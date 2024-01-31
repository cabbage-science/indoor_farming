###########################################
#### LOADING PACKAGES AND SOURCE CODES ####
###########################################

library(data.table)
library(readxl)
source("http://zzlab.net/GAPIT/gapit_functions.txt")
library(tidyverse)
library(qqman)
library(stringr)

my_error_handler <- function(e) {
  message("An error occurred: ", e$message)
  # Return to top-level (console) without halting execution
  return(invisible(NULL))
}

# Set this function as the global error handler
options(error = my_error_handler)

#########################################
#### LOADING FILES AND CLEANING DATA ####
#########################################
# Loading of covariance and kinship data files would be performed in the loop later.

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

###########################################################################
#### mGWAS using GLM for all samples, for 3 different covariance files ####
###########################################################################
no_cv <- paste0(c(3,5,10),"PC")
#no_cv <- "3PC"
kinship_type <- c("centered", "dominance_centered", "dominance_normalized", "normalized")
#kinship_type <- "centered"
GWAS_Bonn_corr_threshold <- -log10(0.05 / 274183) ; Suggested_threshold <- 5

# Creating empty data frame to collect lambda values
lambda_df <- data.frame(NULL)

for (j in 2:ncol(sig_metabolite_intensities)) {
#for (j in 5:5) {
  myY <- sig_metabolite_intensities[,c(1,j)]
  for (filename_cv in no_cv) {
    #
    myCV <- read.table(paste0("../data/kale_LCMS_ABTS_371_covariates_",filename_cv,".txt"), header = TRUE) %>%
      rename("sample_ID" = "Taxa") %>%
      mutate(sample_ID = str_extract(sample_ID, "(?<=BAM-sorted\\/)[^\\/]+(?=-sorted.bam)")) %>%
      mutate(sample_ID = gsub("\\-",":",toupper(sample_ID)))
    
    # Setting WD to results to not flood src directory
    setwd("../results")
    
    # Running GAPIT function for GLM
    try({
      myGAPIT <- GAPIT(
        Y=myY,
        GD=myGD,
        GM=myGM,
        CV=myCV,
        model="GLM")
    })
    # Loading association analysis SNP files from results directory for plotting
    # And also lambda value calculations
    filename <- paste0("GAPIT.Association.GWAS_Results.GLM.",met_ids[j-1],".csv")
    GWAS_matrix <- fread(filename,header = F, sep = ",") %>%
      mutate(V1 = paste(V1,V2,sep=",")) %>%
      select(1,3:5) %>%
      rename("SNP" = "V1",
             "Chr" = "V3",
             "Pos" = "V4",
             "P_value" = "V5")
    
    # Plot manhattan plot
    png(paste0("../images/",met_ids[j-1],"_GLM_",filename_cv,"_log2_MP.png"), width = 1000, height = 600)
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
      main = paste0("GLM with ",filename_cv,"s for metabolite ID ",met_ids[j-1]," (Kale)"),
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
    png(paste0("../images/",met_ids[j-1],"_GLM_",filename_cv,"_log2_QQ.png"), width = 600, height = 600)
    QQ_plot <- qq(GWAS_matrix$P_value, main="", col = "blue", cex = 2, pch = 1,
                  cex.main = 1.8, cex.axis = 1.2, cex.lab = 1.3, las = 1)
    text(x = par("usr")[2] - 1.2, y = par("usr")[3] + 0.4, labels = paste0(" λ = ", lambda), 
         pos = 3, adj = 1, cex = 3)
    dev.off()
    
    # Appending to data frame to collect lambda values
    temp_vector <- c(met_ids[j-1], filename_cv, NA, lambda)
    lambda_df <- rbind(lambda_df, temp_vector)
    
    # Creating new GWAS result files since i cant be bothered to customize the GWAS output
    # data file in the GAPIT source code
    filename_temp <- paste0("GAPIT.Association.GWAS_Results.GLM.",met_ids[j-1],".csv")
    temp_GWAS_results <- fread(filename_temp, header = F, sep = ",")
    write.csv(temp_GWAS_results, paste0("GLM_",met_ids[j-1],"_",filename_cv,"_results.csv"), row.names = FALSE)
    
    # Reiterate for each of the 4 kinship matrices for each no. of PCs
    for (kinship in kinship_type) {
      #
      myCV <- read.table(paste0("../data/kale_LCMS_ABTS_371_covariates_",filename_cv,".txt"), header = TRUE) %>%
        rename("sample_ID" = "Taxa") %>%
        mutate(sample_ID = str_extract(sample_ID, "(?<=BAM-sorted\\/)[^\\/]+(?=-sorted.bam)")) %>%
        mutate(sample_ID = gsub("\\-",":",toupper(sample_ID)))
      
      #
      myKI <- read.table(paste0("../data/kale_LCMS_ABTS_371_kinship_",kinship,".txt"), header = FALSE) %>%
        mutate(V1 = str_extract(V1, "(?<=BAM-sorted\\/)[^\\/]+(?=-sorted.bam)")) %>%
        mutate(V1 = gsub("\\-",":",toupper(V1)))
      
      # Running GAPIT function for GLM
      try({
        myGAPIT <- GAPIT(
          Y=myY,
          GD=myGD,
          GM=myGM,
          CV=myCV,
          KI=myKI,
          model="MLM")
      })
      # Loading association analysis SNP files from results directory for plotting
      # And also lambda value calculations
      filename <- paste0("GAPIT.Association.GWAS_Results.MLM.",met_ids[j-1],".csv")
      GWAS_matrix <- fread(filename,header = F, sep = ",") %>%
        mutate(V1 = paste(V1,V2,sep=",")) %>%
        select(1,3:5) %>%
        rename("SNP" = "V1",
               "Chr" = "V3",
               "Pos" = "V4",
               "P_value" = "V5")
      
      # Plot manhattan plot
      png(paste0("../images/",met_ids[j-1],"_MLM_",filename_cv,"_",kinship,"_log2_MP.png"), width = 1000, height = 600)
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
        main = paste0("MLM with ",filename_cv,"s for metabolite ID ",met_ids[j-1]," (Kale)"),
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
      png(paste0("../images/",met_ids[j-1],"_MLM_",filename_cv,"_",kinship,"_log2_QQ.png"), width = 600, height = 600)
      QQ_plot <- qq(GWAS_matrix$P_value, main="", col = "blue", cex = 2, pch = 1,
                    cex.main = 1.8, cex.axis = 1.2, cex.lab = 1.3, las = 1)
      text(x = par("usr")[2] - 1.2, y = par("usr")[3] + 0.4, labels = paste0(" λ = ", lambda), 
           pos = 3, adj = 1, cex = 3)
      dev.off()
      
      # Appending to data frame to collect lambda values
      temp_vector <- c(met_ids[j-1], filename_cv, kinship, lambda)
      lambda_df <- rbind(lambda_df, temp_vector)
      
      # Creating new GWAS result files since i cant be bothered to customize the GWAS output
      # data file in the GAPIT source code
      filename_temp <- paste0("GAPIT.Association.GWAS_Results.MLM.",met_ids[j-1],".csv")
      temp_GWAS_results <- fread(filename_temp, header = F, sep = ",")
      write.csv(temp_GWAS_results, paste0("MLM_",met_ids[j-1],"_",filename_cv,"_",kinship,"_results.csv"), row.names = FALSE)
    }
  }
}

# Giving column names to lambda df
colnames(lambda_df) <- c("metabolite_ID","no_PCs","kinship_type","lambda_value")

write.table(lambda_df, "../src/lambda_value_for_each_iteration.txt", quote = FALSE, row.names = FALSE, sep = "\t")
