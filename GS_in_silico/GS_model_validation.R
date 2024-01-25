##########################
#### LOADING PACKAGES ####
##########################

#install.packages("BGLR")
#install.packages("BWGS")
library(BWGS)
library(rrBLUP)
library(data.table)
library(tidyverse)

######################
#### LOADING DATA ####
######################

myY_raw <- read.table("../data/phenotype_grp1.txt", header = TRUE) %>%
  mutate(sample_ID = str_extract(Taxa, "(?<=BAM-sorted\\/)[^\\/]+(?=-sorted.bam)")) %>%
  select(sample_ID, everything()) %>% select(-2)
 
myGM_all <- data.frame(fread("../data/GAPIT_GWAS_SUPER_AC.csv")) %>%
#myGM_all <- data.frame(fread("../data/GAPIT_GWAS_MLM_AC.csv")) %>%
#myGM_all <- data.frame(fread("../data/GAPIT_GWAS_GLM_AC.csv")) %>%
  unite(col = SNP_ID, c("SNP","SNP2"), sep = ",") %>%
  select(SNP_ID, Chr, Pos, P.value) %>%
  mutate(SNP_ID = str_replace(SNP_ID, ":",".")) %>%
  mutate(SNP_ID = str_replace(SNP_ID, ",","."))

# Sample_ID column does not provide information on sample ID here, because the numerical .txt file does not have it.
# Hence this script only works for such data. To avoid being unsure in future, please
# include sample names in all files - genotype and phenotype files.
myGD_raw <- data.frame(fread("../data/final_grp1_numeric.txt", col.names = c("sample_ID",myGM_all$SNP_ID)))

######################################################
#### SETTING PARAMETERS AND CREATING NEW MATRICES ####
######################################################

### SETTING PARAMETERS
# Specify GWAS model that produced association testing data (important for optimal SNP selection)
GWAS_model <- "SUPER"
# Setting arbitrary number of SNPs for prediction
# MLM = 2800, GLM = 2100, SUPER = 1400
no_snps <- 1400
# Setting model for genomic selection (GBLUP, LASSO, BA, BB, BRR)
gs_method <- "GBLUP"
# Setting proportion of samples to be used for validation
proportion_validation <- 0.1

### PREPARING NEW MATRICES
# Creating a vector with sample IDs as row names
myY <- setNames(myY_raw$Antioxidant_capacity, myY_raw$sample_ID)

# Creating a genotype map data frame of selected SNPs
myGM <- myGM_all %>%
  arrange(P.value) %>%
  head(no_snps)

# Creating a matrix of numerical genotype data with sample IDs as row names
myGD <- as.matrix(myGD_raw[,myGM$SNP_ID])
row.names(myGD) <- myY_raw$sample_ID

###########################################
#### RUNNING THE BGWS.PREDICT FUNCTION ####
###########################################

#set.seed(123) 

# Get the number of rows (samples) in genotype data 
total_rows <- nrow(myGD)
# Calculate the number of samples for validation (not used for model building)
validation_rows <- round(proportion_validation * total_rows)
# Create vector of random numbers to select samples for validation
validation_indices <- sample(1:total_rows, validation_rows)

# Create training and validation sets for both genotype and phenotype
geno_train <- myGD[-validation_indices,] ; geno_valid <- myGD[validation_indices,]
pheno_train <- myY[-validation_indices] ; pheno_valid <- myY[validation_indices]

# Perform model building on training data set, and use to predict GEBV/phenotype
# based on validation set's genotype data
predicted_data <- bwgs.predict (geno_train, pheno_train, geno_valid,
                                geno.reduct.method="NULL",reduct.size="NULL", pval="NULL", r2="NULL", MAP="NULL",
                                MAXNA = 1, MAF = 0,
                                predict.method = gs_method)

# Creating a dataframe (required for ggplot) of actual and predicted phenotype values
pheno_df <- data.frame(pheno_valid = pheno_valid, pheno_pred = predicted_data[,1])

# GGPLOT
ggplot(pheno_df, aes(x = pheno_valid, y = pheno_pred)) +
  geom_smooth(method = "lm", se = TRUE, color = "red", lty = "dashed") +
  geom_point(shape = 21, size = 3) +
  geom_text(aes(label = paste0("R = ", round(cor(pheno_df$pheno_pred, pheno_df$pheno_valid), 3))),
            x = max(pheno_df$pheno_valid), y = min(pheno_df$pheno_pred),
            hjust = 1, vjust = 0, color = "black",size = 6) +
  labs(title = paste0("Validation of ",gs_method," GS model (",no_snps," SNPs) using ",
                      proportion_validation*100,"% of samples for validation"),
       x = "Actual Antioxidant Capacity (µmol TE/g)",
       y = "Predicted Antioxidant Capacity (µmol TE/g)",
       caption = paste0("SNP results from the ", GWAS_model, " model were used.")) +
  theme_minimal() +
  theme(plot.title = element_text(color="black", size=13, face="bold"),
        axis.title.x = element_text(color="black", size=11, margin = margin(t = 10)),
        axis.title.y = element_text(color="black", size=11, margin = margin(r = 10)))

############################################################
#### CALCULATE CORRELATION COEFFICIENTS FOR DIFFERENT ######
#### PROPORTION OF SAMPLES USED FOR PREDICTION/TRAINING ####
############################################################

corr_vector <- NULL

# Loop to test for 10% to 90% of samples used for prediction, in increments of 10%
for(k in 1:9) {
  proportion_validation <- k/10
  total_rows <- nrow(myGD)
  # Calculate the number of samples for validation (not used for model building)
  validation_rows <- round(proportion_validation * total_rows)
  # Create vector of random numbers to select samples for validation
  validation_indices <- sample(1:total_rows, validation_rows)
  
  # Create training and validation sets for both genotype and phenotype
  geno_train <- myGD[-validation_indices,] ; geno_valid <- myGD[validation_indices,]
  pheno_train <- myY[-validation_indices] ; pheno_valid <- myY[validation_indices]
  
  # Perform model building on training data set, and use to predict GEBV/phenotype
  # based on validation set's genotype data
  predicted_data <- bwgs.predict (geno_train, pheno_train, geno_valid,
                                  geno.reduct.method="NULL",reduct.size="NULL", pval="NULL", r2="NULL", MAP="NULL",
                                  MAXNA = 1, MAF = 0,
                                  predict.method = gs_method)
  
  # Calculating correlation coefficient
  corr_coeff <- round(cor(predicted_data[,1], pheno_valid), 3)
  corr_vector <- c(corr_vector, corr_coeff)
}

###########################################################
#### CALCULATE AVERAGE CORRELATION COEFFICIENTS FOR #######
#### DIFFERENT PROPORTION OF SAMPLES USED FOR N ROUNDS ####
#### OF PREDICTION/TRAINING PER PROPORTION OF SAMPLES #####
###########################################################

# Arbitrary number of rounds you want
no_rounds <- 100
# Creating an empty data frame for appending of values
corr_df <- data.frame(matrix(0,no_rounds,9))
colnames(corr_df) <- paste0(10*seq(1:9),"%")

# Loop to test for 10% to 90% of samples used for prediction, in increments of 10%, for n rounds
for(k in 1:9) {
  for(j in 1:no_rounds) {
    proportion_validation <- k/10
    total_rows <- nrow(myGD)
    # Calculate the number of samples for validation (not used for model building)
    validation_rows <- round(proportion_validation * total_rows)
    # Create vector of random numbers to select samples for validation
    validation_indices <- sample(1:total_rows, validation_rows)
    
    # Create training and validation sets for both genotype and phenotype
    geno_train <- myGD[-validation_indices,] ; geno_valid <- myGD[validation_indices,]
    pheno_train <- myY[-validation_indices] ; pheno_valid <- myY[validation_indices]
    
    # Perform model building on training data set, and use to predict GEBV/phenotype
    # based on validation set's genotype data
    suppressMessages({
      invisible(capture.output({
        predicted_data <- bwgs.predict (geno_train, pheno_train, geno_valid,
                                        geno.reduct.method="NULL",reduct.size="NULL", pval="NULL", r2="NULL", MAP="NULL",
                                        MAXNA = 1, MAF = 0,
                                        predict.method = gs_method)
      }))
    })
    # Calculating correlation coefficient
    corr_coeff <- round(cor(predicted_data[,1], pheno_valid), 3)
    corr_df[j,k] <- corr_coeff
    if(j %% 20 == 0) {
      print(paste0("Validation for ",j," rounds completed!"))
    }
  }
  print(paste0("Validation for ",k*10,"% of samples used for validation for ",gs_method," method completed!"))
}

new_corr_df <- pivot_longer(corr_df, `10%`:`90%`, names_to = "percentage_validation", values_to = "correlation_coefficient")

ggplot(new_corr_df, aes(x = percentage_validation, y = correlation_coefficient)) +
  stat_boxplot(geom = "errorbar",
               width = 0.2) +
  geom_boxplot(outlier.shape = 1, outlier.size = 3, size = 0.7) +
  labs(title = paste0("Prediction accuracy for GS validation (",gs_method,")"), 
       caption = paste0("For each percentage of samples used, ",no_rounds," rounds of validation were performed"),
       x = "Percentage of samples used for validation", 
       y = "Correlation coefficient") +
  theme_minimal() +
  theme(plot.title = element_text(color="black", size=14, face="bold"),
        axis.title.x = element_text(color="black", size=12, margin = margin(t = 16)),
        axis.title.y = element_text(color="black", size=12, margin = margin(r = 10)))

###########################################################
#### CALCULATE AVERAGE CORRELATION COEFFICIENTS FOR #######
#### DIFFERENT PROPORTION OF SAMPLES USED FOR N ROUNDS ####
#### OF PREDICTION/TRAINING PER PROPORTION OF SAMPLES #####
#### FOR 3 DIFFERENT MODELS ###############################
###########################################################

predict_models <- c("GBLUP", "LASSO", "BRR")

# Arbitrary number of rounds you want (per percentage of samples per GS method)
no_rounds <- 20
# Creating an empty data frame for appending of values
average_corr_df <- data.frame(matrix(0,9,1+length(predict_models)))
colnames(average_corr_df) <- c("percentage_samples_validation",predict_models)
average_corr_df$percentage_samples_validation <- paste0(10*seq(1:9),"%")

# Loop to test for 10% to 90% of samples used for prediction, in increments of 10%, for n rounds
for(gs_method in predict_models) {
  for(k in 1:9) {
    corr_vector <- NULL
    for(j in 1:no_rounds) {
      proportion_validation <- k/10
      total_rows <- nrow(myGD)
      # Calculate the number of samples for validation (not used for model building)
      validation_rows <- round(proportion_validation * total_rows)
      # Create vector of random numbers to select samples for validation
      validation_indices <- sample(1:total_rows, validation_rows)
      
      # Create training and validation sets for both genotype and phenotype
      geno_train <- myGD[-validation_indices,] ; geno_valid <- myGD[validation_indices,]
      pheno_train <- myY[-validation_indices] ; pheno_valid <- myY[validation_indices]
      
      # Perform model building on training data set, and use to predict GEBV/phenotype
      # based on validation set's genotype data
      suppressMessages({
        invisible(capture.output({
          predicted_data <- bwgs.predict (geno_train, pheno_train, geno_valid,
                                          geno.reduct.method="NULL",reduct.size="NULL", pval="NULL", r2="NULL", MAP="NULL",
                                          MAXNA = 1, MAF = 0,
                                          predict.method = gs_method)
        }))
      })
      # Calculating correlation coefficient
      corr_coeff <- round(cor(predicted_data[,1], pheno_valid), 8)
      corr_vector <- c(corr_vector,corr_coeff)
      if(j %% 4 == 0) {
        print(paste0("Validation for ",j," rounds completed!"))
      }
    }
    average_corr_df[k,gs_method] <- round(mean(corr_vector, na.rm = TRUE),3)
    print(paste0("Validation for ",k*10,"% of samples used for validation for ",gs_method," method completed!"))
  }
}

write.table(average_corr_df, paste0("average_corr_different_percentage_and_models_",GWAS_model,".txt"), quote = FALSE, row.names = FALSE, sep = "\t")

average_corr_df_final <- pivot_longer(average_corr_df, "GBLUP":"BRR", names_to = "name_of_model",
                                      values_to = "correlation_coefficient")

ggplot(average_corr_df_final, aes(x = percentage_samples_validation, y = correlation_coefficient)) +
  geom_point(aes(color = name_of_model)) +
  geom_line(aes(group = name_of_model, color = name_of_model)) +
  labs(title = paste0("Comparison of prediction accuracy between ",length(predict_models)," methods for internal validation"), 
       caption = paste0("For each percentage of samples used, ",no_rounds," rounds of validation were performed"),
       x = "Percentage of samples used for validation", 
       y = "Correlation coefficient",
       color = "Name of Model") +
  theme_minimal() +
  theme(plot.title = element_text(color="black", size=14, face="bold"),
        axis.title.x = element_text(color="black", size=12, margin = margin(t = 16)),
        axis.title.y = element_text(color="black", size=12, margin = margin(r = 10)))


###########################################################
#### CALCULATE AVERAGE CORRELATION COEFFICIENTS FOR #######
#### DIFFERENT PROPORTION OF SAMPLES USED FOR N ROUNDS ####
#### OF PREDICTION/TRAINING PER PROPORTION OF SAMPLES #####
#### FOR DIFFERENT GWAS MODELS ############################
###########################################################

# Arbitrary number of rounds you want
no_rounds <- 100
# Creating an empty data frame for appending of values
corr_df <- data.frame(matrix(0,no_rounds,9))
colnames(corr_df) <- paste0(10*seq(1:9),"%")

# Loop to test for 10% to 90% of samples used for prediction, in increments of 10%, for n rounds
for(k in 1:9) {
  for(j in 1:no_rounds) {
    proportion_validation <- k/10
    total_rows <- nrow(myGD)
    # Calculate the number of samples for validation (not used for model building)
    validation_rows <- round(proportion_validation * total_rows)
    # Create vector of random numbers to select samples for validation
    validation_indices <- sample(1:total_rows, validation_rows)
    
    # Create training and validation sets for both genotype and phenotype
    geno_train <- myGD[-validation_indices,] ; geno_valid <- myGD[validation_indices,]
    pheno_train <- myY[-validation_indices] ; pheno_valid <- myY[validation_indices]
    
    # Perform model building on training data set, and use to predict GEBV/phenotype
    # based on validation set's genotype data
    suppressMessages({
      invisible(capture.output({
        predicted_data <- bwgs.predict (geno_train, pheno_train, geno_valid,
                                        geno.reduct.method="NULL",reduct.size="NULL", pval="NULL", r2="NULL", MAP="NULL",
                                        MAXNA = 1, MAF = 0,
                                        predict.method = gs_method)
      }))
    })
    # Calculating correlation coefficient
    corr_coeff <- round(cor(predicted_data[,1], pheno_valid), 3)
    corr_df[j,k] <- corr_coeff
    if(j %% 20 == 0) {
      print(paste0("Validation for ",j," rounds completed!"))
    }
  }
  print(paste0("Validation for ",k*10,"% of samples used for validation for ",gs_method," method completed!"))
}

# Remember to change the parameters in the earlier lines (line 39) for each GWAS model!!!!!!
# And also the file reading in line 20 to 22!
new_corr_df_mlm <- pivot_longer(corr_df, `10%`:`90%`, names_to = "percentage_validation", 
                                values_to = "correlation_coefficient") ; new_corr_df_mlm$model = "MLM"
new_corr_df_glm <- pivot_longer(corr_df, `10%`:`90%`, names_to = "percentage_validation", 
                                values_to = "correlation_coefficient") ; new_corr_df_glm$model = "GLM"
new_corr_df_super <- pivot_longer(corr_df, `10%`:`90%`, names_to = "percentage_validation", 
                                  values_to = "correlation_coefficient") ; new_corr_df_super$model = "SUPER"
new_corr_df <- rbind(new_corr_df_mlm, new_corr_df_glm, new_corr_df_super) %>% na.omit()
new_corr_df$model <- factor(new_corr_df$model, levels = c("MLM","GLM","SUPER"))

ggplot(new_corr_df, aes(x = percentage_validation, y = correlation_coefficient)) +
  geom_boxplot(aes(fill = model), outlier.shape = 1, outlier.size = 3, size = 0.7) +
  labs(title = paste0("Prediction accuracy for GS validation (",gs_method,") - comparison between GWAS models"), 
       caption = paste0("For each percentage of samples used, ",no_rounds," rounds of validation were performed"),
       x = "Percentage of samples used for validation", 
       y = "Correlation coefficient",
       fill = "GWAS model") +
  theme_minimal() +
  theme(plot.title = element_text(color="black", size=14, face="bold"),
        axis.title.x = element_text(color="black", size=12, margin = margin(t = 16)),
        axis.title.y = element_text(color="black", size=12, margin = margin(r = 10)))
