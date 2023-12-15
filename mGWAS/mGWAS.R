###########################################
#### LOADING PACKAGES AND SOURCE CODES ####
###########################################

getwd()
library(data.table)
library(readxl)
source("http://zzlab.net/GAPIT/GAPIT.library.R")
source("http://zzlab.net/GAPIT/gapit_functions.txt")
library(tidyverse)

##################################################
#### CONVERTING HAPMAP DIPLOID INTO NUMERICAL ####
##################################################

# Loading in test data (only applicable for conversion into numerical)
myG_ <- data.frame(fread("../data/mdp_genotype_test.hmp.txt", head=FALSE))

# Loading in hapmap diploid data
myG <- data.frame(fread("../data/kale_371_sequenced_LCMS.hmp.txt", head = FALSE))
# Ensuring consistency with test data
temp_names <- as.vector(as.matrix(myG_[1,1:11]))
myG[1,1:11] <- temp_names
# Conversion, and output a numerical txt file
myGAPIT <- GAPIT(G=myG, output.numerical=TRUE)

# Miscellaneous codes: Creating a csv file of sample IDs, which would be further
# Manipulated in excel to obtain sample names as K4:100:1:1 format (example) (DONE)
# Resulting excel file would be loaded in as "sample_ids" under "LOADING FILES"
ids <- as.vector(as.matrix(myG[1,12:382]))
id_df <- data.frame(ids)
write.csv(id_df, "genotyped_LCMS_sample_IDs.csv", row.names = FALSE, quote = FALSE)

#######################
#### LOADING FILES ####
#######################

# Loading in sample ID csv file
sample_ids <- read_xlsx("../data/sample_ids.xlsx", col_names = TRUE)

# Loading genotype map, modifying variables for correct data structure (following the GAPIT sample data)
myGM <- read.table("../data/kale_genotype_map.txt",head=TRUE) %>%
  mutate(Chromosome = as.integer(str_sub(Chromosome,3,3)))
# Creating a vector for myGD renaming
SNP <- myGM$SNP

# Loading numerical data (myGD)
myGD <- data.frame(fread("../data/kale_LCMS_numeric.txt", head=TRUE))
# Renaming sample names to simpler version
myGD$taxa <- sample_ids$Sample_ID
# Renaming column names to match genotype map 
colnames(myGD) <- c("taxa",SNP)

# Loading LCMS peak intensity data (phenotype - myY)
#sig_metabolite_intensities_genotyped <- read.csv("../data/significant_metabolite_peaks_genotyped.csv", header = TRUE)
sig_metabolite_intensities_genotyped <- read.csv("../data/significant_metabolite_peaks_genotyped_log2.csv", header = TRUE)

# Loading PCA data
myCV <- data.frame(fread("../PCA/kale_371_LCMS_eigenvectors.txt", header = TRUE))

###############################
#### mGWAS for ALL samples ####
###############################

# Choosing specific metabolites to include in myY
myY <- sig_metabolite_intensities_genotyped[,c(1,2)]

myGAPIT <- GAPIT(
  Y=myY,
  GD=myGD,
  GM=myGM,
  CV=myCV,
  model="GLM"
)

#############################################
#### mGWAS for group 1 samples (n = 311) ####
#############################################

# IMPORTANT: Please run all codes "LOADING FILES" first
grp1_sample_ids <- read.csv("../data/grp1_sample_ids.csv", header = TRUE)

# Creating myY. Remember to load in metabolite data from above
myY <- sig_metabolite_intensities_genotyped %>%
  # Filtering for grp1 samples only
  semi_join(grp1_sample_ids, by = "Sample_ID") %>%
  # Select specific metabolites 
  select(1,2)

# Creating a myGD for group 1. myGM would be the same.
myGD <- myGD %>%
  # Filtering for grp1 samples only
  semi_join(grp1_sample_ids, by = c("taxa" = "Sample_ID"))

myGAPIT <- GAPIT(
  Y=myY,
  GD=myGD,
  GM=myGM,
  CV=myCV,
  model="GLM"
)
