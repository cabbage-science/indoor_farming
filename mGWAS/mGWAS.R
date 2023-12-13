###########################################
#### LOADING PACKAGES AND SOURCE CODES ####
###########################################

library(tidyverse)
library(data.table)
library(readxl)
source("http://zzlab.net/GAPIT/GAPIT.library.R")
source("http://zzlab.net/GAPIT/gapit_functions.txt")

##################################################
#### CONVERTING HAPMAP DIPLOID INTO NUMERICAL ####
##################################################

# Loading in test data (only applicable for conversion into numerical)
myG_ <- data.frame(fread("../data/mdp_genotype_test.hmp.txt", head=FALSE))

# Loading in hapmap diploid data
myG <- data.frame(fread("../data/kale_371_sequenced_LCMS.hmp.txt", head = FALSE))
# Ensuring consistency with test data
lmao <- as.vector(as.matrix(myG_[1,1:11]))
myG[1,1:11] <- lmao
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
sig_metabolite_intensities_genotyped <- read.csv("../data/significant_metabolite_peaks_genotyped.csv", header = TRUE)
# Choosing specific metabolites to include in myY
myY <- sig_metabolite_intensities_genotyped %>%
  select(1,2)

###############
#### mGWAS ####
###############
myGAPIT <- GAPIT(
  Y=myY,
  GD=myGD,
  GM=myGM,
  PCA.total=5,
  model="GLM"
)
