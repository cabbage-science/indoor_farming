###########################
#### LOADING LIBRARIES ####
###########################

getwd()
library(caret) #5 fold cross validation for different models
library(tidyverse) 
library(readxl)
library(pls)

####################################
#### LOADING AND PREPARING DATA ####
####################################

# Loading metabolite information
metabolites <- read.csv("data/kale_compound_ids_mzrt.csv", header = TRUE)

# Loading antioxidant capacity phenotype data (ABTS)
df2 <- read_excel("data/kale_ABTS_data.xlsx", col_names = TRUE)
# Removing repeated samples
AC_values <- df2[!duplicated(df2$Sample_ID),]

# Loading LCMS peak data (unscaled, with some sample repeats)
df1 <- read.csv("data/kale_peaks_noQC.csv", header = TRUE) %>%
  # Removing redundant columns
  select(-(1:5)) %>%
  # Rename column
  rename(Sample_ID = OriginalID) %>%
  # Configuring the sample IDs
  mutate(Sample_ID = gsub("\\.",":",toupper(Sample_ID))) %>%
  # Filtering only for samples that have been assayed with ABTS
  semi_join(AC_values, by = "Sample_ID")

# Removing repeated samples
kale_peaks_unscaled <- df1[!duplicated(df1$Sample_ID),]

# Creating temporary data frames
kale_ids <- kale_peaks_unscaled %>% select(1)
kale_peaks_unscaled_peaks <- kale_peaks_unscaled %>% select(-1)

# Creating pareto scaling function
paretoscale <- function(data) {
  x <- data
  x.centered <- apply(x, 2, function(x) x - mean(x))
  x.sc <- data.frame(apply(x.centered, 2, function(x) x/sqrt(sd(x))))
  return(x.sc)
}

# Applying pareto scaling function on unscaled peak intensities
kale_peaks_scaled_peaks <- paretoscale(kale_peaks_unscaled_peaks)

# Re-creating the original dataframe, but scaled
kale_peaks_scaled <- cbind(kale_ids, kale_peaks_scaled_peaks)

# Clearing environment
rm(df1, df2, kale_peaks_unscaled_peaks, kale_peaks_scaled_peaks, kale_peaks_unscaled)

#######################################
#### PLS MODELLING AND VIP SCORING ####
#######################################

# Merging the phenotype df and LCMS df, to include only samples that undergone LCMS
merged_df <- right_join(AC_values, kale_peaks_scaled, by = "Sample_ID") %>%
  # Removing Sample ID column
  select(-1)

# Removing sample IDs for kale_peaks_scaled df
kale_peaks_scaled_noid <- select(kale_peaks_scaled, -1)

# Running PLS model and calculating VIP scores for each metabolite
library(mixOmics) #select function from tidyverse will not work after this library is activated
PLS_model <- pls(kale_peaks_scaled_noid, merged_df$AC)
VIP_values <- vip(PLS_model)

# Remember to check off the mixomics package here!!!
# Filtering for metabolites with VIP scores more than 2
VIP_values_large <- data.frame(VIP_values) %>% 
  filter(comp1 >= 2)

# Creating a data frame of metabolite_IDs with VIP more than 2
metabolite_info <- data.frame(rownames(VIP_values_large)) %>%
  # Renaming rows for easier reference downstream
  rename(metabolite_ID = rownames.VIP_values_large.) %>%
  # Preparing metabolite IDs for merging left join below
  mutate(metabolite_ID = as.integer(str_replace(metabolite_ID, "align_", ""))) %>%
  # Extracting information from metabolite information data frame
  left_join(metabolites, by = c("metabolite_ID" = "Alignment.ID"))

