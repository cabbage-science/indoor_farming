###########################
#### LOADING LIBRARIES ####
###########################

library(caret) #5 fold cross validation for different models
library(tidyverse) 
library(readxl)
library(ggplot2)
library(reshape2)

####################################
#### LOADING AND PREPARING DATA ####
####################################

# Loading metabolite peak intensity data
sig_metabolites_kale_ID <- read.csv("data/significant_metabolite_peaks_log2.csv", header = TRUE)

# Removing kale ID names for calculation of correlation coefficients 
sig_metabolite_peaks <- sig_metabolites_kale_ID %>%
  select(-1)

##########################################
#### OPTIMIZATION - MULTICOLLINEARITY ####
##########################################

# Calculation of coefficients to produce a correlation matrix
correlation_matrix <- as.data.frame(cor(sig_metabolite_peaks)) %>%
  # Providing IDs of metabolites for ggplotting
  mutate(metabolite_ID2 = row.names(correlation_matrix))

# Converting into tidy data format 
correlation_long_symmetric <- pivot_longer(correlation_matrix, 1:68 , names_to = "metabolite_ID1", values_to = "correlation_coeff")

# Need to remove the correlation coefficients of repeated comparisons so that the resulting
# Heatmap will not be symmetric

# Plotting
ggplot(correlation_long, aes(x = metabolite_ID1, y = metabolite_ID2, fill = correlation_coeff)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

#####################################################
#### LASSO - 5 FOLD CROSS VALIDATION USING CARET ####
#####################################################


###################################################
#### LASSO - STEPWISE INCREMENT FOR PREDICTION ####
###################################################

