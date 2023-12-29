###########################
#### LOADING LIBRARIES ####
###########################

getwd()
library(mixOmics) #select function from tidyverse will not work after this library is activated
library(tidyverse) 
library(readxl)
library(ggplot2)
#library(pls)

####################################
#### LOADING AND PREPARING DATA ####
####################################

# Loading metabolite information
metabolites <- read.csv("data/arugula_compound_ids_mzrt.csv", header = TRUE)

# Loading antioxidant capacity phenotype data (ABTS)
df2 <- read_excel("data/arugula_ABTS_data.xlsx", col_names = TRUE)
# Removing repeated sample IDs
AC_values <- df2[!duplicated(df2$Sample_ID),]

# Loading LCMS peak data (unscaled, with some sample repeats)
df1 <- read.csv("data/arugula_peaks_noQC.csv", header = TRUE) %>%
  # Removing redundant columns
  select(-6) %>% select(-(1:4)) %>%
  # Rename column
  rename(Sample_ID = OriginalID) %>%
  # Configuring the sample IDs
  mutate(Sample_ID = gsub("\\.",":",toupper(Sample_ID))) %>%
  # Filtering only for samples that have been assayed with ABTS
  semi_join(AC_values, by = "Sample_ID")

# Removing repeated samples
arugula_peaks_unscaled <- df1[!duplicated(df1$Sample_ID),]

# Creating temporary data frames
arugula_ids <- arugula_peaks_unscaled %>% select(1)

# Clearing environment
rm(df1, df2)

#######################################
#### PLS MODELLING AND VIP SCORING ####
#######################################

# Log transforming the metabolite peak intensities
arugula_peaks_scaled <- mutate_if(arugula_peaks_unscaled, is.numeric, function(x) log2(x)) ; rm(arugula_peaks_unscaled)

# Merging the phenotype df and LCMS df, to include only samples that undergone LCMS
AC_merged <- left_join(arugula_peaks_scaled, AC_values, by = "Sample_ID") %>%
  # Selecting only the phenotype column(s)
  select("AC")
#saveRDS(AC_merged, file = "AC_merged.rds")

# Removing sample IDs for kale_peaks_scaled df
arugula_peaks_scaled_noid <- select(arugula_peaks_scaled, -1)

# Running PLS model and calculating VIP scores for each metabolite
PLS_model <- pls(arugula_peaks_scaled_noid, AC_merged)
VIP_values <- vip(PLS_model)

# Filtering for metabolites with VIP scores more than 2
VIP_values_large <- data.frame(VIP_values) %>% 
  filter(comp1 >= 2)
#write.csv(VIP_values_large, "metabolite_VIP_score.csv", quote = FALSE) #Then add name for first column in csv file

# Filtering for metabolites with VIP scores more than 1.5
#VIP_values_large <- data.frame(VIP_values) %>% 
#  filter(comp1 >= 1.5)

# Creating a data frame of metabolite_IDs with VIP more than threshold specified
metabolite_info <- data.frame(rownames(VIP_values_large)) %>%
  # Renaming rows for easier reference downstream
  rename(metabolite_ID = rownames.VIP_values_large.) %>%
  # Extracting information from metabolite information data frame
  left_join(metabolites, by = "metabolite_ID")

rm(PLS_model, VIP_values, AC_values)

#######################################################
#### LINEAR REGRESSION OF PEAK DATA WITH PHENOTYPE ####
#######################################################

# Selecting peak intensity data of significant metabolites, by extracting names of metabolites
sig_metabolite_data <- arugula_peaks_scaled_noid[,rownames(VIP_values_large)]

# Initializing empty vector
phenotype_corr_values <- NULL

# Calculating the phenotypic correlation with peak intensity for each of the significant metabolites
for(k in 1:ncol(sig_metabolite_data)) {
  phenotype_corr_values <- c(phenotype_corr_values, cor(AC_merged, sig_metabolite_data[,k]))
}

# Creating a data frame to store + visualize correlation coefficients
metabolite_corr <- data.frame(rownames(VIP_values_large), phenotype_corr_values)

# Extract only metabolites with correlation > 0, with antioxidant capacity 
#metabolite_corr <- filter(metabolite_corr, phenotype_corr_values >= 0) ; sig_metabolite_data <- kale_peaks_scaled_noid[,metabolite_corr$rownames.VIP_values_large.]

# Combining kale IDs andsignificant metabolite peak intensities into single df for export - downstream GWAS phenotype file
sig_metabolites <- cbind(arugula_ids, sig_metabolite_data)

# Writing csv file
write.csv(sig_metabolites, file = "significant_metabolite_peaks_log2_arugula.csv", quote = FALSE, row.names = FALSE)
write.csv(metabolite_corr, file = "significant_metabolite_correlation_arugula.csv", quote = FALSE, row.names = FALSE)

# Writing metabolite information file
metabolite_info_export <- metabolite_info %>%
  mutate(metabolite_ID = str_c("align_", metabolite_ID))
write.csv(metabolite_info_export, file = "significant_metabolite_info_arugula.csv", quote = FALSE, row.names = FALSE)


##########################################################
#### GGPLOT2 FOR PHENOTYPIC CORRELATION VISUALIZATION ####
##########################################################

# Combining phenotype and significant metabolite peak intensities into single df for ggplot
combined_df <- cbind(AC_merged, sig_metabolite_data)


##########################################################
#### GGPLOT2 FOR PHENOTYPIC CORRELATION VISUALIZATION ####
##########################################################

# Removing objects to clear environment just for visualization
rm(kale_ids, arugula_peaks_scaled_noid, arugula_peaks_scaled, metabolites, VIP_values_large, sig_metabolites)

# Creating a metabolite names vector for convenience
metabolite_names <- colnames(sig_metabolite_data)

#setwd("C:/Users/ETHAN/OneDrive/Documents/MR_B_BY/NUS/INDOOR_FARMING/LC-MS/phenotypic_correlation_new")

for (k in 1:ncol(sig_metabolite_data)) {
  temp_met_name <- metabolite_names[k]
  combined_df <- cbind(AC_merged, sig_metabolite_data[,temp_met_name]) %>%
    rename("x" = "sig_metabolite_data[, temp_met_name]")
  p <- ggplot(data = combined_df, aes(x = x, y = AC)) +
    theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
          axis.title.y = element_text(margin = margin(r = 20)),
          axis.title.x = element_text(margin = margin(t = 20)),
          axis.title = element_text(size = 14),              
          axis.text = element_text(size = 10),
          panel.grid.major = element_line(color = "lightgray", linewidth = 0.3),
          panel.grid.minor = element_line(color = "lightgray", linewidth = 0.3),
          panel.background = element_rect(fill = "white")) + 
    geom_smooth(method = "lm", se = FALSE, color = "red",alpha = 0.7, lty = "dotted") +
    geom_point(shape = 21, size = 3, color = "black") +
    labs(title = paste0("Relationship between metabolite ID ", temp_met_name, " and antioxidant capacity"), 
         x = "Peak intensity (log2 transformed)", 
         y = "ABTS Antioxidant capacity (µmol TE/g)") +
    ylim(0,180) +
    annotate("text", x = -Inf, y = -Inf, 
             # Extracting correlation coefficient
             label = paste0("R = ",round(metabolite_corr$phenotype_corr_values[k],3)), 
             hjust = -0.4, vjust = -2, color = "black", size = 8)
  ggsave(paste0(temp_met_name,"_phenotypic_correlation_arugula.png"), plot = p, width = 10, height = 8, dpi = 800)
}
