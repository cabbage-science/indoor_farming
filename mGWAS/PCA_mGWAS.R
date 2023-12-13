##########################
#### LOADING PACKAGES ####
##########################

getwd()
library(tidyverse)
library(data.table)
library(readxl)
library(ggplot2)

###############################
#### EIGENVALUE PROCESSING ####
###############################

kale_eigenvalues <- fread("kale_371_LCMS_eigenvalues.txt", header = TRUE) %>%
  # Multiplying cumulative proportion by 100 to reflect % in plot
  mutate(`cumulative proportion` = `cumulative proportion`*100) %>%
  # For plotting of correct positions for PCs (x-axis)
  mutate(PC_no = PC + 1) %>%
  # Renaming column names for easier ggploting
  rename("proportion" = "proportion of total",
         "cumulative_proportion" = "cumulative proportion") %>%
  # Multiplying proportion by 100 to reflect % in plot
  mutate(proportion = proportion*100)

# Eigenvalue plot
ggplot(data = kale_eigenvalues, aes(x = PC_no, y = proportion)) +
  geom_point(shape=19, color="black", size=3) + 
  xlim(0,30) + 
  theme_minimal() + 
  labs(title = "Allelic variation controlled by top 30 PCs (371 samples)",
       x = "Principal Component", y = "Allelic variation (%)") +
  theme(plot.title = element_text(color="black", size=14, face="bold", hjust = 0.5),
        axis.title.x = element_text(color="black", size=12, face="bold", margin = margin(t = 10)),
        axis.title.y = element_text(color="black", size=12, face="bold", margin = margin(r = 10)),
        axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=10)) +
  scale_y_continuous(breaks = seq(0, 12, by = 2))

# Cumulative proportion plot
ggplot(data = kale_eigenvalues, aes(x = PC_no, y = cumulative_proportion)) +
  geom_point(shape=19, color="black", size=3) + 
  xlim(0,30) + ylim(0,40) +
  theme_minimal() + 
  labs(title = "Cumulative proportion controlled by top 30 PCs (371 samples)",
       x = "Principal Component", y = "Cumulative proportion (%)") +
  theme(plot.title = element_text(color="black", size=14, face="bold", hjust = 0.5),
        axis.title.x = element_text(color="black", size=12, face="bold", margin = margin(t = 10)),
        axis.title.y = element_text(color="black", size=12, face="bold", margin = margin(r = 10)),
        axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=10))

################################
#### EIGENVECTOR PROCESSING ####
################################

# Loading sample IDs in sequence
sample_ids <- read_xlsx("pca_sample_ids.xlsx", col_names = TRUE)

# Loading eigenvectors txt file
kale_eigenvectors <- fread("kale_371_LCMS_eigenvectors.txt", header = TRUE)
# Changing the sample IDs in the eigenvectors df
kale_eigenvectors[,1] <- sample_ids$Sample_ID

# Filtering for samples clustering at bottom right - lragest group of samples to be group 1
kale_371_grp1 <- kale_eigenvectors %>%
  filter(PC1 > -25, PC2 < 20)

# Plotting PC1 and PC2
ggplot(kale_eigenvectors, aes(x = PC1, y = PC2)) +
  geom_point(shape = 21, size = 3, color = "black") +
  labs(title = "PCA Plot (PC1,PC2)", 
       x = "PC1 (10.73%)", 
       y = "PC2 (4.33%)") +  
  theme_minimal() +
  xlim(-150, 40) +  ylim(-60, 60) +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 15),              
    axis.text = element_text(size = 12),
    axis.title.y = element_text(margin = margin(r = 15)),
    axis.title.x = element_text(margin = margin(t = 15)))

# Isolating the sample IDs of group 1 samples
grp1_sample_ids <- kale_371_grp1[,1] %>% rename(Sample_ID = Taxa)
# Exporting IDs of samples in group 1 into a csv file
write.csv(grp1_sample_ids, "grp1_sample_ids.csv", quote = FALSE, row.names = FALSE)

