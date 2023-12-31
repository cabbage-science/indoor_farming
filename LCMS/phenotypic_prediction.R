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

# Loading metabolite peak intensity data (after removing highly correlated metabolites)
sig_metabolite_peaks_removed <- read.csv("data/significant_metabolite_peaks_log2_removed.csv", header = TRUE) %>%
  select(-1)

# Loading antioxidant capacity data
AC_merged <- readRDS("data/AC_merged.rds")

##########################################
#### OPTIMIZATION - MULTICOLLINEARITY ####
##########################################

# Calculation of coefficients to produce a correlation matrix
correlation_matrix <- as.data.frame(cor(sig_metabolite_peaks)) 
# Providing IDs of metabolites for ggplotting
correlation_matrix <- mutate(correlation_matrix, metabolite_ID2 = row.names(correlation_matrix))

# Converting into tidy data format 
correlation_long_symmetric <- pivot_longer(correlation_matrix, 1:68 , names_to = "metabolite_ID1", values_to = "correlation_coeff")

# Need to remove the correlation coefficients of repeated comparisons so that the resulting
# Heatmap will not be symmetric

# Plotting
ggplot(correlation_long_symmetric, aes(x = metabolite_ID1, y = metabolite_ID2, fill = correlation_coeff)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

# Removing correlated metabolites - checking. Removal done manually. Refer to section below

#####################################################
#### LASSO - 5 FOLD CROSS VALIDATION USING CARET ####
#####################################################

# Set parameters for model training. Number = n-fold cross validation
repeats = 1 ; number = 5
# For randomness. Can add a "#" in front for reproducibility
seeds <- sample(1:10000000, size = (repeats*number) + 1, replace = FALSE)

# Defining parameters for model training 
train_control <- trainControl(method = "repeatedcv", repeats=repeats,
                              number = number,
                              seeds = seeds)

# Creating a merged dataframe for AC values and all metabolite peaks for all samples
merged_df <- cbind(AC_merged,sig_metabolite_peaks)

# Creating the predictive model
modelLASSO <- train(AC ~., data = merged_df,
                          method = "lasso", 
                          trControl = train_control)

# Creating a list for caret input for extractPrediction function
bothModels <- list(lasso = modelLASSO)

# Creating a df of actual and predicted values from the model above
actual_predicted <- data.frame(extractPrediction(bothModels, testX = sig_metabolite_peaks))

ggplot(actual_predicted, aes(x=pred, y=obs)) + 
  geom_point(shape = 21, colour = "black", size = 2, stroke = 1) +
  geom_abline(intercept=0, slope=1, colour = "red") +
  labs(x='Predicted Brix Values', y='Actual Brix Values') +
  theme_minimal() +
  theme(text=element_text(size=20))

# Calculate average R squared value from the model (k repeats and n folds)
mean(modelLASSO$resample$Rsquared)

######################################################
#### LASSO - STEPWISE INCREMENT FOR AC PREDICTION ####
######################################################

# Naming phenotype
phenotype <- "Antioxidant Capacity"
# Naming model
predict_model <- "LASSO"

# Loading metabolite IDs and sorting metabolites in decreasing VIP score
VIP_values <- read.csv("data/metabolite_VIP_score.csv") %>%
  # Removing comp2 column
  select(1,2) %>%
  arrange(desc(comp1))
# Creating a vector for selection of metabolites in for loop later on
metabolite_vector <- VIP_values$metabolite_ID

# Set parameters for model training. Number = n-fold cross validation
repeats = 100 ; number = 5

# Creating a merged dataframe for AC values and all metabolite peaks for all samples
merged_df <- cbind(AC_merged,sig_metabolite_peaks)

# Initializing an empty vector for for-loop
R2_values <- NULL


# For loop to increase no. of metabolites stepwise
for (k in 2:ncol(sig_metabolite_peaks)) {
  # For randomness. Can add a "#" in front for reproducibility
  seeds <- sample(1:10000000, size = (repeats*number) + 1, replace = FALSE)
  
  # Defining parameters for model training 
  train_control <- trainControl(method = "repeatedcv", 
                                repeats = repeats,
                                number = number,
                                seeds = seeds)
  
  # Creating the data frame consisting of AC values and selected metabolites
  temp_df <- merged_df %>% select(1, metabolite_vector[1:k])
  
  # Creating the predictive model
  temp_modelLASSO <- train(AC ~., data = temp_df,
                      method = "lasso", 
                      trControl = train_control)
  
  # Catenating R2 values to the vector
  R2_values <- c(R2_values,mean(temp_modelLASSO$resample$Rsquared, na.rm = TRUE))
  print(paste0("Model building for ",k," metabolites complete"))
  
  # Generating a data frame of observed and predicted values using caret
  actual_predicted <- data.frame(extractPrediction(list(lasso = temp_modelLASSO), testX = temp_df[,-1]))
  
  # Plotting
  p <- ggplot(actual_predicted, aes(x=pred, y=obs)) + 
    geom_abline(intercept=0, slope=1, colour = "red") +
    geom_point(shape = 21, colour = "black", size = 2, stroke = 1) +
    labs(x='Predicted Antioxidant capacity Values', y='Observed Antioxidant capacity Values',
         title = paste0("Prediction for ",phenotype," using ",k," metabolites (",predict_model," model)"),
         subtitle = "Metabolites used for prediction were arranged in decreasing VIP score",
         caption = "Red line indicates the expected linear relationship between observed and predicted values") +
    theme(plot.title = element_text(size = 18, face = "bold"),
          plot.caption = element_text(size = 10),
          plot.subtitle = element_text(size = 12),
          axis.title.y = element_text(margin = margin(r = 20)),
          axis.title.x = element_text(margin = margin(t = 20)),
          axis.title = element_text(size = 14),              
          axis.text = element_text(size = 10),
          panel.grid.major = element_line(color = "lightgray", linewidth = 0.3),
          panel.grid.minor = element_line(color = "lightgray", linewidth = 0.3),
          panel.background = element_rect(fill = "white")) +
    theme(text=element_text(size=20)) +
    annotate("text", x = Inf, y = -Inf, 
             # Extracting correlation coefficient
             label = bquote("R"^2 == .(round(R2_values[k-1],3))), 
             hjust = 1.4, vjust = -1.5, color = "black", size = 7)
  ggsave(paste0(k,"_metabolites_",phenotype,"_prediction.png"), plot = p, width = 10, height = 7, dpi = 600)
} ; print(R2_values)

# Preparing csv file for export
No_of_metabolites <- seq(2,1+length(R2_values),by = 1)
R2_df <- data.frame(No_of_metabolites, R2_values)
write.csv(R2_df, "R2_metabolite_AC_prediction.csv", quote = FALSE, row.names = FALSE)

# Base R plotting for improvements in R2 with stepwise increments of metabolites
plot(No_of_metabolites, R2_values)

#####################################################
#### LASSO - STEPWISE INCREMENT FOR AC PREDICTION ###
#### AFTER ACCOUNTING FOR CORRELATED METABOLITES ####
#### AND REMOVING THEM ##############################
#####################################################

# Manual bootstrapping - removal of metabolites
correlation_over_0.8_1 <- correlation_long_symmetric %>%
  filter(correlation_coeff > 0.8, correlation_coeff < 1)

remove1 <- paste0("align_",c(3687, 5783, 4282, 4437, 4439, 6516, 6517, 6518))
correlation_long_1 <- correlation_long_symmetric %>%
  filter(!metabolite_ID2 %in% remove1) %>%
  filter(!metabolite_ID1 %in% remove1)

correlation_over_0.8_2 <- correlation_long_1 %>%
  filter(correlation_coeff > 0.8, correlation_coeff < 1)

remove2 <- paste0("align_",c(2849, 2854, 3155, 3245, 3246, 4344, 4346, 5013, 5014, 
                             5015, 5018, 5165, 5168, 5173, 6593, 7386, 7387))
correlation_long_2 <- correlation_long_1 %>%
  filter(!metabolite_ID2 %in% remove2) %>%
  filter(!metabolite_ID1 %in% remove2)

correlation_over_0.8_3 <- correlation_long_2 %>%
  filter(correlation_coeff > 0.8, correlation_coeff < 1)

remove3 <- paste0("align_", c(5171, 5172, 5335, 6070, 8994, 9339, 11895, 12291, 14368, 19245))
correlation_long_3 <- correlation_long_2 %>%
  filter(!metabolite_ID2 %in% remove3) %>%
  filter(!metabolite_ID1 %in% remove3)

correlation_over_0.8_4 <- correlation_long_3 %>%
  filter(correlation_coeff > 0.8, correlation_coeff < 1) # no more metabolites significantly correlated

# correlation_long_3 is the final dataframe with 33 metabolites.
# Plotting to visualize and confirm that correlations < 0.8
ggplot(correlation_long_3, aes(x = metabolite_ID1, y = metabolite_ID2, fill = correlation_coeff)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

#### PREPARING NEW CSV FILE CONTAINING ONLY PEAK INTENSITY DATA OF METABOLITES AFTER REMOVAL
# Creating new df with removed metabolites
#sig_metabolite_peaks_new <- sig_metabolite_peaks %>%
  #select(unique(correlation_long_3$metabolite_ID2))
#sig_metabolite_peaks_new <- cbind(sig_metabolites_kale_ID[,1], sig_metabolite_peaks_new)
#write.csv(sig_metabolite_peaks_new, "significant_metabolite_peaks_log2_removed.csv", quote = FALSE, row.names = FALSE)

### -------------------------------------------------------------------------- ###
# ANALYSIS STARTS HERE

# Creating vector of metabolite IDs of those that remained after pruning/removal
remaining_metabolites <- colnames(sig_metabolite_peaks_removed)

# Naming phenotype
phenotype <- "Antioxidant Capacity"
# Naming model
predict_model <- "LASSO"

# Loading metabolite IDs and sorting metabolites in decreasing VIP score
VIP_values <- read.csv("data/metabolite_VIP_score.csv") %>%
  # Removing comp2 column
  select(1,2) %>%
  filter(metabolite_ID %in% remaining_metabolites) %>%
  arrange(desc(comp1))

# Set parameters for model training. Number = n-fold cross validation
repeats = 100 ; number = 5

# Creating a merged dataframe for AC values and all metabolite peaks for all samples
merged_df <- cbind(AC_merged,sig_metabolite_peaks_removed)

# Initializing an empty vector for for-loop
R2_values <- NULL


# For loop to increase no. of metabolites stepwise
for (k in 2:ncol(sig_metabolite_peaks_removed)) {
  # For randomness. Can add a "#" in front for reproducibility
  seeds <- sample(1:10000000, size = (repeats*number) + 1, replace = FALSE)
  
  # Defining parameters for model training 
  train_control <- trainControl(method = "repeatedcv", 
                                repeats = repeats,
                                number = number,
                                seeds = seeds)
  
  # Creating the data frame consisting of AC values and selected metabolites
  temp_df <- merged_df %>% select(1, remaining_metabolites[1:k])
  
  # Creating the predictive model
  temp_modelLASSO <- train(AC ~., data = temp_df,
                           method = "lasso", 
                           trControl = train_control)
  
  # Catenating R2 values to the vector
  R2_values <- c(R2_values,mean(temp_modelLASSO$resample$Rsquared, na.rm = TRUE))
  print(paste0("Model building for ",k," metabolites complete"))
  
  # Generating a data frame of observed and predicted values using caret
  actual_predicted <- data.frame(extractPrediction(list(lasso = temp_modelLASSO), testX = temp_df[,-1]))
  
  # Plotting
  p <- ggplot(actual_predicted, aes(x=pred, y=obs)) + 
    geom_abline(intercept=0, slope=1, colour = "red") +
    geom_point(shape = 21, colour = "black", size = 2, stroke = 1) +
    labs(x='Predicted Antioxidant capacity Values', y='Observed Antioxidant capacity Values',
         title = paste0("Prediction for ",phenotype," using ",k," metabolites (",predict_model," model)"),
         subtitle = "Metabolites used for prediction were arranged in decreasing VIP score",
         caption = "Red line indicates the expected linear relationship between observed and predicted values") +
    theme(plot.title = element_text(size = 18, face = "bold"),
          plot.caption = element_text(size = 10),
          plot.subtitle = element_text(size = 12),
          axis.title.y = element_text(margin = margin(r = 20)),
          axis.title.x = element_text(margin = margin(t = 20)),
          axis.title = element_text(size = 14),              
          axis.text = element_text(size = 10),
          panel.grid.major = element_line(color = "lightgray", linewidth = 0.3),
          panel.grid.minor = element_line(color = "lightgray", linewidth = 0.3),
          panel.background = element_rect(fill = "white")) +
    theme(text=element_text(size=20)) +
    annotate("text", x = Inf, y = -Inf, 
             # Extracting correlation coefficient
             label = bquote("R"^2 == .(round(R2_values[k-1],3))), 
             hjust = 1.4, vjust = -1.5, color = "black", size = 7)
  ggsave(paste0(k,"_metabolites_",phenotype,"_prediction.png"), plot = p, width = 10, height = 7, dpi = 600)
} ; print(R2_values)

# Preparing csv file for export
No_of_metabolites <- seq(2,1+length(R2_values),by = 1)
R2_df <- data.frame(No_of_metabolites, R2_values)
write.csv(R2_df, "R2_metabolite_AC_prediction_removed.csv", quote = FALSE, row.names = FALSE)

# Base R plotting for improvements in R2 with stepwise increments of metabolites
plot(No_of_metabolites, R2_values)

#########################################################
#### LASSO - CORRELATED METABOLITES REMOVED AND #########
#### TOGETHER WITH POSITIVELY CORRELATED METABOLITES ####
#########################################################

# Only run this if you want to keep optimal number of metabolites
sig_metabolite_peaks_removed <- sig_metabolite_peaks_removed %>%
  select(1:7)

# Two of the positively correlated metabolites are highly correlated (27373 and 27380)
# Can consider removing one of them during prediction

positively_correlated_metabolites <- read.csv("data/positively_correlated_metabolite_data.csv", header = TRUE) %>%
  select(-align_27380) %>%
  select(-Sample_ID)

merged_df_new <- cbind(AC_merged, sig_metabolite_peaks_removed, positively_correlated_metabolites)

# Set parameters for model training. Number = n-fold cross validation
repeats = 100 ; number = 5
# For randomness. Can add a "#" in front for reproducibility
seeds <- sample(1:10000000, size = (repeats*number) + 1, replace = FALSE)

# Defining parameters for model training 
train_control <- trainControl(method = "repeatedcv", repeats=repeats,
                              number = number,
                              seeds = seeds)

# Creating the predictive model
modelLASSO <- train(AC ~., data = merged_df_new,
                    method = "lasso", 
                    trControl = train_control)

# Creating a list for caret input for extractPrediction function
bothModels <- list(lasso = modelLASSO)

# Creating a df of actual and predicted values from the model above
actual_predicted <- data.frame(extractPrediction(bothModels, testX = merged_df_new[,-1]))

p <- ggplot(actual_predicted, aes(x=pred, y=obs)) + 
  geom_abline(intercept=0, slope=1, colour = "red") +
  geom_point(shape = 21, colour = "black", size = 2, stroke = 1) +
  labs(x='Predicted Antioxidant capacity Values', y='Observed Antioxidant capacity Values',
       title = paste0("Prediction for ",phenotype," (",predict_model," model)"),
       subtitle = "Optimal number of negatively correlated (n = 6) and positively correlated (n = 4) metabolites were used",
       caption = "Red line indicates the expected linear relationship between observed and predicted values") +
  theme(plot.title = element_text(size = 18, face = "bold"),
        plot.caption = element_text(size = 10),
        plot.subtitle = element_text(size = 11),
        axis.title.y = element_text(margin = margin(r = 20)),
        axis.title.x = element_text(margin = margin(t = 20)),
        axis.title = element_text(size = 14),              
        axis.text = element_text(size = 10),
        panel.grid.major = element_line(color = "lightgray", linewidth = 0.3),
        panel.grid.minor = element_line(color = "lightgray", linewidth = 0.3),
        panel.background = element_rect(fill = "white")) +
  theme(text=element_text(size=20)) +
  annotate("text", x = Inf, y = -Inf, 
           # Extracting correlation coefficient
           label = bquote("R"^2 == .(round(mean(modelLASSO$resample$Rsquared),3))), 
           hjust = 1.4, vjust = -1.5, color = "black", size = 7)
ggsave("positive_and_optimal_metabolites_AC.png", plot = p, width = 10, height = 7, dpi = 600)

# Calculate average R squared value from the model (k repeats and n folds)
mean(modelLASSO$resample$Rsquared)


###############################################
#### COMPARING DIFFERENT PREDICTION MODELS ####
###############################################

# Set parameters for model training. Number = n-fold cross validation
repeats = 20 ; number = 5
# For randomness. Can add a "#" in front for reproducibility
seeds <- sample(1:10000000, size = (repeats*number) + 1, replace = FALSE)

# Defining parameters for model training 
train_control <- trainControl(method = "repeatedcv", repeats=repeats,
                              number = number,
                              seeds = seeds)

# Creating a merged dataframe for AC values and all metabolite peaks for all samples
merged_df <- cbind(AC_merged,sig_metabolite_peaks)

# Creating the LASSO model
model_LASSO <- train(AC ~., data = merged_df,
                    method = "lasso", 
                    trControl = train_control)

# Creating the PLS model
model_PLS <- train(AC ~., data = merged_df,
                    method = "pls", 
                    trControl = train_control)

# Creating the Bayesian Ridge Regression model (bridge) / Bayesian GLM model (bayesglm)
model_bayesian <- train(AC ~., data = merged_df,
                    method = "bayesglm", 
                    trControl = train_control)

# Creating the backward model
model_leapBackward <- train(AC ~., data = merged_df,
                    method = "leapBackward", 
                    trControl = train_control)

# Creating a list for caret input for extractPrediction function
results <- resamples(list(PLS=model_PLS, LASSO=model_LASSO, Backward = model_leapBackward, Bayesian = model_bayesian))
# To view R2 boxplots
bwplot(results, xlim = c(0,1))
bwplot(results)


# Calculate average R squared value from the model (k repeats and n folds)
mean(model_PLS$resample$Rsquared)
mean(model_LASSO$resample$Rsquared)
mean(model_bayesian$resample$Rsquared)
mean(model_leapBackward$resample$Rsquared)
