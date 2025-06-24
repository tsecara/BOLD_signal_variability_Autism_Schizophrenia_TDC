# Load necessary packages
#install.packages("car")
library(car)
library(dplyr)
library(broom)
#install.packages("pbapply")
library(pbapply) # For progress bar

# Going to do regional comparisons across groups to see if there are any differences 
###First for EA:
mssd_val <- read.csv("/projects/tsecara/SPINS_ASD_Project2/data/combined/regional_MSSD/neurocombat/combined_mssd_mean_EA_RScleanFINAL_agereduced.csv")
data <- read.csv("/projects/tsecara/SPINS_ASD_Project2/data/combined/demographics/demo_EA_agereduced.csv")
region_merged <- merge(data, mssd_val, by = 'record_id')

# Define function to perform ANOVA and TukeyHSD
perform_anova <- function(df, response_col, predictor_cols, tukey_if_significant = TRUE) {
  # Build the formula for the ANOVA model
  formula <- as.formula(paste(response_col, '~', paste(predictor_cols, collapse = ' + ')))
  
  # Run the ANOVA model
  model <- aov(formula, data = df)
  
  # Get the ANOVA table using car::Anova
  anova_results <- car::Anova(model, type = 'II')
  
  # Check if 'group' is significant
  group_significant <- anova_results$'Pr(>F)'[rownames(anova_results) == 'group'] < 0.05
  
  # If 'group' is significant, run TukeyHSD and store results
  if (group_significant & tukey_if_significant) {
    tukey_results <- TukeyHSD(aov(as.formula(paste(response_col, '~ group')), data = df))
    tukey_results <- tidy(tukey_results) %>%
      mutate(response_col = response_col)
  } else {
    tukey_results <- NULL
  }
  
  # Return both ANOVA and TukeyHSD results
  list(anova_results = tidy(anova_results) %>% mutate(response_col = response_col),
       tukey_results = tukey_results)
}

# Define the predictor columns
predictor_cols <- c('group', 'age', 'sex', 'avg_fd')

# Use pbapply to include a progress bar and process networks one by one
results_list <- pblapply(names(mssd_val)[-1], function(region_col) {
  perform_anova(region_merged, region_col, predictor_cols)
})

# Combine results into a final dataframe
anova_df <- bind_rows(lapply(results_list, `[[`, 'anova_results'))
tukey_df <- bind_rows(lapply(results_list, `[[`, 'tukey_results'))

# Apply FDR correction using the Benjamini-Hochberg method
tukey_df$fdr_adj_p.value <- p.adjust(tukey_df$adj.p.value, method = "fdr")

#Filter out the values:
filtered_tukey_df <- tukey_df[tukey_df$fdr_adj_p.value < 0.0250, ]
filtered_tukey_df_ASD <- filtered_tukey_df[filtered_tukey_df$contrast == "Control-ASD",]
filtered_tukey_df_SSD <- filtered_tukey_df[filtered_tukey_df$contrast == "SSD-Control",]
filtered_tukey_df_trans <- filtered_tukey_df[filtered_tukey_df$contrast == "SSD-ASD",]

####### Creating p-scaler files:
#####ASD First 
####### Postive values indicate the regions where Controls have higher MSSD than ASD, and Negative values indicuate where the ASD has higher MSSD values than controls 
# Extract column names excluding "record_id"
column_names <- colnames(mssd_val)[colnames(mssd_val) != "record_id"]

# Create a new dataframe with one column titled "p-val"
ASD_regional_pval <- data.frame(p.val = rep(NA, length(column_names)))
row.names(ASD_regional_pval) <- column_names

# Convert p.val column to character to accommodate signs
ASD_regional_pval$p.val <- as.character(ASD_regional_pval$p.val)

# Update ASD_regional_pval with values from filtered_tukey_df_ASD
for (i in seq_len(nrow(filtered_tukey_df_ASD))) {
  response <- filtered_tukey_df_ASD$response_col[i]
  fdr_value <- filtered_tukey_df_ASD$fdr_adj_p.value[i]
  estimate <- filtered_tukey_df_ASD$estimate[i]
  
  # Determine the sign based on the estimate
  sign <- ifelse(estimate > 0, "+", "-")
  
  # Match and update
  if (response %in% rownames(ASD_regional_pval)) {
    ASD_regional_pval[response, "p.val"] <- paste0(sign, sprintf("%.10f", fdr_value))
  }
}

# Set unmatched rows to "0"
ASD_regional_pval[is.na(ASD_regional_pval$p.val), "p.val"] <- "0"

# Print the updated structure to verify
str(ASD_regional_pval)
ASD_regional_pval$p.val <- as.numeric(ASD_regional_pval$p.val)
range(ASD_regional_pval$p.val)

# Convert the p.val column to numeric, ignoring signs for now
ASD_regional_pval$p.val_numeric <- as.numeric(gsub("[+-]", "", ASD_regional_pval$p.val))

# Apply log10 transformation, but only to non-zero values
ASD_regional_pval$p.val_log10 <- ifelse(ASD_regional_pval$p.val_numeric > 0, 
                                        -log10(ASD_regional_pval$p.val_numeric), 
                                        0)

# Restore the sign to the log-transformed values
ASD_regional_pval$p.val_log10 <- ifelse(grepl("-", ASD_regional_pval$p.val), 
                                        -ASD_regional_pval$p.val_log10, 
                                        ASD_regional_pval$p.val_log10)

# Keep only the p.val_log10 column
ASD_regional_pval <- ASD_regional_pval["p.val_log10"]
ASD_regional_pval <- as.data.frame(ASD_regional_pval[-c(1:32),])

write.table(ASD_regional_pval, file = "/projects/tsecara/SPINS_ASD_Project2/data/visualizations/P-scalers/ASD_TDC_regional_EA_adj.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

#####SSD VS Control 
####### Postive values indicate the regions where Controls have higher MSSD than SSD, and Negative values indicuate where the SSDs has higher MSSD values than controls 
# Extract column names excluding "record_id"
column_names <- colnames(mssd_val)[colnames(mssd_val) != "record_id"]

# Create a new dataframe with one column titled "p-val"
SSD_regional_pval <- data.frame(p.val = rep(NA, length(column_names)))
row.names(SSD_regional_pval) <- column_names

# Convert p.val column to character to accommodate signs
SSD_regional_pval$p.val <- as.character(SSD_regional_pval$p.val)

# Update SSD_regional_pval with values from filtered_tukey_df_SSD
for (i in seq_len(nrow(filtered_tukey_df_SSD))) {
  response <- filtered_tukey_df_SSD$response_col[i]
  fdr_value <- filtered_tukey_df_SSD$fdr_adj_p.value[i]
  estimate <- filtered_tukey_df_SSD$estimate[i]
  
  # Determine the sign based on the estimate
  sign <- ifelse(estimate > 0, "-", "+")
  
  # Match and update
  if (response %in% rownames(SSD_regional_pval)) {
    SSD_regional_pval[response, "p.val"] <- paste0(sign, sprintf("%.4f", fdr_value))
  }
}

# Set unmatched rows to "0"
SSD_regional_pval[is.na(SSD_regional_pval$p.val), "p.val"] <- "0"

#Changing back to numeric 
SSD_regional_pval$p.val <- as.numeric(SSD_regional_pval$p.val)
range(SSD_regional_pval$p.val)

# Convert the p.val column to numeric, ignoring signs for now
SSD_regional_pval$p.val_numeric <- as.numeric(gsub("[+-]", "", SSD_regional_pval$p.val))

# Apply log10 transformation, but only to non-zero values
SSD_regional_pval$p.val_log10 <- ifelse(SSD_regional_pval$p.val_numeric > 0, 
                                        -log10(SSD_regional_pval$p.val_numeric), 
                                        0)

# Restore the sign to the log-transformed values
SSD_regional_pval$p.val_log10 <- ifelse(grepl("-", SSD_regional_pval$p.val), 
                                        -SSD_regional_pval$p.val_log10, 
                                        SSD_regional_pval$p.val_log10)

# Keep only the p.val_log10 column
SSD_regional_pval <- SSD_regional_pval["p.val_log10"]
SSD_regional_pval <- as.data.frame(SSD_regional_pval[-c(1:32),])

write.table(SSD_regional_pval, file = "/projects/tsecara/SPINS_ASD_Project2/data/visualizations/P-scalers/SSD_TDC_regional_EA_adj.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

#####SSD vs ASD
####### Postive values indicate the regions where SSD have higher MSSD than ASD, and Negative values indicuate where the ASD has higher MSSD values than SSD 
# Extract column names excluding "record_id"
column_names <- colnames(mssd_val)[colnames(mssd_val) != "record_id"]

# Create a new dataframe with one column titled "p-val"
trans_regional_pval <- data.frame(p.val = rep(NA, length(column_names)))
row.names(trans_regional_pval) <- column_names

# Convert p.val column to character to accommodate signs
trans_regional_pval$p.val <- as.character(trans_regional_pval$p.val)

# Update trans_regional_pval with values from filtered_tukey_df_trans
for (i in seq_len(nrow(filtered_tukey_df_trans))) {
  response <- filtered_tukey_df_trans$response_col[i]
  fdr_value <- filtered_tukey_df_trans$fdr_adj_p.value[i]
  estimate <- filtered_tukey_df_trans$estimate[i]
  
  # Determine the sign based on the estimate
  sign <- ifelse(estimate > 0, "+", "-")
  
  # Match and update
  if (response %in% rownames(trans_regional_pval)) {
    trans_regional_pval[response, "p.val"] <- paste0(sign, sprintf("%.10f", fdr_value))
  }
}

# Set unmatched rows to "0"
trans_regional_pval[is.na(trans_regional_pval$p.val), "p.val"] <- "0"

trans_regional_pval$p.val <- as.numeric(trans_regional_pval$p.val)
range(trans_regional_pval$p.val)

# Convert the p.val column to numeric, ignoring signs for now
trans_regional_pval$p.val_numeric <- as.numeric(gsub("[+-]", "", trans_regional_pval$p.val))

# Apply log10 transformation, but only to non-zero values
trans_regional_pval$p.val_log10 <- ifelse(trans_regional_pval$p.val_numeric > 0, 
                                        -log10(trans_regional_pval$p.val_numeric), 
                                        0)

# Restore the sign to the log-transformed values
trans_regional_pval$p.val_log10 <- ifelse(grepl("-", trans_regional_pval$p.val), 
                                        -trans_regional_pval$p.val_log10, 
                                        trans_regional_pval$p.val_log10)

# Keep only the p.val_log10 column
trans_regional_pval <- trans_regional_pval["p.val_log10"]
trans_regional_pval <- as.data.frame(trans_regional_pval[-c(1:32),])

write.table(trans_regional_pval, file = "/projects/tsecara/SPINS_ASD_Project2/data/visualizations/P-scalers/SSD_ASD_regional_EA_adj.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)


##################################################RESTING STATE #####################################################
###Now repeat for the resting state data 
mssd_val <- read.csv("/projects/tsecara/SPINS_ASD_Project2/data/combined/regional_MSSD/neurocombat/combined_mssd_mean_RS_FINAL_agereduced.csv")
data <- read.csv("/projects/tsecara/SPINS_ASD_Project2/data/combined/demographics/demo_RS_agereduced.csv")
data <- data[(data$record_id %in% mssd_val$record_id), ]

region_merged <- merge(data, mssd_val, by = 'record_id')

# Define function to perform ANOVA and TukeyHSD
perform_anova <- function(df, response_col, predictor_cols, tukey_if_significant = TRUE) {
  # Build the formula for the ANOVA model
  formula <- as.formula(paste(response_col, '~', paste(predictor_cols, collapse = ' + ')))
  
  # Run the ANOVA model
  model <- aov(formula, data = df)
  
  # Get the ANOVA table using car::Anova
  anova_results <- car::Anova(model, type = 'II')
  
  # Check if 'group' is significant
  group_significant <- anova_results$'Pr(>F)'[rownames(anova_results) == 'group'] < 0.05
  
  # If 'group' is significant, run TukeyHSD and store results
  if (group_significant & tukey_if_significant) {
    tukey_results <- TukeyHSD(aov(as.formula(paste(response_col, '~ group')), data = df))
    tukey_results <- tidy(tukey_results) %>%
      mutate(response_col = response_col)
  } else {
    tukey_results <- NULL
  }
  
  # Return both ANOVA and TukeyHSD results
  list(anova_results = tidy(anova_results) %>% mutate(response_col = response_col),
       tukey_results = tukey_results)
}

# Define the predictor columns
predictor_cols <- c('group', 'age', 'sex', 'avg_fd')

# Use pbapply to include a progress bar and process networks one by one
results_list <- pblapply(names(mssd_val)[-1], function(region_col) {
  perform_anova(region_merged, region_col, predictor_cols)
})

# Combine results into a final dataframe
anova_df <- bind_rows(lapply(results_list, `[[`, 'anova_results'))
tukey_df <- bind_rows(lapply(results_list, `[[`, 'tukey_results'))

# Apply FDR correction using the Benjamini-Hochberg method
tukey_df$fdr_adj_p.value <- p.adjust(tukey_df$adj.p.value, method = "fdr")

#Filter out the values:
filtered_tukey_df <- tukey_df[tukey_df$fdr_adj_p.value < 0.025, ]
filtered_tukey_df_ASD <- filtered_tukey_df[filtered_tukey_df$contrast == "Control-ASD",]
filtered_tukey_df_SSD <- filtered_tukey_df[filtered_tukey_df$contrast == "SSD-Control",]
filtered_tukey_df_trans <- filtered_tukey_df[filtered_tukey_df$contrast == "SSD-ASD",]

####### Creating p-scaler files:
#####ASD First 
####### Postive values indicate the regions where Controls have higher MSSD than ASD, and Negative values indicuate where the ASD has higher MSSD values than controls 
# Extract column names excluding "record_id"
column_names <- colnames(mssd_val)[colnames(mssd_val) != "record_id"]

# Create a new dataframe with one column titled "p-val"
ASD_regional_pval <- data.frame(p.val = rep(NA, length(column_names)))
row.names(ASD_regional_pval) <- column_names

# Convert p.val column to character to accommodate signs
ASD_regional_pval$p.val <- as.character(ASD_regional_pval$p.val)

# Update ASD_regional_pval with values from filtered_tukey_df_ASD
for (i in seq_len(nrow(filtered_tukey_df_ASD))) {
  response <- filtered_tukey_df_ASD$response_col[i]
  fdr_value <- filtered_tukey_df_ASD$fdr_adj_p.value[i]
  estimate <- filtered_tukey_df_ASD$estimate[i]
  
  # Determine the sign based on the estimate
  sign <- ifelse(estimate > 0, "+", "-")
  
  # Match and update
  if (response %in% rownames(ASD_regional_pval)) {
    ASD_regional_pval[response, "p.val"] <- paste0(sign, sprintf("%.10f", fdr_value))
  }
}

# Set unmatched rows to "0"
ASD_regional_pval[is.na(ASD_regional_pval$p.val), "p.val"] <- "0"

# Print the updated structure to verify
str(ASD_regional_pval)
ASD_regional_pval$p.val <- as.numeric(ASD_regional_pval$p.val)
range(ASD_regional_pval$p.val)

# Convert the p.val column to numeric, ignoring signs for now
ASD_regional_pval$p.val_numeric <- as.numeric(gsub("[+-]", "", ASD_regional_pval$p.val))

# Apply log10 transformation, but only to non-zero values
ASD_regional_pval$p.val_log10 <- ifelse(ASD_regional_pval$p.val_numeric > 0, 
                                        -log10(ASD_regional_pval$p.val_numeric), 
                                        0)

# Restore the sign to the log-transformed values
ASD_regional_pval$p.val_log10 <- ifelse(grepl("-", ASD_regional_pval$p.val), 
                                        -ASD_regional_pval$p.val_log10, 
                                        ASD_regional_pval$p.val_log10)

# Keep only the p.val_log10 column
ASD_regional_pval <- ASD_regional_pval["p.val_log10"]
ASD_regional_pval <- as.data.frame(ASD_regional_pval[-c(1:32),])

write.table(ASD_regional_pval, file = "/projects/tsecara/SPINS_ASD_Project2/data/visualizations/P-scalers/ASD_TDC_regional_RS_adj.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

#####SSD VS Control 
####### Postive values indicate the regions where Controls have higher MSSD than SSD, and Negative values indicuate where the SSDs has higher MSSD values than controls 
# Extract column names excluding "record_id"
column_names <- colnames(mssd_val)[colnames(mssd_val) != "record_id"]

# Create a new dataframe with one column titled "p-val"
SSD_regional_pval <- data.frame(p.val = rep(NA, length(column_names)))
row.names(SSD_regional_pval) <- column_names

# Convert p.val column to character to accommodate signs
SSD_regional_pval$p.val <- as.character(SSD_regional_pval$p.val)

# Update SSD_regional_pval with values from filtered_tukey_df_SSD
for (i in seq_len(nrow(filtered_tukey_df_SSD))) {
  response <- filtered_tukey_df_SSD$response_col[i]
  fdr_value <- filtered_tukey_df_SSD$fdr_adj_p.value[i]
  estimate <- filtered_tukey_df_SSD$estimate[i]
  
  # Determine the sign based on the estimate
  sign <- ifelse(estimate > 0, "-", "+")
  
  # Match and update
  if (response %in% rownames(SSD_regional_pval)) {
    SSD_regional_pval[response, "p.val"] <- paste0(sign, sprintf("%.20f", fdr_value))
  }
}

# Set unmatched rows to "0"
SSD_regional_pval[is.na(SSD_regional_pval$p.val), "p.val"] <- "0"

#Changing back to numeric 
SSD_regional_pval$p.val <- as.numeric(SSD_regional_pval$p.val)
range(SSD_regional_pval$p.val)

# Convert the p.val column to numeric, ignoring signs for now
SSD_regional_pval$p.val_numeric <- as.numeric(gsub("[+-]", "", SSD_regional_pval$p.val))

# Apply log10 transformation, but only to non-zero values
SSD_regional_pval$p.val_log10 <- ifelse(SSD_regional_pval$p.val_numeric > 0, 
                                        -log10(SSD_regional_pval$p.val_numeric), 
                                        0)

# Restore the sign to the log-transformed values
SSD_regional_pval$p.val_log10 <- ifelse(grepl("-", SSD_regional_pval$p.val), 
                                        -SSD_regional_pval$p.val_log10, 
                                        SSD_regional_pval$p.val_log10)

# Keep only the p.val_log10 column
SSD_regional_pval <- SSD_regional_pval["p.val_log10"]
SSD_regional_pval <- as.data.frame(SSD_regional_pval[-c(1:32),])

write.table(SSD_regional_pval, file = "/projects/tsecara/SPINS_ASD_Project2/data/visualizations/P-scalers/SSD_TDC_regional_RS_adj.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

#####SSD vs ASD
####### Postive values indicate the regions where SSD have higher MSSD than ASD, and Negative values indicuate where the ASD has higher MSSD values than SSD 
# Extract column names excluding "record_id"
column_names <- colnames(mssd_val)[colnames(mssd_val) != "record_id"]

# Create a new dataframe with one column titled "p-val"
trans_regional_pval <- data.frame(p.val = rep(NA, length(column_names)))
row.names(trans_regional_pval) <- column_names

# Convert p.val column to character to accommodate signs
trans_regional_pval$p.val <- as.character(trans_regional_pval$p.val)

# Update trans_regional_pval with values from filtered_tukey_df_trans
for (i in seq_len(nrow(filtered_tukey_df_trans))) {
  response <- filtered_tukey_df_trans$response_col[i]
  fdr_value <- filtered_tukey_df_trans$fdr_adj_p.value[i]
  estimate <- filtered_tukey_df_trans$estimate[i]
  
  # Determine the sign based on the estimate
  sign <- ifelse(estimate > 0, "+", "-")
  
  # Match and update
  if (response %in% rownames(trans_regional_pval)) {
    trans_regional_pval[response, "p.val"] <- paste0(sign, sprintf("%.10f", fdr_value))
  }
}

# Set unmatched rows to "0"
trans_regional_pval[is.na(trans_regional_pval$p.val), "p.val"] <- "0"

trans_regional_pval$p.val <- as.numeric(trans_regional_pval$p.val)
range(trans_regional_pval$p.val)

# Convert the p.val column to numeric, ignoring signs for now
trans_regional_pval$p.val_numeric <- as.numeric(gsub("[+-]", "", trans_regional_pval$p.val))

# Apply log10 transformation, but only to non-zero values
trans_regional_pval$p.val_log10 <- ifelse(trans_regional_pval$p.val_numeric > 0, 
                                          -log10(trans_regional_pval$p.val_numeric), 
                                          0)

# Restore the sign to the log-transformed values
trans_regional_pval$p.val_log10 <- ifelse(grepl("-", trans_regional_pval$p.val), 
                                          -trans_regional_pval$p.val_log10, 
                                          trans_regional_pval$p.val_log10)

# Keep only the p.val_log10 column
trans_regional_pval <- trans_regional_pval["p.val_log10"]
trans_regional_pval <- as.data.frame(trans_regional_pval[-c(1:32),])

write.table(trans_regional_pval, file = "/projects/tsecara/SPINS_ASD_Project2/data/visualizations/P-scalers/SSD_ASD_regional_RS_adj.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)