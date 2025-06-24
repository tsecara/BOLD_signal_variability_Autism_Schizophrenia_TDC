# Load necessary packages
library(car)
library(dplyr)
library(broom)
library(pbapply)

#Calculating network variability effect sizes 
#First for EA task 
mssd_val_network <- read.csv("/projects/tsecara/SPINS_ASD_Project2/data/combined/network_MSSD/EA_MSSD_NETWORK_agereduced.csv")
data <- read.csv("/projects/tsecara/SPINS_ASD_Project2/data/combined/demographics/demo_EA_agereduced.csv")
networked_merged <- merge(data, mssd_val_network, by = 'record_id')
networked_merged <- networked_merged[networked_merged$age <= 35, ]
networked_merged$avg_fd <- as.numeric(networked_merged$avg_fd)

#For resting state
data <- read.csv("/projects/tsecara/SPINS_ASD_Project2/data/combined/demographics/demo_RS_agereduced.csv")
mssd_val_network <- read.csv("/projects/tsecara/SPINS_ASD_Project2/data/combined/network_MSSD/RS_MSSD_NETWORK_agereduced.csv")
networked_merged <- merge(data, mssd_val_network, by = 'record_id')
networked_merged <- networked_merged[networked_merged$age <= 35, ]
networked_merged$avg_fd <- as.numeric(networked_merged$avg_fd)

# Function to compute network group differences & Cohen's d between two groups
compute_cohens_d <- function(df, response_col, group1, group2, group_col = "group") {
  g1_vals <- df[[response_col]][df[[group_col]] == group1]
  g2_vals <- df[[response_col]][df[[group_col]] == group2]
  
  g1_vals <- g1_vals[!is.na(g1_vals)]
  g2_vals <- g2_vals[!is.na(g2_vals)]
  
  n1 <- length(g1_vals)
  n2 <- length(g2_vals)
  sd1 <- sd(g1_vals)
  sd2 <- sd(g2_vals)
  
  pooled_sd <- sqrt(((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2))
  d <- (mean(g1_vals) - mean(g2_vals)) / pooled_sd
  return(d)
}

# Function to perform ANCOVA and post hoc TukeyHSD
perform_anova <- function(df, response_col, predictor_cols, tukey_if_significant = TRUE) {
  # Build formula and run model
  formula <- as.formula(paste(response_col, '~', paste(predictor_cols, collapse = ' + ')))
  model <- aov(formula, data = df)
  
  # Get ANOVA results
  anova_results <- car::Anova(model, type = 'II')
  group_significant <- anova_results$'Pr(>F)'[rownames(anova_results) == 'group'] < 0.05
  
  # Run Tukey HSD if group is significant
  if (group_significant & tukey_if_significant) {
    tukey_raw <- TukeyHSD(aov(as.formula(paste(response_col, '~ group')), data = df))
    tukey_results <- tidy(tukey_raw) %>%
      mutate(response_col = response_col) %>%
      mutate(comparison = rownames(tukey_raw$group)) %>%
      rowwise() %>%
      mutate(
        group1 = strsplit(comparison, "-")[[1]][1],
        group2 = strsplit(comparison, "-")[[1]][2],
        cohens_d = compute_cohens_d(df, response_col, group1, group2)
      )
  } else {
    tukey_results <- NULL
  }
  
  # Return tidy ANOVA and Tukey tables
  list(
    anova_results = tidy(anova_results) %>% mutate(response_col = response_col),
    tukey_results = tukey_results
  )
}

# ----------------------------
# Run full analysis
# Replace with your actual dataframe names:
# - mssd_val_network: includes record_id + network columns
# - networked_merged: includes merged network + group/covariates

# Define predictors
predictor_cols <- c('group', 'age', 'sex', 'avg_fd')

# Run analysis for each network column (skip record_id)
results_list <- pblapply(names(mssd_val_network)[-1], function(network_col) {
  perform_anova(networked_merged, network_col, predictor_cols)
})

# Combine ANOVA and post hoc results
anova_df <- bind_rows(lapply(results_list, `[[`, 'anova_results'))
tukey_df <- bind_rows(lapply(results_list, `[[`, 'tukey_results'))

# Apply FDR correction to post hoc p-values
tukey_df$fdr_adj_p.value <- p.adjust(tukey_df$adj.p.value, method = "fdr")

# Filter significant post hoc results
filtered_tukey_df <- tukey_df %>% filter(fdr_adj_p.value < 0.05)

# Print results
print(anova_df)
print(filtered_tukey_df) #this contains signficant networks with effect size
