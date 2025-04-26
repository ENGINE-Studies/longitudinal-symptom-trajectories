require(reshape2)
require(tidyverse)
require(factoextra)
install.packages("reshape2")
install.packages("dplyr")

# Load the package
library(reshape2)
library(dplyr)

# Install tidyr if not already installed
install.packages("tidyr")

# Load the tidyr package
library(tidyr)

install.packages("tidyverse")

# Load tidyverse, which includes tidyr (or load tidyr directly with library(tidyr))
library(tidyverse)

eeg <- read.csv('/Users/wadimvodovozov/Desktop/ENGINE/analysis/master_connectivity_SPARSE_30-Oct-2024.csv')
eeg2 <- read.csv('//Users/wadimvodovozov/Desktop/ENGINE/analysis/master_connectivity_SPARSE_17-Nov-2024.csv')

# Filter and select required columns
eeg_filtered <- eeg2 %>%
  filter(Population == 'MDD') %>%
  select(Subject, Session, Channel_1, Channel_2, FreqBand, Coherence, Imag_Coherence, WPLI, WPPC)

# Filter subjects with data for all three timepoints
eeg_filtered_complete <- eeg_filtered %>%
  group_by(Subject) %>%
  filter(all(c(1, 2, 3) %in% Session)) %>%  # Ensure data is present for sessions 1, 2, and 3
  ungroup()

# Check how many subjects you have now
length(unique(eeg_filtered_complete$Subject))

# Melt to long format, creating "variable" column for measurement types
m <- melt(eeg_filtered_complete, id.vars = c('Subject', 'Session', 'Channel_1', 'Channel_2', 'FreqBand'))

# Check structure of `m` to ensure it includes all variables
str(m)
head(m) # Inspect to confirm it includes Coherence, Imag_Coherence, WPLI, WPPC

# Reshape to wide format
eeg_wide <- m %>%
  pivot_wider(
    names_from = c(Channel_1, Channel_2, FreqBand, variable),
    values_from = value
  )

## Troubleshooting TBD: pivot wider produced NAs; dropping these for now but figure out why this happens
#na_idx <- which(is.na(eeg_wide), arr.ind = T)
#eeg_wide <- eeg_wide[, -unique(na_idx[, 2])]

# split by timepoint
eeg_wide_by_tp <- split(eeg_wide, eeg_wide$Session)

# Compute PCA for each timepoint (TP1, TP2, TP3)
eeg_names <- grep('Coherence|Imag_Coherence|WPLI|WPPC', names(eeg_wide_by_tp[[1]]), value = TRUE)

# PCA for each timepoint
pc_tp1 <- prcomp(eeg_wide_by_tp[[1]][, eeg_names], scale = TRUE)
pc_tp2 <- predict(pc_tp1, newdata = eeg_wide_by_tp[[2]][, eeg_names])
pc_tp3 <- predict(pc_tp1, newdata = eeg_wide_by_tp[[3]][, eeg_names])

pc_tp1 <- prcomp(eeg_wide_by_tp[[1]][, eeg_names], scale = TRUE)
pc_tp2 <- prcomp(eeg_wide_by_tp[[2]][, eeg_names], scale = TRUE)
pc_tp3 <- prcomp(eeg_wide_by_tp[[3]][, eeg_names], scale = TRUE)

pc_tp1$rotation
pc_tp2$rotation
pc_tp3$rotation

# Check the structure of the result to ensure that `x` is present
str(pc_tp2)

# Now access the PCA scores (the `x` component)
pc_tp2_scores <- pc_tp2$x  # These are the PCA scores for the data in session 2


# Reshape PCA data for each timepoint into long format
# For TP1
pc_tp1_long <- data.frame(
  Subject = rep(eeg_wide_by_tp[[1]]$Subject, each = ncol(pc_tp1$x)),  # Repeat Subject for each PC
  Session = rep(1, nrow(eeg_wide_by_tp[[1]]) * ncol(pc_tp1$x)),  # Session 1
  PC = rep(paste0("PC", 1:ncol(pc_tp1$x)), times = nrow(eeg_wide_by_tp[[1]])),  # Repeat PC labels for each participant
  Value = as.vector(pc_tp1$x)  # PC values
)

# For TP2
pc_tp2_long <- data.frame(
  Subject = rep(eeg_wide_by_tp[[2]]$Subject, each = ncol(pc_tp2)),  # Repeat Subject for each PC
  Session = rep(2, nrow(eeg_wide_by_tp[[2]]) * ncol(pc_tp2)),  # Session 2
  PC = rep(paste0("PC", 1:ncol(pc_tp2)), times = nrow(eeg_wide_by_tp[[2]])),  # Repeat PC labels for each participant
  Value = as.vector(pc_tp2)  # PC values
)

# For TP3
pc_tp3_long <- data.frame(
  Subject = rep(eeg_wide_by_tp[[3]]$Subject, each = ncol(pc_tp3)),  # Repeat Subject for each PC
  Session = rep(3, nrow(eeg_wide_by_tp[[3]]) * ncol(pc_tp3)),  # Session 3
  PC = rep(paste0("PC", 1:ncol(pc_tp3)), times = nrow(eeg_wide_by_tp[[3]])),  # Repeat PC labels for each participant
  Value = as.vector(pc_tp3)  # PC values
)

# Combine pc_tp1_long, pc_tp2_long, and pc_tp3_long sequentially for each subject and session
pc_long <- rbind(
  pc_tp1_long %>% mutate(Session = 1),  # Ensure session 1 is labeled properly
  pc_tp2_long %>% mutate(Session = 2),  # Ensure session 2 is labeled properly
  pc_tp3_long %>% mutate(Session = 3)   # Ensure session 3 is labeled properly
)

# Arrange by Subject and Session, ensuring the data is ordered by participant and session
pc_long <- pc_long %>%
  arrange(Subject, Session, PC)

# Load the dplyr package if not already loaded
library(dplyr)

# Remove rows where Subject is 'sub-3043' or 'sub-3045'
pc_long <- pc_long %>%
  filter(!Subject %in% c('sub-3043', 'sub-3045'))

# View the final long-format data
head(pc_long)

# Ensure 'PC' column is a factor with correct levels
pc_long$PC <- factor(pc_long$PC, levels = paste0("PC", 1:23), ordered = TRUE)

# Check if the order is correct before merging
head(pc_long)

# Create age data for each subject
age_data <- data.frame(
  Subject = unique(pc_long$Subject),
  age = c(16, 17, 14, 16, 17, 15, 17, 15, 17, 14, 16, 16, 17, 14, 17, 15, 16, 14, 17, 15, 15)
)

# First merge PC data with ratings_long on 'Subject' and 'Session'
pc_ratings_merged <- merge(pc_long, ratings_long, by = c("Subject", "Session"))

# Now merge with age_data to include age information
pc_merged_data <- merge(pc_ratings_merged, age_data, by = "Subject")


# Check the structure of the merged data to ensure correct merge
head(pc_merged_data)

# Now, we sort by Subject, Session, and PC to ensure correct order
pc_merged_data <- pc_merged_data %>%
  arrange(Subject, Session, PC)

# Check if the rows are sorted correctly
head(pc_merged_data)
tail(pc_merged_data)

# Optionally, check the full structure of the data
str(pc_merged_data)

# Apply lme
# Load required libraries
library(lme4)
library(glmmLasso)

# Fit the linear mixed-effects model
# Create the formula dynamically for PCs and interactions
pc_vars <- paste0("PC", 1:23)  # This will create a vector of PC1 to PC23

library(tidyr)

# Reshape the data to wide format if necessary
pc_merged_data_wide <- pc_merged_data %>%
  pivot_wider(
    names_from = PC, 
    values_from = Value
  )

# Create a vector of PC variable names (PC1 to PC23)
pc_vars <- paste0("PC", 1:23)

# Create interaction terms for each PC with Session
interaction_terms <- paste0(pc_vars, " * Session")

# Include PC variables, age, Session, and the interaction terms in the model
fixed_effects <- c(pc_vars, "age", "Session", interaction_terms)

# Construct the formula for the linear mixed-effects model
formula <- as.formula(paste("cdrs_total ~", paste(fixed_effects, collapse = " + "), "+ (1 | Subject)"))

library(lme4)

# Fit the linear mixed-effects model
model <- lmer(
  formula,
  data = pc_merged_data_wide,
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 10000))
)

# Check the model summary
summary(model)

# Capture the summary output
summary_output <- capture.output(summary(model))

# Specify the file path and name
file_path <- "new_PCA_model_summary.txt"

# Write the summary to a text file
writeLines(summary_output, file_path)

print(paste("Model summary has been written to:", file_path))

## USE LASSO MODEL

# Fit LASSO model using glmmLasso
pc_merged_data$Subject <- as.factor(pc_merged_data$Subject)

lasso_model <- glmmLasso(
  formula,
  data = pc_merged_data,
  lambda = 0.1
)

summary(lasso_model)

library(glmnet)
# create variables for models
x <- model.matrix(cdrs_total ~ PC1 + PC2 + PC3 + PC4 + age + Session, data = pc_merged_data)[, -1]
y <- pc_merged_data$cdrs_total

# Lasso 
fit_lasso = glmnet(x, y, alpha = 1)

# Ridge
fit_ridge = glmnet(x, y, alpha = 0)

# Elastic net
fit_enet = glmnet(x, y, alpha = 0.5)

# Display the model summary
print(fit_lasso)

# Extract coefficients
coef_matrix = coef(fit_lasso)

# Perform cross validation to select otimal lambda
cv_lasso = cv.glmnet(x, y, alpha = 1)
plot(cv_lasso)

# Get optimal lambda values (indicative of strength of penalty applied)
best_lambda = cv_lasso$lambda.min  # lambda that gives minimum cross-validation error
best_lambda_1se = cv_lasso$lambda.1se  # largest lambda within 1 standard error of minimum

# Extract coefficients for the best model
best_coef = coef(cv_lasso, s = "lambda.min")

# Make predictions based on chosen lambda
predictions = predict(fit_lasso, newx = x, s = best_lambda)

# Add predictions to the original dataset
pc_merged_data$pred_cdrs = as.vector(predictions)
pc_merged_data_wide$pred_cdrs = as.vector(predictions)

# Now you can compare actual vs predicted
head(pc_merged_data[, c("cdrs_total", "pred_cdrs")])

# Calculate correlation between actual and predicted
correlation = cor(pc_merged_data$cdrs_total, pc_merged_data$pred_cdrs)
print(paste("Correlation between actual and predicted cdrs_total:", round(correlation, 3)))

# Calculate Mean Absolute Error (MAE)
mae = mean(abs(pc_merged_data$cdrs_total - pc_merged_data$pred_cdrs))
print(paste("Mean Absolute Error:", round(mae, 3)))

# Calculate Root Mean Squared Error (RMSE)
rmse = sqrt(mean((pc_merged_data$cdrs_total - pc_merged_data$pred_cdrs)^2))
print(paste("Root Mean Squared Error:", round(rmse, 3)))

# Create a scatter plot of actual vs predicted values
plot(pc_merged_data$cdrs_total, pc_merged_data$pred_cdrs, 
     xlab = "Actual cdrs_total", ylab = "Predicted cdrs_total",
     main = "Actual vs Predicted cdrs_total")
abline(a = 0, b = 1, col = "red")  # Add a 45-degree line for reference

# Specify the file path and name
file_path <- "new_pc_merged_data.csv"

# Write the data frame to a CSV file
write.csv(pc_merged_data, file = file_path, row.names = FALSE)

print(paste("Data has been written to:", file_path))

correlation <- cor(pc_merged_data$cdrs_total, pc_merged_data$pred_cdrs)
print(paste("Correlation between actual and predicted cdrs_total:", round(correlation, 3)))

# Calculate R-squared
r_squared <- correlation^2
print(paste("R-squared:", round(r_squared, 3)))


# Function to extract top 5 loadings for each component from a PCA rotation matrix
get_top_loadings <- function(rotation_matrix, n_top = 5) {
  # Initialize list to store top loadings for each component
  top_loadings_list <- list()
  
  # Loop through each component in the rotation matrix
  for (component in 1:ncol(rotation_matrix)) {
    # Extract loadings for the current component
    loadings <- rotation_matrix[, component]
    
    # Get the absolute values of loadings and sort them in descending order
    top_indices <- order(abs(loadings), decreasing = TRUE)[1:n_top]
    
    # Store the top n loadings and the corresponding channel pairs
    top_loadings_list[[paste0("Component_", component)]] <- data.frame(
      Channel_Pair = rownames(rotation_matrix)[top_indices],
      Loading_Value = loadings[top_indices]
    )
  }
  return(top_loadings_list)
}

# Apply the function to each time point
top_loadings_tp1 <- get_top_loadings(pc_tp1$rotation)
top_loadings_tp2 <- get_top_loadings(pc_tp2$rotation)
top_loadings_tp3 <- get_top_loadings(pc_tp3$rotation)

# Print the results for each time point
print("Top 5 Loadings - Time Point 1")
top_loadings_tp1

print("Top 5 Loadings - Time Point 2")
top_loadings_tp2

print("Top 5 Loadings - Time Point 3")
top_loadings_tp3

# Function to save top loadings as a text file
save_top_loadings_to_txt <- function(top_loadings_list, file_prefix) {
  # Loop through each component and save the results as text files
  for (component_name in names(top_loadings_list)) {
    # Extract the top loadings data frame for the current component
    component_data <- top_loadings_list[[component_name]]
    
    # Define the file name based on the component and time point
    filename <- paste0(file_prefix, "_", component_name, ".txt")
    
    # Write the data frame to a text file
    write.table(component_data, file = filename, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  }
}

# Apply the function to save the results for each time point
save_top_loadings_to_txt(top_loadings_tp1, "top_loadings_tp1")
save_top_loadings_to_txt(top_loadings_tp2, "top_loadings_tp2")
save_top_loadings_to_txt(top_loadings_tp3, "top_loadings_tp3")

# Message indicating the files have been saved
cat("Top loadings have been saved as .txt files.\n")

# Function to save the strongest loadings plot as a PNG file
save_strongest_loadings_plot_base <- function(strongest_loadings, title, filename) {
  # Open a PNG device to save the plot
  png(filename, width = 800, height = 600)  # Adjust size as needed
  
  # Create the barplot of loading values
  barplot(
    strongest_loadings$Loading_Value,
    names.arg = strongest_loadings$Channel_Pair,
    main = title,
    las = 2,  # Rotate x-axis labels for readability
    col = "lightblue",
    border = "black",
    cex.names = 0.8,
    ylab = "Loading Value",
    xlab = "Channel Pair"
  )
  
  # Close the PNG device
  dev.off()
}

# Save the strongest loadings for each time point as PNG files
save_strongest_loadings_plot_base(strongest_loadings_tp1, "Top Strongest Loadings - Time Point 1", "strongest_loadings_tp1.png")
save_strongest_loadings_plot_base(strongest_loadings_tp2, "Top Strongest Loadings - Time Point 2", "strongest_loadings_tp2.png")
save_strongest_loadings_plot_base(strongest_loadings_tp3, "Top Strongest Loadings - Time Point 3", "strongest_loadings_tp3.png")


