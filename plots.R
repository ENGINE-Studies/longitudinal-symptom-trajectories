# Forest Plot of PC coefficient from LME model

# Load required libraries
library(ggplot2)
library(broom.mixed)
library(dplyr)

# Extract fixed effects and confidence intervals from the model
model_fixed <- broom.mixed::tidy(model, effects = "fixed", conf.int = TRUE)

# Filter for only PC coefficients (not interaction terms or other covariates)
pc_effects <- model_fixed %>%
  filter(grepl("^PC\\d+$", term)) %>%
  arrange(estimate) %>%
  mutate(term = factor(term, levels = term))  # order for plotting

# Plot using ggplot2
ggplot(pc_effects, aes(x = estimate, y = term)) +
  geom_point() +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Forest Plot of PC Coefficients from LME Model",
       x = "Coefficient Estimate",
       y = "Principal Component") +
  theme_minimal(base_size = 14)


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
#top_loadings_tp2 <- get_top_loadings(pc_tp2$rotation)
#top_loadings_tp3 <- get_top_loadings(pc_tp3$rotation)

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
#save_top_loadings_to_txt(top_loadings_tp2, "top_loadings_tp2")
#save_top_loadings_to_txt(top_loadings_tp3, "top_loadings_tp3")

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

# Load required libraries
library(dplyr)
library(tidyr)
library(ggplot2)

# Step 1: Reshape and clean data
cdrs_long <- cdrs_data %>%
  pivot_longer(
    cols = c(CDRS_total_1, CDRS_total_2, CDRS_total_3),
    names_to = "Session",
    values_to = "CDRS_Total"
  ) %>%
  mutate(
    Session = recode(Session,
                     CDRS_total_1 = "Baseline",
                     CDRS_total_2 = "Week 4",
                     CDRS_total_3 = "Week 16")
  )

# Ensure consistent factor levels (avoids duplication of sessions)
cdrs_long$Session <- factor(cdrs_long$Session,
                            levels = c("Baseline", "Week 4", "Week 16"))

# Step 2: Calculate group means and standard errors
cdrs_summary <- cdrs_long %>%
  group_by(Session) %>%
  summarise(
    mean_cdrs = mean(CDRS_Total, na.rm = TRUE),
    se_cdrs = sd(CDRS_Total, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

# Step 3: Create plot
cdrs_plot <- ggplot(cdrs_long, aes(x = Session, y = CDRS_Total, group = Subject)) +
  geom_line(aes(color = Subject), alpha = 0.4) +
  geom_point(aes(color = Subject), size = 2) +
  
  # Group mean line
  geom_line(data = cdrs_summary,
            aes(x = Session, y = mean_cdrs, group = 1),
            color = "black", linewidth = 1.2,
            inherit.aes = FALSE) +
  
  # Error bars
  geom_errorbar(data = cdrs_summary,
                aes(x = Session,
                    ymin = mean_cdrs - se_cdrs,
                    ymax = mean_cdrs + se_cdrs),
                width = 0.2,
                color = "black", linewidth = 0.8,
                inherit.aes = FALSE) +
  
  # Label above Week 16 mean + SE
  geom_text(data = cdrs_summary %>% filter(Session == "Week 16"),
            aes(x = 3.35, y = 34, label = "Mean CDRS ± SE"),
            inherit.aes = FALSE,
            color = "black", fontface = "italic", size = 3) +
  
  labs(title = "CDRS Total Scores Over Time",
       x = "Session",
       y = "CDRS Total Score") +
  
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  )

# Step 4: Print the plot
cdrs_plot


# Save the plot as JPEG
ggsave("/Users/wadimvodovozov/Desktop/ENGINE/Figures/NEW_CDRS_total_over_time.jpg",
       plot = cdrs_plot,
       width = 8, height = 6, dpi = 300,
       device = "jpeg")

## run LME on CDRS scores over time only

library(lme4)
library(lmerTest)

# Convert data to long format
cdrs_long <- cdrs_data %>%
  pivot_longer(cols = starts_with("CDRS_total"),
               names_to = "Session",
               values_to = "CDRS_Total") %>%
  mutate(Session = recode(Session,
                          CDRS_total_1 = "Baseline",
                          CDRS_total_2 = "Week4",
                          CDRS_total_3 = "Week16"),
         Session = factor(Session, levels = c("Baseline", "Week4", "Week16")))

# Fit linear mixed-effects model
lmm <- lmer(CDRS_Total ~ Session + (1 | Subject), data = cdrs_long)
summary(lmm)

# Load emmeans if not already loaded
install.packages("emmeans")  # only if not yet installed
library(emmeans)

# Get estimated marginal means
emm <- emmeans(lmm, ~ Session)

# Contrast Week 16 vs Week 4
contrast(emm, method = list("Week16 vs Week4" = c(0, -1, 1)))

library(lme4)
library(broom.mixed)
library(ggplot2)
library(dplyr)

# Extract fixed effects with confidence intervals
model_fixed <- broom.mixed::tidy(model, effects = "fixed", conf.int = TRUE)

# Filter only interaction terms (e.g., PC21:SessionWeek4)
interaction_effects <- model_fixed %>%
  filter(grepl("^PC(21|23):Session", term)) %>%
  mutate(Significant = TRUE) %>%
  bind_rows(
    model_fixed %>%
      filter(grepl("^PC\\d+:Session", term) & !grepl("^PC(21|23):Session", term)) %>%
      mutate(Significant = FALSE)
  ) %>%
  arrange(estimate) %>%
  mutate(term = factor(term, levels = term))  # Keep order for plotting

# Create the forest plot
interaction_forest_plot <- ggplot(interaction_effects, aes(x = estimate, y = term, color = Significant)) +
  geom_point(size = 2) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  scale_color_manual(values = c("black", "blue")) +
  labs(title = "Forest Plot of PC × Session Interaction Terms",
       x = "Coefficient Estimate",
       y = "Interaction Term") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")  # Hide legend for clean plot

# Save the plot
ggsave(
  filename = "/Users/wadimvodovozov/Desktop/ENGINE/Figures/PC_interaction_forest_plot_highlighted.jpg",
  plot = interaction_forest_plot,
  width = 8,
  height = 6,
  dpi = 300
)

library(tidyverse)

# Extract loadings for PC21 and PC23
loadings_df <- as.data.frame(pc_tp1$rotation[, c("PC21", "PC23")])
loadings_df$Variable <- rownames(loadings_df)

# Create a new column that standardizes channel pair order (e.g., Fz_Pz and Pz_Fz → Fz_Pz)
loadings_df <- loadings_df %>%
  mutate(Dyad = Variable) %>%
  separate(Dyad, into = c("Ch1", "Ch2", "Freq", "Measure"), sep = "_", fill = "right", remove = FALSE) %>%
  mutate(
    Ch_pair = if_else(Ch1 < Ch2, paste(Ch1, Ch2, sep = "_"), paste(Ch2, Ch1, sep = "_")),
    UniqueID = paste(Ch_pair, Freq, Measure, sep = "_")
  ) %>%
  distinct(UniqueID, .keep_all = TRUE)  # Keep only one version of each dyad

# Select top 5 absolute loadings per PC from unique dyads
top_pc21 <- loadings_df %>%
  slice_max(order_by = abs(PC21), n = 5) %>%
  select(Variable, Loading = PC21) %>%
  mutate(PC = "PC21")

top_pc23 <- loadings_df %>%
  slice_max(order_by = abs(PC23), n = 5) %>%
  select(Variable, Loading = PC23) %>%
  mutate(PC = "PC23")

# Combine and reformat for plotting
top_loadings <- bind_rows(top_pc21, top_pc23)

# Plot
pc_loading_plot <- ggplot(top_loadings, aes(x = reorder(Variable, Loading), y = Loading, fill = PC)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  labs(title = "Top 5 Loadings for PC21 and PC23",
       x = "Original Variable",
       y = "Loading Coefficient") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),  # centered, non-bold title
    legend.title = element_blank()
  )

# Save plot as JPEG
ggsave(
  filename = "/Users/wadimvodovozov/Desktop/ENGINE/Figures/Top5_Unique_PC21_PC23_Loadings.jpg",
  plot = pc_loading_plot,
  width = 8,
  height = 6,
  dpi = 300
)


