# Load required libraries
library(tidyverse)
library(reshape2)
library(dplyr)
library(tidyr)

# Read EEG file
eeg2 <- read.csv('//Users/wadimvodovozov/Desktop/ENGINE/analysis/master_connectivity_SPARSE_17-Nov-2024.csv')

# --- Prepare MDD Dataset ---
eeg_mdd <- eeg2 %>%
  filter(Population == 'MDD') %>%
  select(Subject, Session, Channel_1, Channel_2, FreqBand, Coherence, Imag_Coherence, WPLI, WPPC)

# Only include MDD participants with Sessions 1, 2, and 3
eeg_mdd_complete <- eeg_mdd %>%
  group_by(Subject) %>%
  filter(all(c(1, 2, 3) %in% Session)) %>%
  ungroup()

# Melt MDD data into long format
m_mdd <- melt(eeg_mdd_complete, id.vars = c('Subject', 'Session', 'Channel_1', 'Channel_2', 'FreqBand'))

# Reshape to wide format for MDD
eeg_mdd_wide <- m_mdd %>%
  pivot_wider(
    names_from = c(Channel_1, Channel_2, FreqBand, variable),
    values_from = value
  )

# Split by Session
eeg_mdd_wide_by_tp <- split(eeg_mdd_wide, eeg_mdd_wide$Session)

# Get feature names
eeg_names <- grep('Coherence|Imag_Coherence|WPLI|WPPC', names(eeg_mdd_wide_by_tp[[1]]), value = TRUE)

# PCA on Session 1
pc_tp1 <- prcomp(eeg_mdd_wide_by_tp[[1]][, eeg_names], scale. = TRUE)

# --- Prepare HC Dataset ---
eeg_hc <- eeg2 %>%
  filter(Population == 'HC') %>%
  select(Subject, Session, Channel_1, Channel_2, FreqBand, Coherence, Imag_Coherence, WPLI, WPPC)

# Only keep Session 1 for HC
eeg_hc_baseline_raw <- eeg_hc %>%
  filter(Session == 1)

# Melt HC into long format 
m_hc <- melt(eeg_hc_baseline_raw, id.vars = c('Subject', 'Session', 'Channel_1', 'Channel_2', 'FreqBand'))

# pivot_wider for HC
eeg_hc_baseline <- m_hc %>%
  pivot_wider(
    names_from = c(Channel_1, Channel_2, FreqBand, variable),
    values_from = value
  )

# Align HC columns to MDD PCA
missing_features <- setdiff(eeg_names, names(eeg_hc_baseline))
for (feature in missing_features) {
  eeg_hc_baseline[[feature]] <- 0  # Fill missing columns with zeros
}
eeg_hc_baseline <- eeg_hc_baseline[, eeg_names]  # Match correct column order

# --- Project data ---
# MDD Session 3
eeg_mdd_session3 <- eeg_mdd_wide_by_tp[["3"]]
pc_mdd_session3_scores <- predict(pc_tp1, newdata = eeg_mdd_session3[, eeg_names])

# HC Baseline
pc_hc_baseline_scores <- predict(pc_tp1, newdata = eeg_hc_baseline)

# --- Combine scores ---
combined_scores <- rbind(
  data.frame(Group = "MDD_Session3", pc_mdd_session3_scores),
  data.frame(Group = "HC_Baseline", pc_hc_baseline_scores)
)

# Reshape for testing
combined_long <- combined_scores %>%
  pivot_longer(
    cols = starts_with("PC"),
    names_to = "Component",
    values_to = "Score"
  )

# --- Statistical tests ---
t_test_results <- combined_long %>%
  group_by(Component) %>%
  summarize(
    t_test_p = t.test(Score ~ Group)$p.value,
    mean_MDD = mean(Score[Group == "MDD_Session3"]),
    mean_HC = mean(Score[Group == "HC_Baseline"])
  ) %>%
  arrange(t_test_p)

# View t-test results
print(t_test_results)

# Visualization 
library(ggplot2)

ggplot(combined_long, aes(x = Group, y = Score, fill = Group)) +
  geom_boxplot() +
  facet_wrap(~Component, scales = "free_y") +
  theme_minimal() +
  labs(title = "Comparison of PCA Scores: MDD Week16 vs HC Baseline",
       x = "", y = "PCA Score") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Print full t-test results table without truncation
print(t_test_results, n = Inf)
