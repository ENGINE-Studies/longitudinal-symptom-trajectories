# Load required libraries
library(tidyverse)
library(reshape2)
library(dplyr)
library(tidyr)

# --- Read EEG data ---
eeg2 <- read.csv('//Users/wadimvodovozov/Desktop/ENGINE/analysis/master_connectivity_SPARSE_17-Nov-2024.csv')

# --- Prepare MDD Dataset ---
eeg_mdd <- eeg2 %>%
  filter(Population == 'MDD') %>%
  select(Subject, Session, Channel_1, Channel_2, FreqBand, Coherence, Imag_Coherence, WPLI, WPPC)

# Only include MDD participants with complete sessions 1,2,3
mdd_subjects_complete <- eeg_mdd %>%
  group_by(Subject) %>%
  filter(all(c(1,2,3) %in% Session)) %>%
  pull(Subject) %>%
  unique()

# Filter MDD baseline data (Session 1) for only complete subjects
eeg_mdd_baseline <- eeg_mdd %>%
  filter(Session == 1, Subject %in% mdd_subjects_complete)

# Melt MDD baseline data
m_mdd_baseline <- melt(eeg_mdd_baseline, id.vars = c('Subject', 'Session', 'Channel_1', 'Channel_2', 'FreqBand'))

# Reshape to wide format
eeg_mdd_baseline_wide <- m_mdd_baseline %>%
  pivot_wider(
    names_from = c(Channel_1, Channel_2, FreqBand, variable),
    values_from = value
  )

# --- Prepare HC Dataset ---
eeg_hc <- eeg2 %>%
  filter(Population == 'HC') %>%
  select(Subject, Session, Channel_1, Channel_2, FreqBand, Coherence, Imag_Coherence, WPLI, WPPC)

eeg_hc_baseline_raw <- eeg_hc %>%
  filter(Session == 1)

# Melt and reshape HC
m_hc <- melt(eeg_hc_baseline_raw, id.vars = c('Subject', 'Session', 'Channel_1', 'Channel_2', 'FreqBand'))

eeg_hc_baseline_wide <- m_hc %>%
  pivot_wider(
    names_from = c(Channel_1, Channel_2, FreqBand, variable),
    values_from = value
  )

# --- Perform PCA on MDD baseline ---
eeg_names <- grep('Coherence|Imag_Coherence|WPLI|WPPC', names(eeg_mdd_baseline_wide), value = TRUE)

pc_mdd_baseline <- prcomp(eeg_mdd_baseline_wide[, eeg_names], scale. = TRUE)

# ✅ No slicing anymore -> retain all PCs
n_subjects_mdd <- length(unique(eeg_mdd_baseline$Subject))

# Align HC columns
missing_features <- setdiff(eeg_names, names(eeg_hc_baseline_wide))
for (feature in missing_features) {
  eeg_hc_baseline_wide[[feature]] <- 0
}
eeg_hc_baseline_wide <- eeg_hc_baseline_wide[, eeg_names]

# --- Project Data ---
# Project HC baseline into full MDD PCA space (no slicing!)
pc_hc_projected <- predict(pc_mdd_baseline, newdata = eeg_hc_baseline_wide)

# Scores for MDD baseline
pc_mdd_scores <- pc_mdd_baseline$x

# --- Combine Scores ---
combined_baseline_scores <- rbind(
  data.frame(Group = "MDD_Baseline", pc_mdd_scores),
  data.frame(Group = "HC_Baseline", pc_hc_projected)
)

# --- Reshape long ---
combined_baseline_long <- combined_baseline_scores %>%
  pivot_longer(
    cols = starts_with("PC"),
    names_to = "Component",
    values_to = "Score"
  )

# --- Statistical tests ---
t_test_results_baseline <- combined_baseline_long %>%
  group_by(Component) %>%
  summarize(
    t_test_p = t.test(Score ~ Group)$p.value,
    mean_MDD = mean(Score[Group == "MDD_Baseline"]),
    mean_HC = mean(Score[Group == "HC_Baseline"])
  ) %>%
  arrange(t_test_p)

# --- Final formatted printout ---
t_test_results_baseline %>%
  mutate(Row = row_number()) %>%
  select(Row, Component, t_test_p, mean_MDD, mean_HC) %>%
  print(n = Inf)

# --- Visualization ---
combined_baseline_long$Component <- factor(
  combined_baseline_long$Component,
  levels = paste0("PC", 1:(n_subjects_mdd))  # not n_subjects_mdd-1 !!
)

ggplot(combined_baseline_long, aes(x = Group, y = Score, fill = Group)) +
  geom_boxplot() +
  facet_wrap(~Component, scales = "free_y") +
  theme_minimal() +
  labs(title = "Comparison of PCA Scores: MDD Baseline vs HC Baseline",
       x = "", y = "PCA Score") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
