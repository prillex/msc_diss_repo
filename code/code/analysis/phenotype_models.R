# Matt Prill
# MSc Applied Data Science Dissertation
# S. aureus Phenotype Modelling


# Libraries ----
library(randomForest)
library(caret)         # for model evaluation
library(dplyr)         # for data manipulation

# Seed ----
set.seed(1234)

# Data ----
load("data/merged_data/uncorrelated_matrix.Rdata")  # uncorrelated matrix
meta <- read.csv("data/merged_data/combined_metadata_filtered.csv")  # Meta

# Subset to just CC22 data
meta <- meta[match(rownames(uncorrelated_matrix), meta$sampleID), ]  #  Matches them up again

# Create a filter for samples that do NOT start with 'A' followed by digits (i.e. cc22)
keep <- !grepl("^A\\d+$", meta$sampleID)

# Subset both the matrix and meta
peak_matrix_cc22 <- uncorrelated_matrix[keep, ]
meta_cc22 <- meta[keep, ]




# Recursive Feature Selection ----
library(reticulate)
use_python("C:/Users/prill/AppData/Local/Programs/Python/Python313/python.exe", required = TRUE)

# Prepare X and y for toxicity binary
X <- as.matrix(peak_matrix_cc22) * 1  # numeric matrix 0/1
y <- ifelse(meta_cc22$tox_cat == "high", 1L, 0L)  # binary target (high toxicity=1)

# Push X and y to Python environment explicitly
py$X <- X
py$y <- y

# Run RFECV in Python via reticulate
py_run_string("
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.feature_selection import RFECV
from sklearn.model_selection import StratifiedKFold

X = np.array(r.X)
y = np.array(r.y)

model = LogisticRegression(penalty='l2', solver='liblinear', max_iter=1000)  # logit = default
cv = StratifiedKFold(n_splits=10, shuffle=True, random_state=42)

selector = RFECV(estimator=model, step=1, cv=cv, scoring='roc_auc', min_features_to_select=5)
selector.fit(X, y)

r.selected_feature_mask = selector.support_.tolist()
r.n_features_selected = selector.n_features_
r.ranking_list = selector.ranking_.tolist()
")


# Get selected feature names
feature_names <- colnames(peak_matrix_cc22)
selected_feature_names <- feature_names[selected_feature_mask]

cat("Number of features selected:", n_features_selected, "\n")
print(selected_feature_names)


fivethree_peaks <- peak_matrix_cc22[, selected_feature_mask]

# Check dimensions
dim(fivethree_peaks)

# Round the m/z values (column names) to nearest integer
rounded_names <- as.character(round(as.numeric(colnames(fivethree_peaks))))

# Step 3: Assign the new names
colnames(fivethree_peaks) <- rounded_names

# Save the selected matrix to an Rdata file
#save(fivethree_peaks, file = "data/merged_data/selected_peaks_toxicity.Rdata")



# Random forest (toxicity) ----
# Data ----
load("data/merged_data/selected_peaks_toxicity.Rdata")

# Toxicity (Continuous) 
df_tox <- data.frame(toxicity = meta_cc22$tox_cat,  fivethree_peaks)

# Convert tox_cat to factor with levels 0 and 1
df_tox$toxicity <- factor(ifelse(df_tox$toxicity == "high", 1, 0), levels = c(0, 1))

table(df_tox$toxicity)


# The Model ----
rf_tox <- randomForest(toxicity ~ ., data = df_tox, importance = TRUE, ntree = 500)

# Summary
print(rf_tox)

# Variable importance plot
varImpPlot(rf_tox)

# Extract importance and convert to data frame
importance_df <- as.data.frame(importance(rf_tox))
importance_df$Feature <- gsub("^X", "", rownames(importance_df))  # remove leading 'X'

# For accuracy:
importance_df <- importance_df %>%
  arrange(desc(MeanDecreaseAccuracy)) %>%
  slice(1:10)  # Top 10 features


custom_colors <- c("3026" = "#6E4318",   # blue
                   "3006" = "darkblue",   # red
                   "2977" = "purple")   # green


# Assign colors: default black, override if in custom_colors
importance_df$bar_color <- ifelse(importance_df$Feature %in% names(custom_colors),
                                  custom_colors[importance_df$Feature],
                                  "grey")

# Plot with ggplot
(rf_importance <- ggplot(importance_df, aes(x = reorder(Feature, MeanDecreaseAccuracy),
                                           y = MeanDecreaseAccuracy,
                                           fill = bar_color)) +
  geom_col() +
  scale_fill_identity() +  # use actual hex values in bar_color
  coord_flip() +
  scale_y_continuous(expand = c(0, 0)) +
  labs(
    title = "Top 10 Important Features for Toxicity Classification",
    x = "Peak\n",
    y = "\nMean Decrease in Accuracy"
  ) +
  theme_classic(base_size = 16) +
  theme(
    axis.text.y = element_text(size = 20, face = "bold", color = "black"),
    axis.text.x = element_text(size = 20, face = "bold", color = "black"),
    axis.title.y = element_text(size = 24, face = "bold", color = "black"),
    axis.title.x = element_text(size = 24, face = "bold", color = "black"),
    plot.title = element_text(size = 18, face = "bold", color = "black"),
    plot.margin = margin(t = 10, r = 15, b = 10, l = 10)
  ))




# Save using ggsave with wrap_plots()
#ggsave("code/code/figures/toxicity/tox_importance.png",
#       plot = rf_importance,
#       width = 14, height = 12, dpi = 300)


# Checking correlated peaks
# Load correlated peak pairs
load("data/merged_data/correlated_peak_pairs.Rdata")

top_features <- as.numeric(importance_df$Feature[1:10])


# Flatten the list into a data frame
# Round the features in the list to, say, 0 decimal places
cor_df <- do.call(rbind, lapply(correlated_pairs, function(x) {
  data.frame(
    feature1 = round(as.numeric(x$feature1)),
    feature2 = round(as.numeric(x$feature2)),
    correlation = x$correlation
  )
}))


correlated_with_top <- cor_df %>%
  filter(feature1 %in% top_features | feature2 %in% top_features)




# Cooccurence  ----

library(ggplot2)
library(UpSetR)  # For Venn-like intersection visualisation

# Prepare data
peaks_upset <- df_tox[, c("X3026", "X3006", "X2977")]
colnames(peaks_upset) <- c("Peak 3026", "Peak 3006", "Peak 2977")

# Assign custom colors for each peak
custom_colors <- c("Peak 3026" = "#6E4318",   # blue
                   "Peak 3006" = "darkblue",   # red
                   "Peak 2977" = "purple")   # green

# Reorder dataframe to match colors
peaks_upset <- peaks_upset[, names(custom_colors)]

# Plot
(co_occurence <- upset(peaks_upset,
                     sets = names(custom_colors),
                     nsets = 3,
                     keep.order = TRUE,
                     order.by = "freq",
                     main.bar.color = "darkgrey",  # keep main bars neutral or edit as needed
                     matrix.color = "darkgrey",     # still only one color for matrix dots
                     sets.bar.color = unname(custom_colors),
                     mainbar.y.label = "Intersection Count",
                     sets.x.label = "Peak Occurrence Count",
                     text.scale = c(3, 2.5, 2.5, 2.5, 3, 3)))

# Save using ggsave with wrap_plots()
#ggsave("code/code/figures/clustering/co_occurence.png",
#       plot = wrap_plots(co_occurence),
#       width = 14, height = 12, dpi = 300)





# Plotting toxicities (mean) ----

library(dplyr)
library(ggplot2)

# Start with your processed df_combined (with combo labels)
df_combined <- df_tox %>%
  mutate(
    toxicity = meta_cc22$toxicity,
    combo_raw = paste0(X3026, X3006, X2977)
  )

# Label combinations
combo_labels <- c(
  "100" = "3026",
  "010" = "3006",
  "001" = "2977",
  "110" = "3026 × 3006",
  "101" = "3026 × 2977",
  "011" = "3006 × 2977",
  "111" = "All"
)

# Remove "000" (no peaks)
df_combined <- df_combined %>%
  filter(combo_raw != "000") %>%
  mutate(combo = combo_labels[combo_raw])

# Summarise mean and SD of toxicity
tox_summary <- df_combined %>%
  group_by(combo) %>%
  summarise(
    mean_tox = mean(toxicity, na.rm = TRUE),
    sd_tox = sd(toxicity, na.rm = TRUE),
    count = n()
  ) %>%
  arrange(desc(mean_tox))

# Plot barplot with error bars
(barplot_tox <- ggplot(tox_summary, aes(x = reorder(combo, -mean_tox), y = mean_tox, fill = combo)) +
    geom_col(width = 0.6, color = "black") +
    geom_errorbar(aes(ymin = mean_tox - sd_tox, ymax = mean_tox + sd_tox),
                  width = 0.2, size = 0.6) +
    labs(
      x = "\nPeak Combination",
      y = "Mean Toxicity (%)\n",
      title = "Mean Toxicity Across Peak Combinations"
    ) +
    scale_fill_manual(values = rep("#6E4318", length(unique(tox_summary$combo)))) +
    theme_classic(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 14, color = "black", face = "bold"),
      axis.text.y = element_text(size = 14, color = "black"),
      axis.title = element_text(size = 18, color = "black", face = "bold"),
      plot.title = element_text(size = 18, face = "bold", color = "black"),
      legend.position = "none"
    ))











# Boxplots ----
df_combined <- df_tox %>%
  mutate(toxicity = meta_cc22$toxicity)  # Add continuous toxicity

# Create a combination label
df_combined$combination <- with(df_combined, paste0(X3026, "_", X3006, "_", X2977))

# Group by that combo and summarise
tox_by_combo <- df_combined %>%
  group_by(combination) %>%
  summarise(
    count = n(),
    mean_toxicity = mean(toxicity, na.rm = TRUE),
    sd_toxicity = sd(toxicity, na.rm = TRUE)
  ) %>%
  arrange(desc(mean_toxicity))

# View result
print(tox_by_combo)



# Combine peaks and toxicity
df_combined <- df_tox %>%
  mutate(
    toxicity = meta_cc22$toxicity,
    combo_raw = paste0(X3026, X3006, X2977)
  )

# Map combinations to descriptive names
combo_labels <- c(
  "000" = "None",
  "100" = "3026",
  "010" = "3006",
  "001" = "2977",
  "110" = "3026 × 3006",
  "101" = "3026 × 2977",
  "011" = "3006 × 2977",
  "111" = "All"
)

df_combined <- df_combined %>%
  mutate(combo = combo_labels[combo_raw])  # Apply labels

# Optional: Filter out rare combos if needed (e.g., n < 2)
combo_counts <- df_combined %>%
  count(combo) %>%
  filter(n >= 1)  # change to 2 if you want to exclude very rare ones

filtered_df <- df_combined %>%
  filter(combo %in% combo_counts$combo)

# Plot
(tox_combo_plot <- ggplot() +
  geom_violin(data = filtered_df,
              aes(x = combo, y = toxicity, fill = combo),
              alpha = 0.3, trim = TRUE, size = 0.2, scale = "width") +
  stat_boxplot(data = filtered_df,
               aes(x = combo, y = toxicity),
               geom = "errorbar", width = 0.4, size = 0.5, coef = Inf) +
  geom_boxplot(data = filtered_df,
               aes(x = combo, y = toxicity, fill = combo),
               width = 0.5, outlier.shape = NA) +
  geom_jitter(data = df_combined,
              aes(x = combo, y = toxicity),
              color = "black", width = 0.2, height = 0, size = 1, alpha = 0.9) +
  labs(
    y = "Toxicity (%)\n",
    x = "\nPeak Combination",
    title = "Toxicity Across Peak Combinations"
  ) +
  theme_classic(base_size = 14) +
  scale_fill_manual(values = rep("#6E4318", length(unique(filtered_df$combo)))) +
  theme(
    axis.title.x = element_text(colour = "black", size = 20),
    axis.title.y = element_text(colour = "black", size = 20),
    axis.text.x = element_text(colour = "black", size = 14, angle = 45, hjust = 1),
    axis.text.y = element_text(colour = "black", size = 14),
    axis.line = element_line(color = "black"),
    plot.title = element_text(size = 18, face = "bold", color = "black"),
    plot.margin = unit(c(20, 15, 15, 15), "pt"),
    panel.background = element_rect(fill = "transparent"),
    plot.background = element_rect(fill = "transparent", color = NA),
  ))




# Without 'none' ----
library(dplyr)
library(ggplot2)

# Add continuous toxicity to df_tox
df_combined <- df_tox %>%
  mutate(toxicity = meta_cc22$toxicity)

# Create a combination code
df_combined$combination <- with(df_combined, paste0(X3026, "_", X3006, "_", X2977))

# Exclude '000' (None) before summary
tox_by_combo <- df_combined %>%
  filter(!(combination == "0_0_0")) %>%
  group_by(combination) %>%
  summarise(
    count = n(),
    mean_toxicity = mean(toxicity, na.rm = TRUE),
    sd_toxicity = sd(toxicity, na.rm = TRUE)
  ) %>%
  arrange(desc(mean_toxicity))

# View result
print(tox_by_combo)


# Prepare combo labels
df_combined <- df_combined %>%
  mutate(
    combo_raw = paste0(X3026, X3006, X2977)
  )

# Combo label map
combo_labels <- c(
  "100" = "3026",
  "010" = "3006",
  "001" = "2977",
  "110" = "3026 × 3006",
  "101" = "3026 × 2977",
  "011" = "3006 × 2977",
  "111" = "All"
)

# Remove '000' and apply labels
df_combined <- df_combined %>%
  filter(combo_raw != "000") %>%
  mutate(combo = combo_labels[combo_raw])

# Optional: Filter out rare combos if needed
combo_counts <- df_combined %>%
  count(combo) %>%
  filter(n >= 1)

filtered_df <- df_combined %>%
  filter(combo %in% combo_counts$combo)

# Plot
(tox_combo_plot <- ggplot() +
    geom_violin(data = filtered_df,
                aes(x = combo, y = toxicity, fill = combo),
                alpha = 0.3, trim = TRUE, size = 0.2, scale = "width") +
    stat_boxplot(data = filtered_df,
                 aes(x = combo, y = toxicity),
                 geom = "errorbar", width = 0.4, size = 0.5, coef = Inf) +
    geom_boxplot(data = filtered_df,
                 aes(x = combo, y = toxicity, fill = combo),
                 width = 0.5, outlier.shape = NA) +
    geom_jitter(data = df_combined,
                aes(x = combo, y = toxicity),
                color = "black", width = 0.2, height = 0, size = 1, alpha = 0.9) +
    labs(
      y = "Toxicity (%)\n",
      x = "\nPeak Combination",
      title = "Toxicity Across Peak Combinations"
    ) +
    theme_classic(base_size = 14) +
    scale_fill_manual(values = rep("#6E4318", length(unique(filtered_df$combo)))) +
    theme(
      axis.title.x = element_text(colour = "black", size = 20),
      axis.title.y = element_text(colour = "black", size = 20),
      axis.text.x = element_text(colour = "black", size = 14, angle = 45, hjust = 1),
      axis.text.y = element_text(colour = "black", size = 14),
      axis.line = element_line(color = "black"),
      plot.title = element_text(size = 18, face = "bold", color = "black"),
      plot.margin = unit(c(20, 15, 15, 15), "pt")
    ))

# One-way ANOVA -----
anova_result <- aov(toxicity ~ combo, data = filtered_df)
summary(anova_result)

levels(filtered_df$combo)



# Zoomed ----
# Ensure combo is a factor with correct levels (same order as in your plot)
levels_order <- unique(filtered_df$combo)
filtered_df$combo <- factor(filtered_df$combo, levels = levels_order)

# Filter combos with at least 5 points for violin and boxplot
filtered_for_plots <- filtered_df %>%
  group_by(combo) %>%
  filter(n() >= 5)

# Recalculate box_stats with consistent factor levels
box_stats <- filtered_for_plots %>%
  group_by(combo) %>%
  summarise(
    ymin = quantile(toxicity, 0.25) - 1.5 * IQR(toxicity),
    ymin_adj = pmax(ymin, 80)  # Ensure whisker extends at least to 80
  ) %>%
  mutate(x = as.numeric(combo))  # Get numeric x position


library(RColorBrewer)

# Define color palette for points only
palette_colors <- brewer.pal(n = max(3, length(levels_order)), name = "Set2")

tox_combo_plot <- ggplot() +
  # Violin plots in grey
  geom_violin(data = filtered_for_plots,
              aes(x = combo, y = toxicity),
              fill = "grey80", alpha = 0.3, trim = TRUE, size = 0.2, scale = "width") +
  
  # Boxplot error bars in grey
  stat_boxplot(data = filtered_for_plots,
               aes(x = combo, y = toxicity),
               geom = "errorbar", width = 0.4, size = 0.5, coef = Inf, color = "grey40") +
  
  # Boxplots in grey
  geom_boxplot(data = filtered_for_plots,
               aes(x = combo, y = toxicity),
               width = 0.5, outlier.shape = NA, fill = "lightgrey", color = "grey40") +
  
  # Extended whisker segment for "All"
  geom_segment(data = filter(box_stats, combo == "All"),
               aes(x = x, xend = x, y = ymin_adj, yend = 80),
               color = "black", size = 0.5) +
  
  # Colored jitter points
  geom_jitter(data = filtered_df,
              aes(x = combo, y = toxicity, color = combo),
              width = 0.1, height = 0, size = 3, alpha = 0.9) +
  
  coord_cartesian(ylim = c(80, 100)) +
  labs(
    y = "Toxicity\n",
    x = "\nPeak Combination",
    title = "Toxicity Across Peak Combinations"
  ) +
  theme_classic(base_size = 14) +
  scale_color_manual(values = palette_colors) +  # for points
  guides(fill = "none") +  # suppress fill legend if any
  theme(
    axis.title.x = element_text(colour = "black", size = 24),
    axis.title.y = element_text(colour = "black", size = 24),
    axis.text.x = element_text(colour = "black", size = 20, angle = 45, hjust = 1),
    axis.text.y = element_text(colour = "black", size = 20),
    axis.line = element_line(color = "black"),
    plot.title = element_text(size = 5, face = "bold", color = "black"),
    plot.margin = unit(c(20, 15, 15, 15), "pt")
  )

print(tox_combo_plot)

# Save using ggsave with wrap_plots()
#ggsave("code/code/figures/toxicity/tox_combo_plot.png",
#       plot = tox_combo_plot,
#       width = 14, height = 12, dpi = 300)

