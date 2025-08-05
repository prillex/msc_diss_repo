# Matt Prill
# MSc Diss
# Final Analysis


# Seed
set.seed(1234)

# NMDS and Clustering Analysis ----

# Load packages ----
library(vegan)
library(ggplot2)
library(dplyr)
library(patchwork)
library(cluster)
library(factoextra)

# Load data ----
load("data/merged_data/filtered_matrix.Rdata")
meta <- read.csv("data/merged_data/combined_metadata_filtered.csv")

# Binary matrix
filtered_matrix_bin <- ifelse(is.na(filtered_matrix), 0, 1)

# Compute Jaccard distance and NMDS ----
dist_jaccard <- vegdist(filtered_matrix_bin, method = "jaccard")
set.seed(1234)
nmds <- metaMDS(dist_jaccard, k = 2, trymax = 100)
nmds_data <- as.data.frame(scores(nmds, display = "sites"))
nmds_data$sampleID <- rownames(nmds_data)
nmds_data <- left_join(nmds_data, meta, by = "sampleID")
paste("Stress =", round(nmds$stress, 3))

# NMDS plotting function ----
plot_nmds <- function(data, group, title, legend_title, colours = c("#6E4318", "darkblue")) {
  data$group_var <- as.factor(data[[group]])
  
  hull_data <- bind_rows(lapply(split(data, data$group_var), function(g) {
    if (nrow(g) >= 3) g[chull(g$NMDS1, g$NMDS2), ] else NULL
  }))
  
  ggplot() +
    geom_polygon(data = hull_data, aes(NMDS1, NMDS2, fill = group_var, group = group_var), alpha = 0.3) +
    geom_point(data = data, aes(NMDS1, NMDS2, colour = group_var), size = 3) +
    scale_colour_manual(values = colours) +
    scale_fill_manual(values = colours) +
    theme_classic() +
    labs(title = title, x = "\nNMDS1", y = "NMDS2\n", colour = legend_title, fill = legend_title) +
    theme(
      legend.position = "right",
      legend.justification = c(0, 1),
      text = element_text(size = 18),
      legend.text = element_text(size = 16),
      legend.title = element_text(size = 18)
    )
}

# Generate NMDS plots ----
p1 <- plot_nmds(nmds_data, "Outcome30days", "A)", "Outcome (30 Days)")
p2 <- plot_nmds(nmds_data, "Collection", "B)", "Collection")
p3 <- plot_nmds(nmds_data, "source", "C)", "Source")


# Silhouette analysis for clustering ----
sil_widths <- sapply(2:5, function(k) {
  mean(silhouette(cutree(hclust(dist_jaccard), k = k), dist_jaccard)[, 3])
})

plot(2:5, sil_widths, type = "b", pch = 19,
     xlab = "Number of Clusters", ylab = "Average Silhouette Width",
     main = "Silhouette Analysis (Jaccard)")

optimal_k <- which.max(sil_widths) + 1
cat("Optimal number of clusters based on silhouette width:", optimal_k, "\n")

df_silhouette <- data.frame(
  k = 2:5,
  width = sil_widths
)

# Final Figure
# Plot
silhouette_plot <- ggplot(df_silhouette, aes(x = k, y = width)) +
  geom_line(color = "black", size = 1.2) +
  geom_point(size = 4, colour = "darkblue") +
  scale_x_continuous(breaks = 2:5) +
  labs(
    x = "\nNumber of Clusters",
    y = "Average Silhouette Width\n"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 36, face = "bold", hjust = 0.5),
    axis.title.x = element_text(colour = "black", size = 40),
    axis.title.y = element_text(colour = "black", size = 40),
    axis.text.x = element_text(colour = "black", size = 30),
    axis.text.y = element_text(colour = "black", size = 30),
  )

# Display plot
silhouette_plot


ggsave(silhouette_plot, filename = "code/code/figures/clustering/silhouette_plot.png",
       width = 16, height = 12, dpi = 300, bg = "white")




# Add clusters to data and plot (and meta) ----
nmds_data$cluster <- as.factor(cutree(hclust(dist_jaccard), k = optimal_k))
meta_with_cluster <- left_join(meta, nmds_data[, c("sampleID", "cluster")], by = "sampleID")
meta_with_cluster <- meta_with_cluster[, -c(1, 2)]  # remove redundant columns
#write.csv(meta_with_cluster, "data/merged_data/meta_with_cluster.csv", row.names = FALSE)  # save



p4 <- ggplot(nmds_data, aes(NMDS1, NMDS2, colour = cluster)) +
  geom_point(size = 3) +
  scale_colour_manual(values = c("darkblue", "#6E4318", "purple" )) +
  theme_classic() +
  labs(title = "D)", x = "\nNMDS1", y = "NMDS2\n", colour = "Cluster") +
  theme(
    legend.position = "right",
    legend.justification = c(0, 1),
    text = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 18)
  )



# Combine and save plots ----
(combined_nmds <- (p1 | p2) / (p3 | p4))

#ggsave("code/code/figures/clustering/nmds_combined_plot.png",
#       plot = combined_nmds, width = 14, height = 12, dpi = 300)

# PERMANOVA ----
cat("\nPERMANOVA Results:\n")
print(adonis2(dist_jaccard ~ source, data = meta, permutations = 999))
print(adonis2(dist_jaccard ~ Collection, data = meta, permutations = 999))
print(adonis2(dist_jaccard ~ Outcome30days, data = meta, permutations = 999))

# ----
# ----
# ----
# Feature Selection ----
# Load packages ----
library(tidyverse)
library(caret)

# Data ----
load("data/merged_data/filtered_matrix.Rdata")
meta <- read.csv("data/merged_data/combined_metadata_filtered.csv")

filtered_matrix_bin <- ifelse(is.na(filtered_matrix), 0, 1)


# Correlated Peak Removal ----
# Identify correlated features (above the threshold)
cor_matrix <- cor(filtered_matrix_bin)
high_cor <- findCorrelation(cor_matrix, cutoff = 0.7)
high_cor_names <- colnames(filtered_matrix_bin)[high_cor]
uncorrelated_matrix <- filtered_matrix_bin[, -high_cor]
# save(uncorrelated_matrix, file = "data/merged_data/uncorrelated_matrix.Rdata")

# RFECV ----
# Load data
load("data/merged_data/uncorrelated_matrix.Rdata")  # uncorrelated matrix
load("data/merged_data/correlated_peak_pairs.Rdata")  # correlated pairs for reference
meta <- read.csv("data/merged_data/meta_with_cluster.csv")

# Make sure  matrix is numeric
X <- as.matrix(uncorrelated_matrix) * 1  # Convert logical to numeric 0/1
y <- ifelse(meta$Outcome30days == "Dead", 1L, 0L)  # Binary target variable

# Load reticulate and specify Python version
library(reticulate)
use_python("C:/Users/prill/AppData/Local/Programs/Python/Python313/python.exe", required = TRUE)

# Run Python RFECV within one py_run_string call
py_run_string("
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.feature_selection import RFECV
from sklearn.model_selection import StratifiedKFold

# Access R variables using r.
X = np.array(r.X)
y = np.array(r.y)

model = LogisticRegression(penalty='l2', solver='liblinear', max_iter=1000)
cv = StratifiedKFold(n_splits=10, shuffle=True, random_state=42)

selector = RFECV(estimator=model, step=1, cv=cv, scoring='roc_auc', min_features_to_select=5)
selector.fit(X, y)

# Assign back to R environment
r.selected_feature_mask = selector.support_.tolist()
r.n_features_selected = selector.n_features_
r.ranking_list = selector.ranking_.tolist()
")

# get  selected features by name
feature_names <- colnames(uncorrelated_matrix)
selected_feature_names <- feature_names[selected_feature_mask]

# Print selected feature names and number of features selected
print(selected_feature_names)
cat("Number of features selected:", n_features_selected, "\n")



# write.csv(reduced_matrix, "data/merged_data/nine_peaks.csv")


# ----
# ----
# ----
# Combining Peaks w metadata for modelling ----
# Data ----
nine_peaks <- read.csv("data/merged_data/nine_peaks.csv", row.names = 1, check.names = FALSE)
meta <- read.csv("data/merged_data/meta_with_cluster.csv")


# Rounding numbers in nine_peaks columns ----
# Step 2: Round the m/z values (column names) to nearest integer
rounded_names <- as.character(round(as.numeric(colnames(nine_peaks))))

# Step 3: Assign the new names
colnames(nine_peaks) <- rounded_names

# Combining peaks, meta and response ----
modelling_data <- cbind(meta, nine_peaks)
#write.csv(modelling_data, "data/merged_data/modelling_data.csv", row.names = FALSE)
modelling_data <- read.csv("data/merged_data/modelling_data.csv", check.names = FALSE)

# Peak frequency (grouped by source) ----
# Select source and peaks
modelling_data <- read.csv("data/merged_data/modelling_data.csv", check.names = FALSE)
peak_cols <- c("2304", "2689", "4318", "4449", "6404", "6628", "6943", "7422", "9631")

modelling_data %>%
  select(source, all_of(peak_cols)) %>%
  group_by(source) %>%
  summarise(across(everything(), ~ mean(.x) * 100)) %>%
  print()
# ----
# ----
# ----
# Toxicity ----
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




# ----
# ----
# ----
# Bayes Models (Cross Validation) ----
# Data ----
# Load and prepare data
modelling_data <- read.csv("data/merged_data/modelling_data.csv", check.names = FALSE)
modelling_data$OutcomeBinary <- factor(modelling_data$Outcome30days == "Dead", levels = c(FALSE, TRUE), labels = c("Alive", "Dead"))

# Re-code binary predictor (Gender)
modelling_data <- modelling_data %>%
  mutate(Gender_Male = ifelse(Gender == "Male", 1, 0))  # Male encoded as 1


# Metadata only (no peaks)
meta_df <- modelling_data %>%
  select(OutcomeBinary, Age, Gender_Male) %>%
  na.omit()


# Meta & Peaks (Full Data)
# Rename peaks to valid names
peak_cols <- grep("^\\d+$", colnames(modelling_data), value = TRUE)  
new_peak_names <- paste0("peak_", peak_cols)
colnames(modelling_data)[colnames(modelling_data) %in% peak_cols] <- new_peak_names
peak_cols <- new_peak_names


# Subset needed variables & create df
vars <- c("OutcomeBinary", "Age", "Gender_Male", peak_cols)

bayes_df <- modelling_data %>%
  select(all_of(vars)) %>%
  na.omit()


# Cross Validation Folds ----
folds_meta_df <- createFolds(meta_df$OutcomeBinary, k = 5, returnTrain = TRUE)
folds_bayes_df <- createFolds(bayes_df$OutcomeBinary, k = 5, returnTrain = TRUE)

# Priors ----
# Uninformative Priors 
uninformative_priors <- c(
  prior(normal(0, 10), class = "b"),             # Wide prior for all coefficients
  prior(normal(0, 10), class = "Intercept")      # Wide prior for intercept
)

# Informative Priors 
informative_priors <- c(
  prior(normal(0.03, 0.01), class = "b", coef = "Age"),
  prior(normal(-1.35, 0.5), class = "b", coef = "Gender_Male"),
  prior(normal(-1.51, 0.061), class = "Intercept")
)

# Bayes 1. (uninformative) ----
cv_uninformative_results <- lapply(folds_meta_df, function(train_idx) {
  train_data <- meta_df[train_idx, ]
  test_data <- meta_df[-train_idx, ]
  
  # Fit the model on the training fold
  bayes_uninformative <- brm(
    formula = OutcomeBinary ~ ., 
    data = train_data,
    family = bernoulli(link = "logit"),
    prior = uninformative_priors,
    chains = 3,
    cores = 3,
    iter = 2000,
    warmup = 200,
    seed = 1234
  )
  
  # Predict on the held-out fold
  test_probs_bayes_uninformative <- fitted(bayes_uninformative, newdata = test_data)[, "Estimate"]
  test_preds_bayes_uninformative <- ifelse(test_probs_bayes_uninformative > 0.5, "Dead", "Alive")
  
  list(
    preds_bayes_uninformative = test_preds_bayes_uninformative,
    probs_bayes_uninformative = test_probs_bayes_uninformative,
    truth_bayes_uninformative = test_data$OutcomeBinary
  )
})

# Combine CV results
cv_uninformative_preds <- unlist(lapply(cv_uninformative_results, \(x) x$preds_bayes_uninformative))
cv_uninformative_probs <- unlist(lapply(cv_uninformative_results, \(x) x$probs_bayes_uninformative))
cv_uninformative_truth <- unlist(lapply(cv_uninformative_results, \(x) as.character(x$truth_bayes_uninformative)))

# Evaluate performance
conf_mat <- caret::confusionMatrix(
  factor(cv_uninformative_preds, levels = c("Alive", "Dead")),
  factor(cv_uninformative_truth, levels = c("Alive", "Dead"))
)
print(conf_mat)

# AUC
truth_bin <- ifelse(cv_uninformative_truth == "Dead", 1, 0)
roc_obj <- pROC::roc(truth_bin, cv_uninformative_probs)
auc_value <- auc(roc_obj)
print(auc_value)

# Optional: Plot ROC
plot(roc_obj, main = "ROC Curve - Uninformative Priors")
beep(8)










# Bayes 2. (informative) ----
cv_informative_results <- lapply(folds_meta_df, function(train_idx) {
  train_data <- meta_df[train_idx, ]
  test_data <- meta_df[-train_idx, ]
  
  # Fit the model on the training fold
  bayes_informative <- brm(
    formula = OutcomeBinary ~ ., 
    data = train_data,
    family = bernoulli(link = "logit"),
    prior = informative_priors,
    chains = 3,
    cores = 3,
    iter = 2000,
    warmup = 200,
    seed = 123
  )
  
  # Predict on the held-out fold
  test_probs_bayes_informative <- fitted(bayes_informative, newdata = test_data)[, "Estimate"]
  test_preds_bayes_informative <- ifelse(test_probs_bayes_informative > 0.5, "Dead", "Alive")
  
  list(
    preds_bayes_informative = test_preds_bayes_informative,
    probs_bayes_informative = test_probs_bayes_informative,
    truth_bayes_informative = test_data$OutcomeBinary
  )
})

# Combine CV results
cv_informative_preds <- unlist(lapply(cv_informative_results, \(x) x$preds_bayes_informative))
cv_informative_probs <- unlist(lapply(cv_informative_results, \(x) x$probs_bayes_informative))
cv_informative_truth <- unlist(lapply(cv_informative_results, \(x) as.character(x$truth_bayes_informative)))

# Confusion Matrix
conf_mat_informative <- caret::confusionMatrix(
  factor(cv_informative_preds, levels = c("Alive", "Dead")),
  factor(cv_informative_truth, levels = c("Alive", "Dead"))
)
print(conf_mat_informative)

# AUC
truth_bin_informative <- ifelse(cv_informative_truth == "Dead", 1, 0)
roc_obj_informative <- pROC::roc(truth_bin_informative, cv_informative_probs)
auc_value_informative <- auc(roc_obj_informative)
print(auc_value_informative)

# Optional: Plot ROC
plot(roc_obj_informative, main = "ROC Curve - Informative Priors")
beep(8)

















# Bayes 3. (Peaks + uninformative) ----
# Cross-Validation for Bayes 3: Peaks + Uninformative Priors
cv_uninformative_peaks_results <- lapply(folds_bayes_df, function(train_idx) {
  train_data <- bayes_df[train_idx, ]
  test_data <- bayes_df[-train_idx, ]
  
  # Fit the model on the training fold
  bayes_uninformative_peaks <- brm(
    formula = OutcomeBinary ~ ., 
    data = train_data,
    family = bernoulli(link = "logit"),
    prior = uninformative_priors,
    chains = 3,
    cores = 3,
    iter = 2000,
    warmup = 200,
    seed = 123
  )
  
  # Predict on the held-out fold
  test_probs_bayes_uninformative_peaks <- fitted(bayes_uninformative_peaks, newdata = test_data)[, "Estimate"]
  test_preds_bayes_uninformative_peaks <- ifelse(test_probs_bayes_uninformative_peaks > 0.5, "Dead", "Alive")
  
  list(
    preds_bayes_uninformative_peaks = test_preds_bayes_uninformative_peaks,
    probs_bayes_uninformative_peaks = test_probs_bayes_uninformative_peaks,
    truth_bayes_uninformative_peaks = test_data$OutcomeBinary
  )
})

# Combine CV results
cv_uninformative_peaks_preds <- unlist(lapply(cv_uninformative_peaks_results, \(x) x$preds_bayes_uninformative_peaks))
cv_uninformative_peaks_probs <- unlist(lapply(cv_uninformative_peaks_results, \(x) x$probs_bayes_uninformative_peaks))
cv_uninformative_peaks_truth <- unlist(lapply(cv_uninformative_peaks_results, \(x) as.character(x$truth_bayes_uninformative_peaks)))

# Confusion Matrix
conf_mat_uninformative_peaks <- caret::confusionMatrix(
  factor(cv_uninformative_peaks_preds, levels = c("Alive", "Dead")),
  factor(cv_uninformative_peaks_truth, levels = c("Alive", "Dead"))
)
print(conf_mat_uninformative_peaks)

# AUC
truth_bin_uninformative_peaks <- ifelse(cv_uninformative_peaks_truth == "Dead", 1, 0)
roc_obj_uninformative_peaks <- pROC::roc(truth_bin_uninformative_peaks, cv_uninformative_peaks_probs)
auc_value_uninformative_peaks <- auc(roc_obj_uninformative_peaks)
print(auc_value_uninformative_peaks)

# Optional: ROC plot
plot(roc_obj_uninformative_peaks, main = "ROC Curve - Uninformative Priors (Peaks)")
beep(8)














# Bayes 4. (Peaks + informative Priors) ----
cv_peaks_results <- lapply(folds_bayes_df, function(train_idx) {
  train_data <- bayes_df[train_idx, ]
  test_data <- bayes_df[-train_idx, ]
  
  # Fit model on training fold
  bayes_peaks <- brm(
    formula = OutcomeBinary ~ ., 
    data = train_data,
    family = bernoulli(link = "logit"),
    prior = informative_priors,
    chains = 3,
    cores = 3,
    iter = 2000,
    warmup = 200,
    seed = 123
  )
  
  # Predict on test fold
  test_probs_peaks <- fitted(bayes_peaks, newdata = test_data)[, "Estimate"]
  test_preds_peaks <- ifelse(test_probs_peaks > 0.5, "Dead", "Alive")
  
  list(
    preds_peaks = test_preds_peaks,
    probs_peaks = test_probs_peaks,
    truth_peaks = test_data$OutcomeBinary
  )
})

# Combine CV predictions and truths
cv_peaks_preds <- unlist(lapply(cv_peaks_results, \(x) x$preds_peaks))
cv_peaks_probs <- unlist(lapply(cv_peaks_results, \(x) x$probs_peaks))
cv_peaks_truth <- unlist(lapply(cv_peaks_results, \(x) as.character(x$truth_peaks)))

# Confusion Matrix
conf_mat_peaks <- caret::confusionMatrix(
  factor(cv_peaks_preds, levels = c("Alive", "Dead")),
  factor(cv_peaks_truth, levels = c("Alive", "Dead"))
)
print(conf_mat_peaks)

# AUC
truth_bin_peaks <- ifelse(cv_peaks_truth == "Dead", 1, 0)
roc_obj_peaks <- pROC::roc(truth_bin_peaks, cv_peaks_probs)
auc_value_peaks <- auc(roc_obj_peaks)
print(auc_value_peaks)

# Optional ROC curve
plot(roc_obj_peaks, main = "ROC Curve - Informative Priors (Peaks)")
beep(8)




# ----
# ----
# ----
# ----
# Summary ----
summary(bayes_uninformative)         # Meta: Uninformative
summary(bayes_informative)           # Meta: Informative
summary(bayes_uninformative_peaks)   # Peaks: Uninformative
summary(bayes_peaks)                 # Peaks: Informative

# Linearity of Age on logitscale check ----
library(splines)

bayes_age_spline <- brm(
  OutcomeBinary ~ ns(Age, df = 3) + Gender_Male,
  data = meta_df,
  family = bernoulli(link = "logit"),
  prior = uninformative_priors,
  chains = 3,
  cores = 3
)

loo_compare(loo(bayes_uninformative), loo(bayes_age_spline))



meta_df$OutcomeNum <- ifelse(meta_df$OutcomeBinary == "Dead", 1, 0)

ggplot(meta_df, aes(x = Age, y = OutcomeNum)) +
  geom_smooth(method = "loess", se = FALSE) +
  labs(title = "Observed Probability of Death vs. Age")

# Convergence Checks ----
# Posterior Predictive Checks
pp_check(bayes_uninformative)         # Meta: Uninformative
pp_check(bayes_informative)           # Meta: Informative
pp_check(bayes_uninformative_peaks)   # Peaks: Uninformative
pp_check(bayes_peaks)                 # Peaks: Informative

# Caterpillar Plots
plot(bayes_uninformative)
plot(bayes_informative)
plot(bayes_uninformative_peaks)
plot(bayes_peaks)
# Table ----
# Load required packages
library(dplyr)
library(caret)
library(pROC)
library(ggplot2)

# Function to compute performance metrics from predictions and truth
evaluate_model <- function(preds, probs, truth) {
  cm <- caret::confusionMatrix(factor(preds, levels = c("Alive", "Dead")),
                               factor(truth, levels = c("Alive", "Dead")))
  
  truth_bin <- ifelse(truth == "Dead", 1, 0)
  roc_obj <- pROC::roc(truth_bin, probs)
  auc_val <- pROC::auc(roc_obj)
  
  data.frame(
    Accuracy = cm$overall["Accuracy"],
    F1_Dead = cm$byClass["F1"],
    AUC = as.numeric(auc_val)
  )
}

# Evaluate all models
results_uninformative     <- evaluate_model(cv_uninformative_preds, cv_uninformative_probs, cv_uninformative_truth)
results_informative       <- evaluate_model(cv_informative_preds, cv_informative_probs, cv_informative_truth)
results_uninformative_peaks <- evaluate_model(cv_uninformative_peaks_preds, cv_uninformative_peaks_probs, cv_uninformative_peaks_truth)
results_peaks             <- evaluate_model(cv_peaks_preds, cv_peaks_probs, cv_peaks_truth)

# Combine into one summary table
model_comparison <- bind_rows(
  results_uninformative,
  results_informative,
  results_uninformative_peaks,
  results_peaks
) %>%
  mutate(Model = c("Meta: Uninformative", "Meta: Informative", "Peaks: Uninformative", "Peaks: Informative")) %>%
  select(Model, everything())

print(model_comparison)



# Metrics Plot ----
# Convert to long format for plotting
model_long <- tidyr::pivot_longer(
  model_comparison,
  cols = c("Accuracy", "F1_Dead", "AUC"),
  names_to = "Metric",
  values_to = "Value"
)

# Plot
ggplot(model_long, aes(x = Model, y = Value, fill = Metric)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(aes(label = round(Value, 2)), 
            position = position_dodge(width = 0.7), vjust = -0.5, size = 3.5) +
  labs(title = "Model Performance Comparison", y = "Score", x = NULL) +
  theme_minimal() +
  scale_fill_brewer(palette = "Set2") +
  theme(axis.text.x = element_text(angle = 15, hjust = 1))

# Credibility Interval Plot (Model 4) ----
# Combine credible intervals
# Load necessary libraries
library(brms)
library(ggplot2)
library(dplyr)
library(tidyr)
library(posterior)
library(stringr)

# labels 
bayes_peaks_ci <- bayes_peaks_ci %>%
  mutate(
    term_clean = term %>%
      str_replace("^b_", "") %>%
      str_replace_all("_", " ") %>%
      str_replace("^peak", "Peak") %>%
      str_replace("^Peak", "Peak") %>%
      str_replace("^Gender male$", "Gender male") %>%
      str_replace("^Intercept$", "Intercept")
  )

# Make factor with order by median
bayes_peaks_ci$term_clean <- factor(bayes_peaks_ci$term_clean, levels = bayes_peaks_ci$term_clean[order(bayes_peaks_ci$median)])

cap_width <- 0.15  # whiskers

(CI_int <- ggplot(bayes_peaks_ci, aes(x = term_clean, y = median)) +
  # vertical lines (CI whiskers)
  geom_segment(aes(x = term_clean, xend = term_clean, y = lower, yend = upper),
               size = 1, color = "darkblue") +
  # horizontal caps on lower end
  geom_segment(aes(x = as.numeric(term_clean) - cap_width,
                   xend = as.numeric(term_clean) + cap_width,
                   y = lower, yend = lower),
               color = "darkblue", size = 1) +
  # horizontal caps on upper end
  geom_segment(aes(x = as.numeric(term_clean) - cap_width,
                   xend = as.numeric(term_clean) + cap_width,
                   y = upper, yend = upper),
               color = "darkblue", size = 1) +
  # points for median
  geom_point(size = 3, color = "darkblue") +
  coord_flip() + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.9) +
  labs(
    title = "95% Credible Intervals for Model Coefficients",
    x = "Predictor\n",
    y = "\nLog-Odds Estimate"
  ) +
  theme(
    axis.title.x = element_text(colour = "black", size = 25),
    axis.title.y = element_text(colour = "black", size = 25),
    axis.text.x = element_text(colour = "black", size = 20),
    axis.text.y = element_text(colour = "black", size = 20),
    axis.line = element_line(color = "black", size = 1.5),         # thicker axis lines
    axis.ticks = element_line(color = "black", size = 1.5),        # thicker ticks
    axis.ticks.length = unit(0.25, "cm"),                          # longer ticks
    plot.title = element_text(size = 18, face = "bold", color = "black"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ))



# Save using ggsave with wrap_plots()
ggsave("code/code/figures/outcome/CI_int.png",
       plot = CI_int,
       width = 14, height = 12, dpi = 300)


# ----
# ----
# ----
# ----
# Instant Models (Saved & Ready to Load) ----
# Save models
#saveRDS(bayes_uninformative, "code/code/analysis/models/bayes_uninformative.rds")
#saveRDS(bayes_informative, "code/code/analysis/models/bayes_informative.rds")
#saveRDS(bayes_uninformative_peaks, "code/code/analysis/models/bayes_uninformative_peaks.rds")
#saveRDS(bayes_peaks, "code/code/analysis/models/bayes_peaks.rds")




# Load models
bayes_uninformative <- readRDS("code/code/analysis/models/bayes_uninformative.rds")
bayes_informative <- readRDS("code/code/analysis/models/bayes_informative.rds")
bayes_uninformative_peaks <- readRDS("code/code/analysis/models/bayes_uninformative_peaks.rds")
bayes_peaks <- readRDS("code/code/analysis/models/bayes_peaks.rds")



# Example Of How To Refit ----
# Summary
summary(bayes_peaks)

# Posterior predictive checks
pp_check(bayes_peaks)
plot(bayes_peaks)  # Caterpillar-style posterior plots

# Predictions
probs_peaks <- fitted(bayes_peaks)[, "Estimate"]
preds_peaks <- ifelse(probs_peaks > 0.5, "Dead", "Alive")
truth_peaks <- bayes_df$OutcomeBinary  # Ensure bayes_df is defined

# Confusion Matrix
caret::confusionMatrix(
  factor(preds_peaks, levels = c("Alive", "Dead")),
  factor(truth_peaks, levels = c("Alive", "Dead"))
)

# AUC
library(pROC)
truth_bin_peaks <- ifelse(truth_peaks == "Dead", 1, 0)
roc_obj_peaks <- roc(truth_bin_peaks, probs_peaks)
auc(roc_obj_peaks)




# Alternative helper function  ----
# Function for evaluation (training) ----
evaluate_model <- function(model, data, label) {
  probs <- fitted(model, newdata = data)[, "Estimate"]
  preds <- ifelse(probs > 0.5, "Dead", "Alive")
  truth <- data$OutcomeBinary
  
  cm <- caret::confusionMatrix(
    factor(preds, levels = c("Alive", "Dead")),
    factor(truth, levels = c("Alive", "Dead"))
  )
  
  truth_bin <- ifelse(truth == "Dead", 1, 0)
  roc_obj <- pROC::roc(truth_bin, probs)
  auc_val <- as.numeric(pROC::auc(roc_obj))
  
  # Handle binary and multiclass F1 scoring
  f1_alive <- tryCatch({
    if (is.matrix(cm$byClass)) cm$byClass["F1", "Class: Alive"] else cm$byClass["F1"]
  }, error = function(e) NA)
  
  f1_dead <- tryCatch({
    if (is.matrix(cm$byClass)) cm$byClass["F1", "Class: Dead"] else cm$byClass["F1"]
  }, error = function(e) NA)
  
  tibble(
    Model = label,
    Accuracy = cm$overall["Accuracy"],
    F1_Alive = f1_alive,
    F1_Dead = f1_dead,
    AUC = auc_val
  )
}


# Evaluate all models (training)
results <- bind_rows(
  evaluate_model(bayes_uninformative, meta_df, "Meta - Uninformative"),
  evaluate_model(bayes_informative, meta_df, "Meta - Informative"),
  evaluate_model(bayes_uninformative_peaks, bayes_df, "Peaks - Uninformative"),
  evaluate_model(bayes_peaks, bayes_df, "Peaks - Informative")
)

print(results)

# Function for evaluation (testing) ----
evaluate_cv_results <- function(preds, probs, truth, label) {
  cm <- caret::confusionMatrix(
    factor(preds, levels = c("Alive", "Dead")),
    factor(truth, levels = c("Alive", "Dead"))
  )
  
  truth_bin <- ifelse(truth == "Dead", 1, 0)
  roc_obj <- pROC::roc(truth_bin, probs)
  auc_val <- as.numeric(pROC::auc(roc_obj))
  
  tibble(
    Model = label,
    Accuracy = cm$overall["Accuracy"],
    F1_Alive = cm$byClass["F1"],
    F1_Dead = cm$byClass["F1"],  # For binary, they are the same
    AUC = auc_val
  )
}

# Model 1
evaluate_cv_results(
  preds = cv_uninformative_preds,
  probs = cv_uninformative_probs,
  truth = cv_uninformative_truth,
  label = "Meta - Uninformative (CV)"
)

# Model 2
evaluate_cv_results(
  preds = cv_informative_preds,
  probs = cv_informative_probs,
  truth = cv_informative_truth,
  label = "Meta - Uninformative (CV)"
)


# Model 3
evaluate_cv_results(
  preds = cv_uninformative_peaks_preds,
  probs = cv_uninformative_peaks_probs,
  truth = cv_uninformative_peaks_truth,
  label = "Meta - Uninformative (CV)"
)

# Model 4
evaluate_cv_results(
  preds = cv_peaks_preds,
  probs = cv_peaks_probs,
  truth = cv_peaks_truth,
  label = "Meta - Uninformative (CV)"
)


# Credibility Intervals ----
# Function to print 95% credibility intervals
print_cred_intervals <- function(model, model_name) {
  cat("\n", "=== ", model_name, " ===", "\n")
  summary_data <- posterior_summary(model, probs = c(0.025, 0.975))
  coef_summary <- summary_data[grep("^b_", rownames(summary_data)), ]  # just the coefficients
  print(round(coef_summary, 3))
}

print_cred_intervals(bayes_peaks, "Peaks Model (Complete)")
# ----
# ----
# ----




