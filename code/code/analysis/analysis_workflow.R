# Matt Prill
# MSc Diss
# Machine learning workflow

# Seed ----
set.seed(1234)

# Libraries ----
library(tidyverse)
library(umap)
library(uwot)

# UMAP ----
# Peak Matrix
load("data/merged_data/filtered_matrix.Rdata")
meta <- read.csv("data/merged_data/combined_metadata_filtered.csv")  # metadata

# run UMAP
filtered_matrix_bin <- ifelse(is.na(filtered_matrix), 0, 1)  # remove NA, and make binary

umap_res <- umap(filtered_matrix_bin)  # this will be a matrix (uwot)

# Clustering ----
# k-means clustering (try k = 3 based on visual clusters)
kmeans_res <- kmeans(umap_res, centers = 3)
cluster_labels <- as.factor(kmeans_res$cluster)

# Plot UMAP with cluster labels
plot(umap_res, col = cluster_labels, pch = 19,
     main = "UMAP with K-means Clusters", xlab = "UMAP1", ylab = "UMAP2")
legend("topright", legend = levels(cluster_labels),
       col = 1:length(levels(cluster_labels)), pch = 19)

# Plot UMAP coloured by Outcome (30-day mortality)
plot(umap_res, col = as.factor(meta$Outcome30days), pch = 19,
     main = "UMAP of Filtered Peaks (Mortality)", xlab = "UMAP1", ylab = "UMAP2")
legend("topright", legend = levels(as.factor(meta$Outcome30days)),
       col = 1:length(unique(meta$Outcome30days)), pch = 19)

# Plot UMAP coloured by MRSA/MSSA
plot(umap_res, col = as.factor(meta$Collection), pch = 19,
     main = "UMAP of Filtered Peaks (MRSA/MSSA)", xlab = "UMAP1", ylab = "UMAP2")
legend("topright", legend = levels(as.factor(meta$Collection)),
       col = 1:length(unique(meta$Collection)), pch = 19)

# Plot UMAP coloured by source
plot(umap_res, col = as.factor(meta$source), pch = 19,
     main = "UMAP of Filtered Peaks (Source)", xlab = "UMAP1", ylab = "UMAP2")
legend("topright", legend = levels(as.factor(meta$source)),
       col = 1:length(unique(meta$source)), pch = 19)




# Seeing which peaks drive data structure










# UMAP but not binary ----

# Replace NA with 0 but keep numeric values (non-binary)
filtered_matrix_clean <- ifelse(is.na(filtered_matrix), 0, filtered_matrix)

# Run UMAP
umap_res <- umap(filtered_matrix_clean)  # returns matrix from uwot

# UMAP Clustering ----
# K-means clustering (try k = 3)
kmeans_res <- kmeans(umap_res, centers = 3)
cluster_labels <- as.factor(kmeans_res$cluster)

# Plot UMAP with cluster labels
plot(umap_res, col = cluster_labels, pch = 19,
     main = "UMAP with K-means Clusters (Non-binary)", xlab = "UMAP1", ylab = "UMAP2")
legend("topright", legend = levels(cluster_labels),
       col = 1:length(levels(cluster_labels)), pch = 19)

# Plot UMAP coloured by Outcome (30-day mortality)
plot(umap_res, col = as.factor(meta$Outcome30days), pch = 19,
     main = "UMAP of Filtered Peaks (Mortality)", xlab = "UMAP1", ylab = "UMAP2")
legend("topright", legend = levels(as.factor(meta$Outcome30days)),
       col = 1:length(unique(meta$Outcome30days)), pch = 19)

# Plot UMAP coloured by MRSA/MSSA
plot(umap_res, col = as.factor(meta$Collection), pch = 19,
     main = "UMAP of Filtered Peaks (MRSA/MSSA)", xlab = "UMAP1", ylab = "UMAP2")
legend("topright", legend = levels(as.factor(meta$Collection)),
       col = 1:length(unique(meta$Collection)), pch = 19)

# Plot UMAP coloured by source
plot(umap_res, col = as.factor(meta$source), pch = 19,
     main = "UMAP of Filtered Peaks (Source)", xlab = "UMAP1", ylab = "UMAP2")
legend("topright", legend = levels(as.factor(meta$source)),
       col = 1:length(unique(meta$source)), pch = 19)
















# NMDS ----
library(vegan)

# Peak Matrix
load("data/merged_data/filtered_matrix.Rdata")
meta <- read.csv("data/merged_data/combined_metadata_filtered.csv")  # metadata


# run UMAP
filtered_matrix_bin <- ifelse(is.na(filtered_matrix), 0, 1)  # remove na, and make binary

# Step 2: Compute Jaccard distance & NMDS (2D)
dist_jaccard <- vegdist(filtered_matrix_bin, method = "jaccard")
nmds_res <- metaMDS(dist_jaccard, k = 2, trymax = 100)

# Step 3: Extract NMDS coordinates
nmds_scores <- as.data.frame(scores(nmds_res, display = "sites"))
nmds_scores$sampleID <- rownames(nmds_scores)

# Step 4: Join with metadata
nmds_data <- left_join(nmds_scores, meta, by = "sampleID")

# Optional: create a combined grouping variable
nmds_data <- nmds_data %>%
  mutate(group_var = paste0(Outcome90days, "-", source))

# Step 5: Get convex hulls
hull_data <- data.frame()
for (grp in unique(nmds_data$group_var)) {
  temp <- nmds_data[nmds_data$group_var == grp, ]
  if (nrow(temp) >= 3) {  # Need ≥3 points for a hull
    hull_pts <- temp[chull(temp$NMDS1, temp$NMDS2), ]
    hull_data <- rbind(hull_data, hull_pts)
  }
}
































# NMDS and Clustering Analysis ----

# --- Load packages ---
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

# Add clusters to data and plot ----
nmds_data$cluster <- as.factor(cutree(hclust(dist_jaccard), k = optimal_k))

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












# Feature Selection (removing co-linear peaks) ----
library(tidyverse)
library(caret)


load("data/merged_data/filtered_matrix.Rdata")
meta <- read.csv("data/merged_data/combined_metadata_filtered.csv")

filtered_matrix_bin <- ifelse(is.na(filtered_matrix), 0, 1)


# Identify correlated features (above the threshold)
cor_matrix <- cor(filtered_matrix_bin)
high_cor <- findCorrelation(cor_matrix, cutoff = 0.7)
high_cor_names <- colnames(filtered_matrix_bin)[high_cor_indices]
uncorrelated_matrix <- filtered_matrix_bin[, -high_cor_indices]
# save(uncorrelated_matrix, file = "data/merged_data/uncorrelated_matrix.Rdata")


# Number of peaks before filtering
n_before <- ncol(filtered_matrix_bin)

# Number of peaks after filtering
n_after <- ncol(uncorrelated_matrix)

# Print summary
cat("Peaks before filtering:", n_before, "\n")
cat("Peaks after filtering:", n_after, "\n")
cat("Number of peaks removed:", n_before - n_after, "\n")


# Document the correlated pairs
get_correlated_pairs <- function(cor_matrix, cutoff = 0.7) {
  cor_pairs <- which(abs(cor_matrix) > cutoff & lower.tri(cor_matrix), arr.ind = TRUE)
  pair_list <- apply(cor_pairs, 1, function(idx) {
    list(
      feature1 = colnames(cor_matrix)[idx[1]],
      feature2 = colnames(cor_matrix)[idx[2]],
      correlation = cor_matrix[idx[1], idx[2]]
    )
  })
  return(pair_list)
}

correlated_pairs <- get_correlated_pairs(cor_matrix, cutoff = 0.7)
print(correlated_pairs[1:5])
#save(correlated_pairs, file = "data/merged_data/correlated_peak_pairs.Rdata")




# Feature Importance Dimensionality Reduction ----
load("data/merged_data/uncorrelated_matrix.Rdata")  # uncorrelated matrix
load("data/merged_data/correlated_peak_pairs.Rdata")  # correlated pairs for reference
meta <- read.csv("data/merged_data/meta_with_cluster.csv")

# library
library(reticulate)

# Ensure matrix is numeric
X <- as.matrix(uncorrelated_matrix) * 1  # force 0/1
y <- ifelse(meta$Outcome30days == "Dead", 1L, 0L)  # 1 = Dead

# Make available to Python
reticulate::use_virtualenv("r-reticulate")
py$X <- X
py$y <- y


py_run_string("
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.feature_selection import RFECV
from sklearn.model_selection import StratifiedKFold

# Define model
model = LogisticRegression(penalty='l2', solver='liblinear', max_iter=1000)

# Define stratified CV, each fold had the same split at the data 79:21 alive:dead
cv = StratifiedKFold(n_splits=10, shuffle=True, random_state=42)

# RFECV
# starts with all features and removes the least important one at each step
# At each step it uses 10-fold CV to test performance
# it score each subset based on AUROC (ability to distinguish Alive vs Dead)
selector = RFECV(estimator=model, step=1, cv=cv, scoring='roc_auc', min_features_to_select=5)
selector.fit(X, y)

# Outputs
optimal_num_features = selector.n_features_
selected_feature_mask = selector.support_
ranking = selector.ranking_
")

ranking <- py$ranking
# Optimal number of features
optimal_features <- py$optimal_num_features
cat("Optimal number of features:", optimal_features, "\n")

# Boolean mask of selected features
selected_mask <- py$selected_feature_mask
selected_indices <- which(selected_mask)

# Get column names of selected features
selected_feature_names <- colnames(uncorrelated_matrix)[selected_indices]
print(selected_feature_names)

reduced_matrix <- uncorrelated_matrix[, selected_feature_names]
nine_peaks <- reduced_matrix 
# write.csv(reduced_matrix, "data/merged_data/nine_peaks.csv")

# Combine with the meta
nine_peaks <- read.csv("data/merged_data/nine_peaks.csv", row.names = 1, check.names = FALSE)
meta <- read.csv("data/merged_data/meta_with_cluster.csv")

modelling_data <- cbind(meta, nine_peaks)
write.csv(modelling_data, "data/merged_data/modelling_data.csv", row.names = FALSE)












# Failed random forest/ bayesian optimisation ----
library(caret)
library(pROC)
library(randomForest)
library(ParBayesianOptimization)
library(dplyr)

set.seed(1234)

# Create outer stratified folds (returns test indices)
folds <- createFolds(model_df$OutcomeBinary, k = 10, returnTrain = FALSE)

# Store results
auc_vals <- c()
f1_dead_vals <- c()
f1_alive_vals <- c()
bal_acc_vals <- c()

for (i in seq_along(folds)) {
  cat(sprintf("Processing Fold %d...\n", i))
  
  test_indices <- folds[[i]]
  train_indices <- setdiff(seq_len(nrow(model_df)), test_indices)
  
  train_data <- model_df[train_indices, ]
  test_data <- model_df[test_indices, ]
  
  # Check for missing class in test fold; skip if so
  if (any(table(test_data$OutcomeBinary) == 0)) {
    warning(sprintf("Fold %d skipped due to missing class in test set", i))
    next
  }
  
  train_data_copy <- train_data  # freeze this fold’s training data
  
  # Optimisation function
  opt_function <- function(mtry, maxnodes, nodesize) {
    mtry <- floor(mtry)
    maxnodes <- floor(maxnodes)
    nodesize <- floor(nodesize)
    
    internal_folds <- createFolds(train_data_copy$OutcomeBinary, k = 5, returnTrain = FALSE)
    sens_vals <- c()
    
    for (j in seq_along(internal_folds)) {
      val_idx <- internal_folds[[j]]
      train_idx_inner <- setdiff(seq_len(nrow(train_data_copy)), val_idx)
      
      train_inner <- train_data_copy[train_idx_inner, ]
      val_inner <- train_data_copy[val_idx, ]
      
      model_inner <- randomForest(
        OutcomeBinary ~ ., 
        data = train_inner,
        mtry = mtry,
        maxnodes = maxnodes,
        nodesize = nodesize,
        ntree = 500
      )
      
      pred_inner <- predict(model_inner, newdata = val_inner, type = "response")
      obs_inner <- factor(val_inner$OutcomeBinary, levels = c("Alive", "Dead"))
      pred_inner <- factor(pred_inner, levels = c("Alive", "Dead"))
      
      cm_inner <- confusionMatrix(pred_inner, obs_inner)
      sens_vals[j] <- cm_inner$byClass["Sensitivity"]
    }
    
    return(list(Score = mean(sens_vals, na.rm = TRUE)))
  }
  
  # Define hyperparameter bounds
  bounds <- list(
    mtry = c(2L, ncol(train_data) - 1L),
    maxnodes = c(5L, 30L),
    nodesize = c(1L, 10L)
  )
  
  # Run Bayesian Optimisation
  opt_result <- bayesOpt(
    FUN = opt_function,
    bounds = bounds,
    initPoints = 5,
    iters.n = 25,
    acq = "ei",
    gsPoints = 100,  # helps model fit better
    verbose = 0
  )
  
  # Extract best hyperparameters
  best_params <- getBestPars(opt_result)
  
  # Train final RF model on entire train_data with best hyperparameters
  final_model <- randomForest(
    OutcomeBinary ~ .,
    data = train_data,
    mtry = floor(best_params$mtry),
    maxnodes = floor(best_params$maxnodes),
    nodesize = floor(best_params$nodesize),
    ntree = 500
  )
  
  # Predict on outer test fold
  preds <- predict(final_model, newdata = test_data, type = "response")
  probs <- as.data.frame(predict(final_model, newdata = test_data, type = "prob"))
  
  obs <- factor(test_data$OutcomeBinary, levels = c("Alive", "Dead"))
  pred <- factor(preds, levels = c("Alive", "Dead"))
  
  # AUROC
  roc_obj <- roc(obs, probs$Dead)
  auc_vals[i] <- as.numeric(auc(roc_obj))
  
  # F1 scores
  f1_dead_vals[i] <- F_meas(pred, obs, relevant = "Dead")
  f1_alive_vals[i] <- F_meas(pred, obs, relevant = "Alive")
  
  # Balanced accuracy
  cm <- confusionMatrix(pred, obs)
  bal_acc_vals[i] <- cm$byClass["Balanced Accuracy"]
}

# Report mean ± SD of metrics
cat("Random Forest Performance (Nested CV with Bayesian Optimisation, 10 folds):\n")
cat(sprintf("AUROC: %.2f ± %.2f\n", mean(auc_vals, na.rm = TRUE), sd(auc_vals, na.rm = TRUE)))
cat(sprintf("F1-score (Dead): %.2f ± %.2f\n", mean(f1_dead_vals, na.rm = TRUE), sd(f1_dead_vals, na.rm = TRUE)))
cat(sprintf("F1-score (Alive): %.2f ± %.2f\n", mean(f1_alive_vals, na.rm = TRUE), sd(f1_alive_vals, na.rm = TRUE)))
cat(sprintf("Balanced Accuracy: %.2f ± %.2f\n", mean(bal_acc_vals, na.rm = TRUE), sd(bal_acc_vals, na.rm = TRUE)))














# XGBoost ----
library(caret)
library(pROC)
library(dplyr)
library(xgboost)  # Make sure this is installed

folds <- createFolds(model_df$OutcomeBinary, k = 10, returnTrain = TRUE)

# Store results
auc_vals <- c()
f1_dead_vals <- c()
f1_alive_vals <- c()
bal_acc_vals <- c()

for (i in seq_along(folds)) {
  cat(sprintf("Processing Fold %d...\n", i))
  
  train_indices <- folds[[i]]
  train_data <- model_df[train_indices, ]
  test_data <- model_df[-train_indices, ]
  
  if (any(table(test_data$OutcomeBinary) == 0)) {
    warning(sprintf("Fold %d skipped: missing class in test set", i))
    next
  }
  
  # Inner CV control for tuning
  inner_ctrl <- trainControl(
    method = "cv",
    number = 5,
    classProbs = TRUE,
    summaryFunction = twoClassSummary,
    verboseIter = FALSE
  )
  
  # Tuning grid for xgboost
  tune_grid <- expand.grid(
    nrounds = c(100, 200),           # number of boosting rounds
    max_depth = c(3, 6, 9),          # tree depth
    eta = c(0.01, 0.1, 0.3),         # learning rate
    gamma = 0,                       # min loss reduction to make split
    colsample_bytree = 1,            # subsample ratio of columns per tree
    min_child_weight = 1,            # min sum of instance weight needed in child
    subsample = 1                   # subsample ratio of training instances
  )
  
  # Train model with inner CV and grid search
  xgb_model <- train(
    OutcomeBinary ~ .,
    data = train_data,
    method = "xgbTree",
    metric = "ROC",
    trControl = inner_ctrl,
    tuneGrid = tune_grid
  )
  
  # Predict on outer test set
  preds <- predict(xgb_model, newdata = test_data)
  probs <- predict(xgb_model, newdata = test_data, type = "prob")
  
  obs <- factor(test_data$OutcomeBinary, levels = c("Alive", "Dead"))
  pred <- factor(preds, levels = c("Alive", "Dead"))
  
  # AUROC
  roc_obj <- roc(obs, probs$Dead)
  auc_vals[i] <- as.numeric(auc(roc_obj))
  
  # F1-scores
  f1_dead_vals[i] <- F_meas(pred, obs, relevant = "Dead")
  f1_alive_vals[i] <- F_meas(pred, obs, relevant = "Alive")
  
  # Balanced Accuracy
  cm <- confusionMatrix(pred, obs)
  bal_acc_vals[i] <- cm$byClass["Balanced Accuracy"]
}

cat("XGBoost Performance (Nested CV, 10 folds):\n")
cat(sprintf("AUROC: %.2f ± %.2f\n", mean(auc_vals, na.rm = TRUE), sd(auc_vals, na.rm = TRUE)))
cat(sprintf("F1-score (Dead): %.2f ± %.2f\n", mean(f1_dead_vals, na.rm = TRUE), sd(f1_dead_vals, na.rm = TRUE)))
cat(sprintf("F1-score (Alive): %.2f ± %.2f\n", mean(f1_alive_vals, na.rm = TRUE), sd(f1_alive_vals, na.rm = TRUE)))
cat(sprintf("Balanced Accuracy: %.2f ± %.2f\n", mean(bal_acc_vals, na.rm = TRUE), sd(bal_acc_vals, na.rm = TRUE)))





# CatBoost ----
library(reticulate)
library(caret)
library(pROC)
library(dplyr)



# Create numeric outcome for Python
modelling_data$OutcomeBinaryNum <- ifelse(modelling_data$OutcomeBinary == "Dead", 1L, 0L)

# Identify peak columns (those that are purely numeric names)
peak_cols <- grep("^\\d+$", colnames(modelling_data), value = TRUE)

# Build input data: peaks + numeric outcome
peaks_input <- modelling_data[, c(peak_cols, "OutcomeBinaryNum")]

# --- 4. Pass data to Python ---
df_py <- r_to_py(peaks_input)
py$df <- df_py


py_run_string("
import pandas as pd
from catboost import CatBoostClassifier, Pool
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score, f1_score, balanced_accuracy_score

# Use DataFrame from R
X = df.drop(columns=['OutcomeBinaryNum'])
y = df['OutcomeBinaryNum']

skf = StratifiedKFold(n_splits=10, shuffle=True, random_state=42)

auc_vals = []
f1_dead_vals = []
f1_alive_vals = []
bal_acc_vals = []

for train_idx, test_idx in skf.split(X, y):
    X_train, X_test = X.iloc[train_idx], X.iloc[test_idx]
    y_train, y_test = y.iloc[train_idx], y.iloc[test_idx]

    train_pool = Pool(X_train, y_train)
    test_pool = Pool(X_test, y_test)

    model = CatBoostClassifier(
        iterations=200,
        learning_rate=0.1,
        depth=6,
        loss_function='Logloss',
        eval_metric='AUC',
        verbose=False,
        random_seed=42
    )

    model.fit(train_pool)

    probs = model.predict_proba(test_pool)[:, 1]
    preds = (probs > 0.5).astype(int)

    auc_vals.append(roc_auc_score(y_test, probs))
    f1_dead_vals.append(f1_score(y_test, preds, pos_label=1))
    f1_alive_vals.append(f1_score(y_test, preds, pos_label=0))
    bal_acc_vals.append(balanced_accuracy_score(y_test, preds))
")


cat("CatBoost Performance (10-fold CV):\n")
cat(sprintf("AUROC: %.2f ± %.2f\n", mean(py$auc_vals), sd(py$auc_vals)))
cat(sprintf("F1-score (Dead): %.2f ± %.2f\n", mean(py$f1_dead_vals), sd(py$f1_dead_vals)))
cat(sprintf("F1-score (Alive): %.2f ± %.2f\n", mean(py$f1_alive_vals), sd(py$f1_alive_vals)))
cat(sprintf("Balanced Accuracy: %.2f ± %.2f\n", mean(py$bal_acc_vals), sd(py$bal_acc_vals)))




# SVM ----
library(caret)
library(pROC)
library(dplyr)
library(e1071)

# find peak columns
peak_cols <- grep("^\\d+$", colnames(modelling_data), value = TRUE)

# create dataset
svm_data <- modelling_data[, c("OutcomeBinary", peak_cols)]
svm_data$OutcomeBinary <- factor(svm_data$OutcomeBinary, levels = c("Alive", "Dead"))

# create 10 folds
folds <- createFolds(svm_data$OutcomeBinary, k = 10, returnTrain = TRUE)

# store results
auc_vals <- c()
f1_dead_vals <- c()
f1_alive_vals <- c()
bal_acc_vals <- c()

# loop over folds
for (i in seq_along(folds)) {
  cat(sprintf("processing fold %d...\n", i))
  
  train_idx <- folds[[i]]
  train_data <- svm_data[train_idx, ]
  test_data <- svm_data[-train_idx, ]
  
  if (any(table(test_data$OutcomeBinary) == 0)) {
    warning(sprintf("fold %d skipped: missing class in test set", i))
    next
  }
  
  ctrl <- trainControl(
    method = "cv",
    number = 5,
    classProbs = TRUE,
    summaryFunction = twoClassSummary
  )
  
  svm_grid <- expand.grid(
    C = 2^(-2:2),
    sigma = 2^(-5:-1)
  )
  
  svm_model <- train(
    OutcomeBinary ~ .,
    data = train_data,
    method = "svmRadial",
    metric = "ROC",
    trControl = ctrl,
    tuneGrid = svm_grid,
    preProcess = c("center", "scale")
  )
  
  probs <- predict(svm_model, newdata = test_data, type = "prob")
  preds <- predict(svm_model, newdata = test_data)
  
  obs <- factor(test_data$OutcomeBinary, levels = c("Alive", "Dead"))
  pred <- factor(preds, levels = c("Alive", "Dead"))
  
  roc_obj <- roc(obs, probs$Dead)
  auc_vals[i] <- as.numeric(auc(roc_obj))
  f1_dead_vals[i] <- F_meas(pred, obs, relevant = "Dead")
  f1_alive_vals[i] <- F_meas(pred, obs, relevant = "Alive")
  bal_acc_vals[i] <- confusionMatrix(pred, obs)$byClass["Balanced Accuracy"]
}

# show results
cat("svm radial performance (10-fold cv):\n")
cat(sprintf("auroc: %.2f ± %.2f\n", mean(auc_vals, na.rm = TRUE), sd(auc_vals, na.rm = TRUE)))
cat(sprintf("f1-score (dead): %.2f ± %.2f\n", mean(f1_dead_vals, na.rm = TRUE), sd(f1_dead_vals, na.rm = TRUE)))
cat(sprintf("f1-score (alive): %.2f ± %.2f\n", mean(f1_alive_vals, na.rm = TRUE), sd(f1_alive_vals, na.rm = TRUE)))
cat(sprintf("balanced accuracy: %.2f ± %.2f\n", mean(bal_acc_vals, na.rm = TRUE), sd(bal_acc_vals, na.rm = TRUE)))




# Bart ----
# Only run if not already installed
library(bartMachine)  # Required me to install JAVA
library(caret)
library(pROC)
library(dplyr)

# for reproducibility
set_bart_machine_num_cores(1)  # you can increase this for speed if your system allows

# subset to OutcomeBinary + peaks only
peak_cols <- grep("^\\d+$", colnames(modelling_data), value = TRUE)
bart_df <- modelling_data[, c("OutcomeBinary", peak_cols)]
bart_df$OutcomeBinary <- factor(bart_df$OutcomeBinary, levels = c("Alive", "Dead"))

# set up 10-fold CV
folds <- createFolds(bart_df$OutcomeBinary, k = 10, returnTrain = TRUE)

auc_vals <- c()
f1_dead_vals <- c()
f1_alive_vals <- c()
bal_acc_vals <- c()

for (i in seq_along(folds)) {
  cat(sprintf("processing fold %d...\n", i))
  
  train_idx <- folds[[i]]
  train_data <- bart_df[train_idx, ]
  test_data <- bart_df[-train_idx, ]
  
  if (any(table(test_data$OutcomeBinary) == 0)) {
    warning(sprintf("fold %d skipped due to missing class in test set", i))
    next
  }
  
  # bartMachine expects numeric outcome (0/1)
  train_data$Y <- ifelse(train_data$OutcomeBinary == "Dead", 1, 0)
  test_data$Y <- ifelse(test_data$OutcomeBinary == "Dead", 1, 0)
  
  train_data$OutcomeBinary <- NULL
  test_data$OutcomeBinary <- NULL
  
  model <- bartMachine::bartMachine(
    X = train_data[, !(names(train_data) %in% "Y")],
    y = train_data$Y,
    num_trees = 50,
    num_burn_in = 200,
    num_iterations_after_burn_in = 1000,
    verbose = FALSE
  )
  
  preds_prob <- predict(model, test_data[, !(names(test_data) %in% "Y")])
  preds_class <- ifelse(preds_prob > 0.5, "Dead", "Alive")
  
  obs <- factor(ifelse(test_data$Y == 1, "Dead", "Alive"), levels = c("Alive", "Dead"))
  pred <- factor(preds_class, levels = c("Alive", "Dead"))
  
  roc_obj <- roc(obs, preds_prob)
  auc_vals[i] <- auc(roc_obj)
  
  f1_dead_vals[i] <- F_meas(pred, obs, relevant = "Dead")
  f1_alive_vals[i] <- F_meas(pred, obs, relevant = "Alive")
  
  cm <- confusionMatrix(pred, obs)
  bal_acc_vals[i] <- cm$byClass["Balanced Accuracy"]
}

cat("bart performance (10-fold cv):\n")
cat(sprintf("auroc: %.2f ± %.2f\n", mean(auc_vals, na.rm = TRUE), sd(auc_vals, na.rm = TRUE)))
cat(sprintf("f1-score (dead): %.2f ± %.2f\n", mean(f1_dead_vals, na.rm = TRUE), sd(f1_dead_vals, na.rm = TRUE)))
cat(sprintf("f1-score (alive): %.2f ± %.2f\n", mean(f1_alive_vals, na.rm = TRUE), sd(f1_alive_vals, na.rm = TRUE)))
cat(sprintf("balanced accuracy: %.2f ± %.2f\n", mean(bal_acc_vals, na.rm = TRUE), sd(bal_acc_vals, na.rm = TRUE)))



# Logistic regression with metadata as predictors ----
# Libraries
library(caret)
library(pROC)
library(dplyr)

# Merge additional predictors
model_df_extended <- model_df
model_df_extended$tox_cat <- modelling_data$tox_cat
model_df_extended$CCI_class <- modelling_data$CCI_class

# Filter to complete cases only (where tox_cat and CCI_class are not missing)
model_df_extended <- model_df_extended %>%
  filter(!is.na(tox_cat) & !is.na(CCI_class))

# Cross-validation control
ctrl <- trainControl( 
  method = "cv",
  number = 10,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = "final"
)

# Train logistic regression model
logit_model_extended <- train(
  OutcomeBinary ~ .,
  data = model_df_extended,
  method = "glm",
  family = "binomial",
  metric = "ROC",
  trControl = ctrl
)

# Evaluate model
preds <- logit_model_extended$pred
preds$obs <- factor(preds$obs, levels = c("Alive", "Dead"))
preds$pred <- factor(preds$pred, levels = c("Alive", "Dead"))

# AUROC
roc_obj <- pROC::roc(preds$obs, preds$Dead)
mean_auc <- pROC::auc(roc_obj)
sd_auc <- sd(logit_model_extended$resample$ROC)

# F1-scores (with NA-safe summarisation)
f1_dead_vals <- preds %>%
  group_by(Resample) %>%
  summarise(f1 = caret::F_meas(pred, obs, relevant = "Dead")) %>%
  pull(f1) %>%
  na.omit()

f1_alive_vals <- preds %>%
  group_by(Resample) %>%
  summarise(f1 = caret::F_meas(pred, obs, relevant = "Alive")) %>%
  pull(f1) %>%
  na.omit()

mean_f1_dead <- mean(f1_dead_vals)
sd_f1_dead <- sd(f1_dead_vals)
mean_f1_alive <- mean(f1_alive_vals)
sd_f1_alive <- sd(f1_alive_vals)

# Balanced Accuracy
bal_acc_vals <- preds %>%
  group_by(Resample) %>%
  summarise(bal_acc = caret::confusionMatrix(pred, obs)$byClass["Balanced Accuracy"]) %>%
  pull(bal_acc)

mean_bal_acc <- mean(bal_acc_vals, na.rm = TRUE)
sd_bal_acc <- sd(bal_acc_vals, na.rm = TRUE)

# Results
cat("Logistic Regression (with tox_cat and CCI_class) Performance (10-fold CV):\n")
cat(sprintf("AUROC: %.2f ± %.2f\n", mean_auc, sd_auc))
cat(sprintf("F1-score (Dead): %.2f ± %.2f\n", mean_f1_dead, sd_f1_dead))
cat(sprintf("F1-score (Alive): %.2f ± %.2f\n", mean_f1_alive, sd_f1_alive))
cat(sprintf("Balanced Accuracy: %.2f ± %.2f\n", mean_bal_acc, sd_bal_acc))










# Bayesian Regression (no peaks) ----

library(brms)
library(dplyr)
library(caret)

# Load and prepare data
modelling_data <- read.csv("data/merged_data/modelling_data.csv", check.names = FALSE)
modelling_data$OutcomeBinary <- factor(modelling_data$Outcome30days == "Dead", levels = c(FALSE, TRUE), labels = c("Alive", "Dead"))

# Recode binary predictors
modelling_data <- modelling_data %>%
  mutate(
    Gender_Male = ifelse(Gender == "Male", 1, 0),
    Collection_MRSA = ifelse(Collection == "MRSA", 1, 0)
  )

# Subset metadata only (no peaks)
meta_df <- modelling_data %>%
  select(OutcomeBinary, Age, Gender_Male, Collection_MRSA, cluster) %>%
  mutate(cluster = factor(cluster)) %>%
  na.omit()

# Priors (same as original model)
priors_meta_only <- c(
  prior(normal(0.03, 0.01), class = "b", coef = "Age"),
  prior(normal(-1.35, 0.5), class = "b", coef = "Gender_Male"),
  prior(normal(0.3, 0.3), class = "b", coef = "Collection_MRSA"),
  prior(normal(1.51, 0.061), class = "Intercept")
)

# Fit model
meta_only_fit <- brm(
  formula = OutcomeBinary ~ Age + Gender_Male + Collection_MRSA + cluster,
  data = meta_df,
  family = bernoulli(link = "logit"),
  prior = priors_meta_only,
  chains = 3,
  cores = 3,
  iter = 2000,
  warmup = 200,
  seed = 123
)

# Model diagnostics
summary(meta_only_fit)
pp_check(meta_only_fit)

# Predictions
meta_probs <- fitted(meta_only_fit)[, "Estimate"]
meta_preds <- ifelse(meta_probs > 0.5, "Dead", "Alive")
meta_truth <- meta_df$OutcomeBinary

# Confusion matrix
caret::confusionMatrix(
  factor(meta_preds, levels = c("Alive", "Dead")),
  factor(meta_truth, levels = c("Alive", "Dead"))
)


# Bayesian Regression (with peaks) ----
library(brms)
library(dplyr)

# Load and prepare data
modelling_data <- read.csv("data/merged_data/modelling_data.csv", check.names = FALSE)
modelling_data$OutcomeBinary <- factor(modelling_data$Outcome30days == "Dead", levels = c(FALSE, TRUE), labels = c("Alive", "Dead"))

# Rename peaks to valid names
peak_cols <- grep("^\\d+$", colnames(modelling_data), value = TRUE)
new_peak_names <- paste0("peak_", peak_cols)
colnames(modelling_data)[colnames(modelling_data) %in% peak_cols] <- new_peak_names
peak_cols <- new_peak_names

# Recode binary predictors manually
modelling_data <- modelling_data %>%
  mutate(
    Gender_Male = ifelse(Gender == "Male", 1, 0),  # 1 = Female
    Collection_MRSA = ifelse(Collection == "MRSA", 1, 0)  # 1 = MRSA
  )

# Subset needed variables
vars <- c("OutcomeBinary", "Age", "Gender_Male", "Collection_MRSA", "cluster", peak_cols)
bayes_df <- modelling_data %>%
  select(all_of(vars)) %>%
  mutate(cluster = factor(cluster)) %>%
  na.omit()

# Check class balance
table(bayes_df$OutcomeBinary)

# Define priors for numeric predictors
priors <- c(
  prior(normal(0.03, 0.01), class = "b", coef = "Age"),
  prior(normal(-1.35, 0.5), class = "b", coef = "Gender_Male"),
  prior(normal(0.3, 0.3), class = "b", coef = "Collection_MRSA"),
  prior(normal(-1.51, 0.061), class = "Intercept")
)

# Fit the model
fit <- brm(
  formula = OutcomeBinary ~ Age + Gender_Male + Collection_MRSA + cluster + ., 
  data = bayes_df,
  family = bernoulli(link = "logit"),
  prior = priors,
  chains = 3,
  cores = 3,
  iter = 2000,
  warmup = 200,
  seed = 123
)


# Evaluate
summary(fit)
pp_check(fit)

# Predict
probs <- fitted(fit)[, "Estimate"]
preds <- ifelse(probs > 0.5, "Dead", "Alive")
truth <- bayes_df$OutcomeBinary

caret::confusionMatrix(factor(preds, levels = c("Alive", "Dead")),
                       factor(truth, levels = c("Alive", "Dead")))

# Alive F1 - 0.91
# Dead F1 - 0.55




# Comparing bayesian models
library(yardstick)
library(dplyr)
library(pROC)

# ----- Prepare data 

# Full model data
df_full <- data.frame(
  truth = factor(truth, levels = c("Alive", "Dead")),
  prediction = factor(preds, levels = c("Alive", "Dead")),
  prob = probs
)

# Metadata-only model data
df_meta <- data.frame(
  truth = factor(meta_truth, levels = c("Alive", "Dead")),
  prediction = factor(meta_preds, levels = c("Alive", "Dead")),
  prob = meta_probs
)

# ----- Calculate metrics 

# Function to get metrics for a given dataframe
get_metrics <- function(df) {
  acc <- accuracy_vec(df$truth, df$prediction)
  f1_alive <- f_meas_vec(df$truth, df$prediction, event_level = "first")  # Alive
  f1_dead  <- f_meas_vec(df$truth, df$prediction, event_level = "second") # Dead
  auc <- roc(as.numeric(df$truth) == 2, df$prob)$auc
  tibble(
    Accuracy = round(acc, 3),
    F1_Alive = round(f1_alive, 3),
    F1_Dead = round(f1_dead, 3),
    AUROC = round(auc, 3)
  )
}

# Get metrics
full_metrics <- get_metrics(df_full)
meta_metrics <- get_metrics(df_meta)

# ----- Combine and print 
comparison_table <- bind_rows(
  Full_Model = full_metrics,
  Metadata_Only = meta_metrics,
  .id = "Model"
)

print(comparison_table)



# ----
# ----
# ----
# Comparing Modelling Frameworks ----
# Data ----
modelling_data <- read.csv("data/merged_data/modelling_data.csv", check.names = FALSE)

# Create OutcomeBinary 
modelling_data$OutcomeBinary <- ifelse(modelling_data$Outcome30days == "Dead", 1, 0)

# Convert to factor for caret
modelling_data$OutcomeBinary <- factor(modelling_data$OutcomeBinary, levels = c(0, 1), labels = c("Alive", "Dead"))

# Identify peak columns
peak_cols <- grep("^\\d+$", colnames(modelling_data), value = TRUE)

# Create modelling dataset
model_df <- modelling_data[, c("OutcomeBinary", peak_cols)]


# Logistic Regression ----

# libraries
library(caret)
library(pROC)

# Cross validation
ctrl <- trainControl( 
  method = "cv",
  number = 10,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,  # for AUROC
  savePredictions = "final"
)


# Train model
logit_model <- train(
  OutcomeBinary ~ .,
  data = model_df,
  method = "glm",
  family = "binomial",
  metric = "ROC",
  trControl = ctrl
)


# Evaluate 
# Extract predictions
preds <- logit_model$pred
preds$obs <- factor(preds$obs, levels = c("Alive", "Dead"))
preds$pred <- factor(preds$pred, levels = c("Alive", "Dead"))

# AUROC
roc_obj <- pROC::roc(preds$obs, preds$Dead)
auc_val <- pROC::auc(roc_obj)
sd_auc <- sd(logit_model$resample$ROC)
mean_auc <- mean(logit_model$resample$ROC)


# F1-scores for both classes
f1_dead_vals <- preds %>%
  group_by(Resample) %>%
  summarise(f1 = caret::F_meas(pred, obs, relevant = "Dead")) %>%
  pull(f1)

f1_alive_vals <- preds %>%
  group_by(Resample) %>%
  summarise(f1 = caret::F_meas(pred, obs, relevant = "Alive")) %>%
  pull(f1)

mean_f1_dead <- mean(f1_dead_vals)
sd_f1_dead <- sd(f1_dead_vals)

mean_f1_alive <- mean(f1_alive_vals)
sd_f1_alive <- sd(f1_alive_vals)

# Balanced Accuracy
bal_acc_vals <- preds %>%
  group_by(Resample) %>%
  summarise(bal_acc = caret::confusionMatrix(pred, obs)$byClass["Balanced Accuracy"]) %>%
  pull(bal_acc)

mean_bal_acc <- mean(bal_acc_vals)
sd_bal_acc <- sd(bal_acc_vals)

# Results
cat("Logistic Regression Performance (10-fold CV):\n")
cat(sprintf("AUROC: %.2f ± %.2f\n", mean_auc, sd_auc))
cat(sprintf("F1-score (Dead): %.2f ± %.2f\n", mean_f1_dead, sd_f1_dead))
cat(sprintf("F1-score (Alive): %.2f ± %.2f\n", mean_f1_alive, sd_f1_alive))
cat(sprintf("Balanced Accuracy: %.2f ± %.2f\n", mean_bal_acc, sd_bal_acc))


# Random Forest ----
library(caret)
library(pROC)
library(dplyr)


folds <- createFolds(model_df$OutcomeBinary, k = 10, returnTrain = TRUE)

# Store results
auc_vals <- c()
f1_dead_vals <- c()
f1_alive_vals <- c()
bal_acc_vals <- c()

for (i in seq_along(folds)) {
  cat(sprintf("Processing Fold %d...\n", i))
  
  train_indices <- folds[[i]]
  train_data <- model_df[train_indices, ]
  test_data <- model_df[-train_indices, ]
  
  if (any(table(test_data$OutcomeBinary) == 0)) {
    warning(sprintf("Fold %d skipped: missing class in test set", i))
    next
  }
  
  # Inner CV control for tuning
  inner_ctrl <- trainControl(
    method = "cv",
    number = 5,
    classProbs = TRUE,
    summaryFunction = twoClassSummary,
    verboseIter = FALSE
  )
  
  # Tuning grid for ranger
  tune_grid <- expand.grid(
    mtry = c(2, 4, 6),            # tuning no. variables at each split
    min.node.size = c(1, 5, 10),  # tuning min node size
    splitrule = "gini"  # classification
  )
  
  # Train model with inner CV and grid search
  rf_model <- train(
    OutcomeBinary ~ .,
    data = train_data,
    method = "ranger",
    metric = "ROC",
    trControl = inner_ctrl,
    tuneGrid = tune_grid,
    num.trees = 500
  )
  
  # Predict on outer test set
  preds <- predict(rf_model, newdata = test_data)
  probs <- predict(rf_model, newdata = test_data, type = "prob")
  
  obs <- factor(test_data$OutcomeBinary, levels = c("Alive", "Dead"))
  pred <- factor(preds, levels = c("Alive", "Dead"))
  
  # AUROC
  roc_obj <- roc(obs, probs$Dead)
  auc_vals[i] <- as.numeric(auc(roc_obj))
  
  # F1-scores
  f1_dead_vals[i] <- F_meas(pred, obs, relevant = "Dead")
  f1_alive_vals[i] <- F_meas(pred, obs, relevant = "Alive")
  
  # Balanced Accuracy
  cm <- confusionMatrix(pred, obs)
  bal_acc_vals[i] <- cm$byClass["Balanced Accuracy"]
}

cat("Ranger Random Forest Performance (Nested CV, 10 folds):\n")
cat(sprintf("AUROC: %.2f ± %.2f\n", mean(auc_vals, na.rm = TRUE), sd(auc_vals, na.rm = TRUE)))
cat(sprintf("F1-score (Dead): %.2f ± %.2f\n", mean(f1_dead_vals, na.rm = TRUE), sd(f1_dead_vals, na.rm = TRUE)))
cat(sprintf("F1-score (Alive): %.2f ± %.2f\n", mean(f1_alive_vals, na.rm = TRUE), sd(f1_alive_vals, na.rm = TRUE)))
cat(sprintf("Balanced Accuracy: %.2f ± %.2f\n", mean(bal_acc_vals, na.rm = TRUE), sd(bal_acc_vals, na.rm = TRUE)))



# XGBoost ----
library(caret)
library(pROC)
library(dplyr)
library(xgboost)  # Make sure this is installed

folds <- createFolds(model_df$OutcomeBinary, k = 10, returnTrain = TRUE)

# Store results
auc_vals <- c()
f1_dead_vals <- c()
f1_alive_vals <- c()
bal_acc_vals <- c()

for (i in seq_along(folds)) {
  cat(sprintf("Processing Fold %d...\n", i))
  
  train_indices <- folds[[i]]
  train_data <- model_df[train_indices, ]
  test_data <- model_df[-train_indices, ]
  
  if (any(table(test_data$OutcomeBinary) == 0)) {
    warning(sprintf("Fold %d skipped: missing class in test set", i))
    next
  }
  
  # Inner CV control for tuning
  inner_ctrl <- trainControl(
    method = "cv",
    number = 5,
    classProbs = TRUE,
    summaryFunction = twoClassSummary,
    verboseIter = FALSE
  )
  
  # Tuning grid for xgboost
  tune_grid <- expand.grid(
    nrounds = c(100, 200),           # number of boosting rounds
    max_depth = c(3, 6, 9),          # tree depth
    eta = c(0.01, 0.1, 0.3),         # learning rate
    gamma = 0,                       # min loss reduction to make split
    colsample_bytree = 1,            # subsample ratio of columns per tree
    min_child_weight = 1,            # min sum of instance weight needed in child
    subsample = 1                   # subsample ratio of training instances
  )
  
  # Train model with inner CV and grid search
  xgb_model <- train(
    OutcomeBinary ~ .,
    data = train_data,
    method = "xgbTree",
    metric = "ROC",
    trControl = inner_ctrl,
    tuneGrid = tune_grid
  )
  
  # Predict on outer test set
  preds <- predict(xgb_model, newdata = test_data)
  probs <- predict(xgb_model, newdata = test_data, type = "prob")
  
  obs <- factor(test_data$OutcomeBinary, levels = c("Alive", "Dead"))
  pred <- factor(preds, levels = c("Alive", "Dead"))
  
  # AUROC
  roc_obj <- roc(obs, probs$Dead)
  auc_vals[i] <- as.numeric(auc(roc_obj))
  
  # F1-scores
  f1_dead_vals[i] <- F_meas(pred, obs, relevant = "Dead")
  f1_alive_vals[i] <- F_meas(pred, obs, relevant = "Alive")
  
  # Balanced Accuracy
  cm <- confusionMatrix(pred, obs)
  bal_acc_vals[i] <- cm$byClass["Balanced Accuracy"]
}

cat("XGBoost Performance (Nested CV, 10 folds):\n")
cat(sprintf("AUROC: %.2f ± %.2f\n", mean(auc_vals, na.rm = TRUE), sd(auc_vals, na.rm = TRUE)))
cat(sprintf("F1-score (Dead): %.2f ± %.2f\n", mean(f1_dead_vals, na.rm = TRUE), sd(f1_dead_vals, na.rm = TRUE)))
cat(sprintf("F1-score (Alive): %.2f ± %.2f\n", mean(f1_alive_vals, na.rm = TRUE), sd(f1_alive_vals, na.rm = TRUE)))
cat(sprintf("Balanced Accuracy: %.2f ± %.2f\n", mean(bal_acc_vals, na.rm = TRUE), sd(bal_acc_vals, na.rm = TRUE)))




# CatBoost ----
library(reticulate)
library(caret)
library(pROC)
library(dplyr)


# Create numeric outcome for Python
modelling_data$OutcomeBinaryNum <- ifelse(modelling_data$OutcomeBinary == "Dead", 1L, 0L)

# Identify peak columns (those that are purely numeric names)
peak_cols <- grep("^\\d+$", colnames(modelling_data), value = TRUE)

# Build input data: peaks + numeric outcome
peaks_input <- modelling_data[, c(peak_cols, "OutcomeBinaryNum")]

# --- 4. Pass data to Python ---
df_py <- r_to_py(peaks_input)
py$df <- df_py


py_run_string("
import pandas as pd
from catboost import CatBoostClassifier, Pool
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score, f1_score, balanced_accuracy_score

# Use DataFrame from R
X = df.drop(columns=['OutcomeBinaryNum'])
y = df['OutcomeBinaryNum']

skf = StratifiedKFold(n_splits=10, shuffle=True, random_state=42)

auc_vals = []
f1_dead_vals = []
f1_alive_vals = []
bal_acc_vals = []

for train_idx, test_idx in skf.split(X, y):
    X_train, X_test = X.iloc[train_idx], X.iloc[test_idx]
    y_train, y_test = y.iloc[train_idx], y.iloc[test_idx]

    train_pool = Pool(X_train, y_train)
    test_pool = Pool(X_test, y_test)

    model = CatBoostClassifier(
        iterations=200,
        learning_rate=0.1,
        depth=6,
        loss_function='Logloss',
        eval_metric='AUC',
        verbose=False,
        random_seed=42
    )

    model.fit(train_pool)

    probs = model.predict_proba(test_pool)[:, 1]
    preds = (probs > 0.5).astype(int)

    auc_vals.append(roc_auc_score(y_test, probs))
    f1_dead_vals.append(f1_score(y_test, preds, pos_label=1))
    f1_alive_vals.append(f1_score(y_test, preds, pos_label=0))
    bal_acc_vals.append(balanced_accuracy_score(y_test, preds))
")


cat("CatBoost Performance (10-fold CV):\n")
cat(sprintf("AUROC: %.2f ± %.2f\n", mean(py$auc_vals), sd(py$auc_vals)))
cat(sprintf("F1-score (Dead): %.2f ± %.2f\n", mean(py$f1_dead_vals), sd(py$f1_dead_vals)))
cat(sprintf("F1-score (Alive): %.2f ± %.2f\n", mean(py$f1_alive_vals), sd(py$f1_alive_vals)))
cat(sprintf("Balanced Accuracy: %.2f ± %.2f\n", mean(py$bal_acc_vals), sd(py$bal_acc_vals)))


# SVM ----
library(caret)
library(pROC)
library(dplyr)
library(e1071)

# find peak columns
peak_cols <- grep("^\\d+$", colnames(modelling_data), value = TRUE)

# create dataset
svm_data <- modelling_data[, c("OutcomeBinary", peak_cols)]
svm_data$OutcomeBinary <- factor(svm_data$OutcomeBinary, levels = c("Alive", "Dead"))

# create 10 folds
folds <- createFolds(svm_data$OutcomeBinary, k = 10, returnTrain = TRUE)

# store results
auc_vals <- c()
f1_dead_vals <- c()
f1_alive_vals <- c()
bal_acc_vals <- c()

# loop over folds
for (i in seq_along(folds)) {
  cat(sprintf("processing fold %d...\n", i))
  
  train_idx <- folds[[i]]
  train_data <- svm_data[train_idx, ]
  test_data <- svm_data[-train_idx, ]
  
  if (any(table(test_data$OutcomeBinary) == 0)) {
    warning(sprintf("fold %d skipped: missing class in test set", i))
    next
  }
  
  ctrl <- trainControl(
    method = "cv",
    number = 5,
    classProbs = TRUE,
    summaryFunction = twoClassSummary
  )
  
  svm_grid <- expand.grid(
    C = 2^(-2:2),
    sigma = 2^(-5:-1)
  )
  
  svm_model <- train(
    OutcomeBinary ~ .,
    data = train_data,
    method = "svmRadial",
    metric = "ROC",
    trControl = ctrl,
    tuneGrid = svm_grid,
    preProcess = c("center", "scale")
  )
  
  probs <- predict(svm_model, newdata = test_data, type = "prob")
  preds <- predict(svm_model, newdata = test_data)
  
  obs <- factor(test_data$OutcomeBinary, levels = c("Alive", "Dead"))
  pred <- factor(preds, levels = c("Alive", "Dead"))
  
  roc_obj <- roc(obs, probs$Dead)
  auc_vals[i] <- as.numeric(auc(roc_obj))
  f1_dead_vals[i] <- F_meas(pred, obs, relevant = "Dead")
  f1_alive_vals[i] <- F_meas(pred, obs, relevant = "Alive")
  bal_acc_vals[i] <- confusionMatrix(pred, obs)$byClass["Balanced Accuracy"]
}

# show results
cat("svm radial performance (10-fold cv):\n")
cat(sprintf("auroc: %.2f ± %.2f\n", mean(auc_vals, na.rm = TRUE), sd(auc_vals, na.rm = TRUE)))
cat(sprintf("f1-score (dead): %.2f ± %.2f\n", mean(f1_dead_vals, na.rm = TRUE), sd(f1_dead_vals, na.rm = TRUE)))
cat(sprintf("f1-score (alive): %.2f ± %.2f\n", mean(f1_alive_vals, na.rm = TRUE), sd(f1_alive_vals, na.rm = TRUE)))
cat(sprintf("balanced accuracy: %.2f ± %.2f\n", mean(bal_acc_vals, na.rm = TRUE), sd(bal_acc_vals, na.rm = TRUE)))
# BART ----

library(bartMachine)  # Required me to install JAVA
library(caret)
library(pROC)
library(dplyr)

# for reproducibility
set_bart_machine_num_cores(1)  # you can increase this for speed if your system allows

# subset to OutcomeBinary + peaks only
peak_cols <- grep("^\\d+$", colnames(modelling_data), value = TRUE)
bart_df <- modelling_data[, c("OutcomeBinary", peak_cols)]
bart_df$OutcomeBinary <- factor(bart_df$OutcomeBinary, levels = c("Alive", "Dead"))

# set up 10-fold CV
folds <- createFolds(bart_df$OutcomeBinary, k = 10, returnTrain = TRUE)

auc_vals <- c()
f1_dead_vals <- c()
f1_alive_vals <- c()
bal_acc_vals <- c()

for (i in seq_along(folds)) {
  cat(sprintf("processing fold %d...\n", i))
  
  train_idx <- folds[[i]]
  train_data <- bart_df[train_idx, ]
  test_data <- bart_df[-train_idx, ]
  
  if (any(table(test_data$OutcomeBinary) == 0)) {
    warning(sprintf("fold %d skipped due to missing class in test set", i))
    next
  }
  
  # bartMachine expects numeric outcome (0/1)
  train_data$Y <- ifelse(train_data$OutcomeBinary == "Dead", 1, 0)
  test_data$Y <- ifelse(test_data$OutcomeBinary == "Dead", 1, 0)
  
  train_data$OutcomeBinary <- NULL
  test_data$OutcomeBinary <- NULL
  
  model <- bartMachine::bartMachine(
    X = train_data[, !(names(train_data) %in% "Y")],
    y = train_data$Y,
    num_trees = 50,
    num_burn_in = 200,
    num_iterations_after_burn_in = 1000,
    verbose = FALSE
  )
  
  preds_prob <- predict(model, test_data[, !(names(test_data) %in% "Y")])
  preds_class <- ifelse(preds_prob > 0.5, "Dead", "Alive")
  
  obs <- factor(ifelse(test_data$Y == 1, "Dead", "Alive"), levels = c("Alive", "Dead"))
  pred <- factor(preds_class, levels = c("Alive", "Dead"))
  
  roc_obj <- roc(obs, preds_prob)
  auc_vals[i] <- auc(roc_obj)
  
  f1_dead_vals[i] <- F_meas(pred, obs, relevant = "Dead")
  f1_alive_vals[i] <- F_meas(pred, obs, relevant = "Alive")
  
  cm <- confusionMatrix(pred, obs)
  bal_acc_vals[i] <- cm$byClass["Balanced Accuracy"]
}

cat("bart performance (10-fold cv):\n")
cat(sprintf("auroc: %.2f ± %.2f\n", mean(auc_vals, na.rm = TRUE), sd(auc_vals, na.rm = TRUE)))
cat(sprintf("f1-score (dead): %.2f ± %.2f\n", mean(f1_dead_vals, na.rm = TRUE), sd(f1_dead_vals, na.rm = TRUE)))
cat(sprintf("f1-score (alive): %.2f ± %.2f\n", mean(f1_alive_vals, na.rm = TRUE), sd(f1_alive_vals, na.rm = TRUE)))
cat(sprintf("balanced accuracy: %.2f ± %.2f\n", mean(bal_acc_vals, na.rm = TRUE), sd(bal_acc_vals, na.rm = TRUE)))

# Bayesian Regression ----
library(brms)
library(dplyr)

# identify peaks
peak_cols <- grep("^\\d+$", colnames(modelling_data), value = TRUE)

# rename peaks to valid R variable names
new_peak_names <- paste0("peak_", peak_cols)
colnames(modelling_data)[colnames(modelling_data) %in% peak_cols] <- new_peak_names

# update peak column list
peak_cols <- new_peak_names

# select model variables
vars <- c("OutcomeBinary", "Age", "cluster", peak_cols)
bayes_df <- modelling_data %>% select(all_of(vars)) %>% na.omit()

# binary outcome for brms
bayes_df$y <- ifelse(bayes_df$OutcomeBinary == "Dead", 1, 0)

# priors
priors <- c(
  prior(normal(0, 2), class = "b"),
  prior(student_t(3, 0, 2.5), class = "Intercept")
)

# fit model
fit <- brm(
  formula = y ~ .,               
  data = bayes_df %>% select(-OutcomeBinary),
  family = bernoulli(link = "logit"),
  prior = priors,
  chains = 3,
  cores = 3,
  iter = 2000,
  warmup = 200,
  seed = 123
)

# check results
summary(fit)
pp_check(fit)

# predictions
probs <- fitted(fit)[, "Estimate"]
preds <- ifelse(probs > 0.5, "Dead", "Alive")
truth <- bayes_df$OutcomeBinary

# performance
caret::confusionMatrix(factor(preds, levels = c("Alive", "Dead")),
                       factor(truth, levels = c("Alive", "Dead")))





