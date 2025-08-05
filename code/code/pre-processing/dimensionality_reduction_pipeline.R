# Matt Prill
# MSc Applied Data Science
# Dimensionality Reduction pipeline

# Libraries ----
library(tidyverse)
library(caret)
# Labelling (Start of Pipeline 2) ----
meta <- read.csv("data/merged_data/combined_metadata_filtered.csv")  # metadata

# Checking no isolates have been lost
length(meta$sampleID)  # same no. samples
nrow(mat)   # same no. samples

all(rownames(mat) == meta$sampleID)  # Same samples but in a different order

meta <- meta[match(rownames(mat), meta$sampleID), ]  #  Matches them up again
all(rownames(mat) == meta$sampleID)  # proof it's fixed

labels <- meta$Outcome30days  # Extract response variable






# Peak Filtering  ----
# Note: Not part of the pipeline, but first steps before analyses

load("data/merged_data/aligned_peaks.RData")  # called mat

# total peaks before filtering
total_peaks <- ncol(mat)
cat("Total peaks before filtering:", total_peaks, "\n")

# minimum frequency threshold (reference peak present in at least 15% (42) of isolates )
threshold <- 0.15
min_freq <- threshold * nrow(mat) 

# no. isolates with a non-NA (i.e., present) value for each peak
peak_counts <- colSums(!is.na(mat))  

# keep only peaks above the threshold
filtered_matrix <- mat[, peak_counts >= min_freq]

# total peaks after filtering (239)
filtered_peaks_no <- ncol(filtered_matrix)
cat("Peaks after frequency filtering (≥15% of isolates):", filtered_peaks_no, "\n")

# save(filtered_matrix, file = "data/merged_data/filtered_matrix.Rdata")


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


