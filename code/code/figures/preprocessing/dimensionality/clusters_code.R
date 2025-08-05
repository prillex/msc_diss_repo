
res <- align_peaks_custom(peak_list, tol = 50)  # tol is the mz within which peaks are the same

# summarise
cat("bins:", length(res$bins), "\n")             # no. distinct features at this point
cat("total peaks:", sum(!is.na(res$mat)), "\n")  # no. peaks assinged to the bins above
print(rowSums(!is.na(res$mat)))

# assign for saving
bins <- res$bins
mat  <- res$mat

meta <- meta[match(rownames(mat), meta$sampleID), ]  #  Matches them up again


threshold <- 0.15
min_freq <- threshold * nrow(mat) 

# no. isolates with a non-NA (i.e., present) value for each peak
peak_counts <- colSums(!is.na(mat))  

# keep only peaks above the threshold
filtered_matrix <- mat[, peak_counts >= min_freq]

filtered_matrix_bin <- ifelse(is.na(filtered_matrix), 0, 1)  # remove na, and make binary

nrow(filtered_matrix_bin)


umap_res <- umap(filtered_matrix_bin)

plot(umap_res$layout)

# For Source
plot(umap_res$layout, col = as.factor(meta$source), pch = 19,
     main = "UMAP of Filtered Peaks (Source) (50 Da Tolerance)", xlab = "UMAP1", ylab = "UMAP2")
legend("topright", legend = levels(as.factor(meta$source)),
       col = 1:length(unique(meta$source)), pch = 19)








filtered_matrix_bin


filter(filtered_matrix_bin, source == "A")