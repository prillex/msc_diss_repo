# NMDS colouring by meta ----
# ---- Load packages ----
library(vegan)
library(ggplot2)
library(dplyr)
library(patchwork)

# ---- Load data ----
load("data/merged_data/filtered_matrix.Rdata")
meta <- read.csv("data/merged_data/combined_metadata_filtered.csv")
filtered_matrix_bin <- ifelse(is.na(filtered_matrix), 0, 1)

# ---- Compute Jaccard distance and NMDS ----
dist_jaccard <- vegdist(filtered_matrix_bin, method = "jaccard")
set.seed(1234)
nmds_res <- metaMDS(dist_jaccard, k = 2, trymax = 100)
nmds_scores <- as.data.frame(scores(nmds_res, display = "sites"))
nmds_scores$sampleID <- rownames(nmds_scores)
nmds_data <- left_join(nmds_scores, meta, by = "sampleID")




plot_nmds_plain <- function(data, title_text = "") {
ggplot(data, aes(x = NMDS1, y = NMDS2)) +
  geom_point(size = 3, colour = "black") +
  theme_classic() +
  labs(title = title_text, x = "\nNMDS1", y = "NMDS2\n") +
  theme(
    axis.title.x = element_text(colour = "black", size = 70),
    axis.title.y = element_text(colour = "black", size = 70),
    axis.text.x = element_text(colour = "black", size = 50),
    axis.text.y = element_text(colour = "black", size = 50),
    axis.line = element_line(colour = "black", linewidth = 2),
    panel.background = element_rect(fill = "transparent"),
    plot.background = element_rect(fill = "transparent", colour = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = unit(c(1, 1, 1, 1), "cm"),
    legend.position = "none"
  )
}

# Create and show the plot
(p0 <- plot_nmds_plain(nmds_data, title_text = ""))












# ---- Function to generate NMDS plot with 2-colour scheme ----
plot_nmds_with_hulls <- function(data, group_var, title_text, legend_title) {
  data$group_var <- as.factor(data[[group_var]])
  
  # Compute convex hulls
  hull_data <- data.frame()
  for (grp in unique(data$group_var)) {
    temp <- data[data$group_var == grp, ]
    if (nrow(temp) >= 3) {
      hull_pts <- temp[chull(temp$NMDS1, temp$NMDS2), ]
      hull_data <- rbind(hull_data, hull_pts)
    }
  }
  
  ggplot() +
    geom_polygon(data = hull_data, aes(x = NMDS1, y = NMDS2, fill = group_var, group = group_var), alpha = 0.3) +
    geom_point(data = data, aes(x = NMDS1, y = NMDS2, colour = group_var), size = 3) +
    scale_colour_manual(values = c("#6E4318", "darkblue")) +
    scale_fill_manual(values = c("#6E4318", "darkblue")) +
    theme_classic() +
    labs(title = title_text,
         x = "\nNMDS1", y = "NMDS2\n", colour = legend_title, fill = legend_title) +
    theme(
      axis.title.x = element_text(colour = "black", size = 70),
      axis.title.y = element_text(colour = "black", size = 70),
      axis.text.x = element_text(colour = "black", size = 50),
      axis.text.y = element_text(colour = "black", size = 50),
      axis.line = element_line(colour = "black", linewidth = 2),
      panel.background = element_rect(fill = "transparent"),
      plot.background = element_rect(fill = "transparent", colour = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = unit(c(1, 1, 1, 1), "cm"),
      legend.key.size = unit(3, 'cm'), #change legend key size
      legend.key.height = unit(2, 'cm'), #change legend key height
      legend.key.width = unit(5, 'cm'), #change legend key width
      legend.title = element_text(size = 20), #change legend title font size
      legend.text = element_text(size = 20)
    )
}



# ---- Plot 1: Colour by source ----
(p1 <- plot_nmds_with_hulls(nmds_data, "source", "", "Source"))

# ---- Plot 2: Colour by collection ----
(p2 <- plot_nmds_with_hulls(nmds_data, "Collection", " ", "Collection"))

# ---- Plot 3: Colour by Outcome30days ----
(p3 <- plot_nmds_with_hulls(nmds_data, "Outcome30days", " ", "Outcome (30 Days)"))

# ---- Display all plots ----
p1
p2
p3

(p2 | p3) / p1

# ---- Clustering ----
library(cluster)
library(factoextra)

# Try 2 to 5 clusters and compute average silhouette widths
sil_widths <- sapply(2:5, function(k) {
  cluster_assignments <- cutree(hclust(dist_jaccard), k = k)
  silhouette_obj <- silhouette(cluster_assignments, dist_jaccard)
  mean(silhouette_obj[, 3])
})

# Plot silhouette scores 
plot(2:5, sil_widths, type = "b", pch = 19,
     xlab = "Number of Clusters", ylab = "Average Silhouette Width",
     main = "Silhouette Analysis (Jaccard)")

# Pick best number of clusters
optimal_k <- which.max(sil_widths) + 1  # because index starts at 1 for k=2
cat("Optimal number of clusters based on silhouette width:", optimal_k, "\n")

# Assign clusters
nmds_data$cluster <- as.factor(cutree(hclust(dist_jaccard), k = optimal_k))

# ---- Plot 4: NMDS coloured by clustering ----
p4 <- ggplot(nmds_data, aes(x = NMDS1, y = NMDS2, colour = cluster)) +
  geom_point(size = 3) +
  theme_classic() +
  labs(title = paste0(""),
       x = "\nNMDS1", y = "NMDS2\n", colour = "Cluster") +
  scale_colour_manual(values = c("darkblue", "#6E4318", "purple", "darkgreen")[1:optimal_k]) +
  theme(
    axis.title.x = element_text(colour = "black", size = 70),
    axis.title.y = element_text(colour = "black", size = 70),
    axis.text.x = element_text(colour = "black", size = 50),
    axis.text.y = element_text(colour = "black", size = 50),
    axis.line = element_line(colour = "black", linewidth = 2),
    panel.background = element_rect(fill = "transparent"),
    plot.background = element_rect(fill = "transparent", colour = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = unit(c(1, 1, 1, 1), "cm"),
    legend.key.size = unit(3, 'cm'), #change legend key size
    legend.key.height = unit(2, 'cm'), #change legend key height
    legend.key.width = unit(5, 'cm'), #change legend key width
    legend.title = element_text(size = 20), #change legend title font size
    legend.text = element_text(size = 20)
  )

# ---- Combine all plots ----
library(patchwork)
(combined_nmds <- (p2 | p3) / (p1 | p4))


ggsave(
  filename = "code/code/figures/presentations/source_nmds_transparent.png",
  plot = p1,
  bg = "transparent",
  width = 20, height = 14, dpi = 300, limitsize = FALSE)


ggsave(
  filename = "code/code/figures/presentations/base_nmds_transparent.png",
  plot = p0,
  bg = "transparent",
  width = 20, height = 14, dpi = 300, limitsize = FALSE)


ggsave(
  filename = "code/code/figures/presentations/collection_nmds_transparent.png",
  plot = p2,
  bg = "transparent",
  width = 20, height = 14, dpi = 300, limitsize = FALSE)

ggsave(
  filename = "code/code/figures/presentations/outcome_nmds_transparent.png",
  plot = p3,
  bg = "transparent",
  width = 20, height = 14, dpi = 300, limitsize = FALSE)

ggsave(
  filename = "code/code/figures/presentations/cluster_nmds_transparent.png",
  plot = p4,
  bg = "transparent",
  width = 20, height = 14, dpi = 300, limitsize = FALSE)






# Pre-processing figs ----
library(ggplot2)
library(patchwork)

# Data (to avoid running the entire pipeline)
# load required data
load("data/merged_data/combined_spectra_raw_filtered.RData")
load("data/merged_data/smoothed_spectra.RData")
load("data/merged_data/baseline_corrected_spectra.RData")  # baseline corrected
load("data/merged_data/scaled_spectra.RData")
load("data/merged_data/peak_list.RData")
load("data/merged_data/filtered_matrix.Rdata")

# isolate A060 data
raw <- combined_spectra_raw_filtered[["A060"]]
smoothed <- smoothed_spectra[["A060"]]
corrected <- baseline_corrected_spectra[["A060"]]
scaled <- scaled_spectra[["A060"]]
peaks <- peak_list[["A060"]]
filtered_peaks  <- filtered_matrix["A060", ]
df_filtered <- data.frame(
  mz = as.numeric(names(filtered_peaks)),
  intensity = as.numeric(filtered_peaks)
)
df_filtered <- na.omit(df_filtered)  # Optional: remove NAs if present

# Define common theme
library(ggplot2)
library(patchwork)

# Define common theme
custom_theme <- theme_minimal() +
  theme(
    axis.title.x = element_text(size = 20, colour = "black"),
    axis.title.y = element_text(size = 20, colour = "black"),
    axis.text.x = element_text(size = 16, colour = "black"),
    axis.text.y = element_text(size = 16, colour = "black"),
    panel.background = element_rect(fill = "transparent", colour = NA),
    plot.background = element_rect(fill = "transparent", colour = NA),
    panel.grid = element_blank(),
    axis.line = element_line(colour = "black", linewidth = 1.2),
    plot.title = element_text(size = 22, hjust = 0.5)
  )

# Set x-axis limits: from 2000 to max of raw spectrum
xlims <- c(2000, max(raw$mass, na.rm = TRUE))

# Raw spectrum
p1 <- ggplot() +
  geom_line(aes(x = raw$mass, y = raw$intensity), colour = "black") +
  labs(title = "Raw", x = NULL, y = "Intensity") +
  xlim(xlims) +
  custom_theme +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank())

# Selected peaks
p2 <- ggplot() +
  geom_segment(
    aes(x = peaks$mz, xend = peaks$mz, y = 0, yend = peaks$intensity),
    colour = "purple", alpha = 0.7
  ) +
  geom_point(
    aes(x = peaks$mz, y = peaks$intensity),
    colour = "purple", size = 1.5
  ) +
  labs(title = "Selected Peaks", x = NULL, y = "Intensity") +
  xlim(xlims) +
  custom_theme +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank())

# Filtered peaks
p3 <- ggplot(df_filtered, aes(x = mz, y = intensity)) +
  geom_segment(aes(x = mz, xend = mz, y = 0, yend = intensity), 
               colour = "darkorange", alpha = 0.7) +
  geom_point(colour = "darkorange", size = 1.5) +
  labs(title = "Filtered Peaks (Post-Frequency Threshold)", x = "m/z", y = "Intensity") +
  xlim(xlims) +
  custom_theme

# Combine vertically (apply theme directly into plot object)
pre_p_plot <- (p1 / p2 / p3) +
  plot_annotation(title = "Preprocessing Pipeline: Isolate A060") &
  theme(
    plot.title = element_text(size = 26, hjust = 0.5),
    plot.background = element_rect(fill = "transparent", colour = NA)
  )

# Save with transparent background
ggsave(
  filename = "code/code/figures/presentations/preprocessing.png",
  plot = pre_p_plot,
  bg = "transparent",  # <- critical
  width = 36, height = 14, dpi = 300, limitsize = FALSE
)
