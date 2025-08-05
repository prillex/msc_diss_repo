set.seed(1234)
# ---- Load packages ----
library(vegan)
library(ggplot2)
library(dplyr)
library(uwot)
library(patchwork)
library(RColorBrewer)

# ---- Load and binarise matrix ----
load("data/merged_data/filtered_matrix.Rdata")
meta <- read.csv("data/merged_data/combined_metadata_filtered.csv")
stopifnot(all(rownames(filtered_matrix) == meta$sampleID))

filtered_matrix_bin <- ifelse(is.na(filtered_matrix), 0, 1)

# ---- Compute Jaccard distance ----
dist_jaccard <- vegdist(filtered_matrix_bin, method = "jaccard")

# ---- NMDS ----
set.seed(1234)
nmds_res <- metaMDS(dist_jaccard, k = 2, trymax = 100)
nmds_scores <- as.data.frame(scores(nmds_res, display = "sites"))
nmds_scores$sampleID <- rownames(nmds_scores)

# ---- UMAP ----
set.seed(1234)
umap_res <- umap(filtered_matrix_bin)
umap_df <- as.data.frame(umap_res)
colnames(umap_df) <- c("UMAP1", "UMAP2")
umap_df$sampleID <- rownames(filtered_matrix_bin)

# ---- Temporary k-means groups for visual separation (no interpretation!) ----
k <- 2  # You can change this based on how many shapes you want
nmds_scores$group <- as.factor(kmeans(nmds_scores[, 1:2], centers = k)$cluster)
umap_df$group <- as.factor(kmeans(umap_df[, 1:2], centers = k)$cluster)

# ---- Generate colour palette ----
palette <- RColorBrewer::brewer.pal(n = k, name = "Set2")
names(palette) <- levels(nmds_scores$group)

# ---- Convex hulls for NMDS ----
hull_nmds <- data.frame()
for (grp in levels(nmds_scores$group)) {
  temp <- nmds_scores[nmds_scores$group == grp, ]
  if (nrow(temp) >= 3) {
    hull_pts <- temp[chull(temp$NMDS1, temp$NMDS2), ]
    hull_nmds <- rbind(hull_nmds, hull_pts)
  }
}

# ---- Convex hulls for UMAP ----
hull_umap <- data.frame()
for (grp in levels(umap_df$group)) {
  temp <- umap_df[umap_df$group == grp, ]
  if (nrow(temp) >= 3) {
    hull_pts <- temp[chull(temp$UMAP1, temp$UMAP2), ]
    hull_umap <- rbind(hull_umap, hull_pts)
  }
}

# ---- NMDS plot ----
p_nmds <- ggplot() +
  geom_polygon(data = hull_nmds, aes(x = NMDS1, y = NMDS2, fill = group, group = group), alpha = 0.3) +
  geom_point(data = nmds_scores, aes(x = NMDS1, y = NMDS2, colour = group), size = 3) +
  scale_colour_manual(values = palette) +
  scale_fill_manual(values = palette) +
  labs(title = "NMDS (Jaccard)", x = "NMDS1", y = "NMDS2") +
  theme_classic() +
  theme(text = element_text(size = 14), legend.position = "none")

# ---- UMAP plot ----
p_umap <- ggplot() +
  geom_polygon(data = hull_umap, aes(x = UMAP1, y = UMAP2, fill = group, group = group), alpha = 0.3) +
  geom_point(data = umap_df, aes(x = UMAP1, y = UMAP2, colour = group), size = 3) +
  scale_colour_manual(values = palette) +
  scale_fill_manual(values = palette) +
  labs(title = "UMAP", x = "UMAP1", y = "UMAP2") +
  theme_classic() +
  theme(text = element_text(size = 14), legend.position = "none")

# ---- Combine plots ----
p_nmds / p_umap










# No colours ----
# ---- Load packages ----
library(vegan)
library(ggplot2)
library(dplyr)
library(uwot)
library(patchwork)

# ---- Load and binarise matrix ----
load("data/merged_data/filtered_matrix.Rdata")
meta <- read.csv("data/merged_data/combined_metadata_filtered.csv")
stopifnot(all(rownames(filtered_matrix) == meta$sampleID))

filtered_matrix_bin <- ifelse(is.na(filtered_matrix), 0, 1)

# ---- Compute Jaccard distance ----
dist_jaccard <- vegdist(filtered_matrix_bin, method = "jaccard")

# ---- NMDS ----
set.seed(1234)
nmds_res <- metaMDS(dist_jaccard, k = 2, trymax = 100)
nmds_scores <- as.data.frame(scores(nmds_res, display = "sites"))
nmds_scores$sampleID <- rownames(nmds_scores)

# ---- UMAP ----
set.seed(1234)
umap_res <- umap(filtered_matrix_bin)
umap_df <- as.data.frame(umap_res)
colnames(umap_df) <- c("UMAP1", "UMAP2")
umap_df$sampleID <- rownames(filtered_matrix_bin)

# ---- NMDS plot (formatted) ----
p_nmds <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(title =  "A: NMDS", x = "\nNMDS1", y = "NMDS2\n") +
  theme_classic() +
  theme(
    text = element_text(size = 18),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    plot.title = element_text(size = 16, face = "bold")
  )

# ---- UMAP plot (formatted) ----
p_umap <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(title = "B: UMAP", x = "\nUMAP1", y = "UMAP2\n") +
  theme_classic() +
  theme(
    text = element_text(size = 18),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    plot.title = element_text(size = 16, face = "bold")
  )

# ---- Combine side-by-side ----
p_nmds | p_umap
















# ---- Load packages ----
library(vegan)
library(ggplot2)
library(dplyr)

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

# ---- Function to generate NMDS plot with hulls for any variable ----
plot_nmds_with_hulls <- function(data, group_var, title_text) {
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
  
  # Plot
  ggplot() +
    geom_polygon(data = hull_data, aes(x = NMDS1, y = NMDS2, fill = group_var, group = group_var), alpha = 0.3) +
    geom_point(data = data, aes(x = NMDS1, y = NMDS2, colour = group_var), size = 3) +
    theme_classic() +
    labs(title = title_text,
         x = "NMDS1", y = "NMDS2", colour = group_var, fill = group_var) +
    theme(
      text = element_text(size = 14),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 12)
    )
}

# ---- Plot 1: Colour by source ----
p1 <- plot_nmds_with_hulls(nmds_data, "source", "NMDS (Jaccard) — Coloured by Source")

# ---- Plot 2: Colour by collection ----
p2 <- plot_nmds_with_hulls(nmds_data, "Collection", "NMDS (Jaccard) — Coloured by Collection")

# ---- Plot 3: Colour by Outcome30days ----
p3 <- plot_nmds_with_hulls(nmds_data, "Outcome30days", "NMDS (Jaccard) — Coloured by 30-day Outcome")

# ---- Display all plots ----
p1
p2
p3









# NMDS with out legends ----


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

# NMDS plotting function ----
plot_nmds <- function(data, group, title, colours = c("#6E4318", "darkblue")) {
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
    labs(title = title, x = "\nNMDS1", y = "NMDS2\n") +
    theme(
      legend.position = "none",
      text = element_text(size = 24),
      axis.text = element_text(size = 26, colour = "black"),
      axis.title = element_text(size = 30, colour = "black"),
      axis.ticks = element_line(colour = "black", linewidth = 1.2),
      axis.line = element_line(colour = "black", linewidth = 1.2)
    )
}

# Generate NMDS plots ----
p1 <- plot_nmds(nmds_data, "Outcome30days", "A)")
p2 <- plot_nmds(nmds_data, "Collection", "B)")
p3 <- plot_nmds(nmds_data, "source", "C)")

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
  scale_colour_manual(values = c("darkblue", "#6E4318", "purple")) +
  theme_classic() +
  labs(title = "D)", x = "\nNMDS1", y = "NMDS2\n") +
  theme(
    legend.position = "none",
    text = element_text(size = 24),
    axis.text = element_text(size = 26, colour = "black"),
    axis.title = element_text(size = 30, colour = "black"),
    axis.ticks = element_line(colour = "black", linewidth = 1.2),
    axis.line = element_line(colour = "black", linewidth = 1.2)
  )

# Combine and save plots ----
nmds_no_legends <- (p1 | p2) / (p3 | p4)

#ggsave("code/code/figures/clustering/nmds_no_legends.png",
#       plot = combined_nmds,
#       width = 20, height = 20, dpi = 600)

# PERMANOVA ----
cat("\nPERMANOVA Results:\n")
print(adonis2(dist_jaccard ~ source, data = meta, permutations = 999))
print(adonis2(dist_jaccard ~ Collection, data = meta, permutations = 999))
print(adonis2(dist_jaccard ~ Outcome30days, data = meta, permutations = 999))










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
      legend.position = "right",
      legend.justification = c(0, 1),
      text = element_text(size = 18),
      legend.text = element_text(size = 16),
      legend.title = element_text(size = 18)
    )
}

# ---- Plot 1: Colour by source ----
p1 <- plot_nmds_with_hulls(nmds_data, "source", "A)", "Source")

# ---- Plot 2: Colour by collection ----
p2 <- plot_nmds_with_hulls(nmds_data, "Collection", "B) ", "Collection")

# ---- Plot 3: Colour by Outcome30days ----
p3 <- plot_nmds_with_hulls(nmds_data, "Outcome30days", "C) ", "Outcome (30 Days)")

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
  labs(title = paste0("D)"),
       x = "\nNMDS1", y = "NMDS2\n", colour = "Cluster") +
  scale_colour_manual(values = c("darkblue", "#6E4318", "purple", "darkgreen")[1:optimal_k]) +
  theme(
    legend.position = "right",
    legend.justification = c(0, 1),
    text = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 18)
  )

# ---- Combine all plots ----
library(patchwork)
(combined_nmds <- (p2 | p3) / (p1 | p4))

#ggsave(
#  filename = "code/code/figures/clustering/nmds_combined_plot.png", 
#  plot = combined_nmds,
#  width = 14,       
#  height = 12,     
#  dpi = 300         # high resolution
#)




ggsave(
filename = "code/code/figures/presentations/combined_nmds_transparent.png",
plot = combined_nmds,
bg = "transparent",
width = 14, height = 14, dpi = 300, limitsize = FALSE)







# Permanova ----
print(adonis2(dist_jaccard ~ source, data = meta, permutations = 999, method = "jaccard"))
print(adonis2(dist_jaccard ~ Collection, data = meta, permutations = 999, method = "jaccard"))
print(adonis2(dist_jaccard ~ Outcome30days, data = meta, permutations = 999, method = "jaccard"))



