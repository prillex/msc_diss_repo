# Matt Prill
# Toxicity figs



# Random forest
# libraries ----
library(tidyverse)
library(randomForest)

# Data ----
load("data/merged_data/selected_peaks_toxicity.Rdata")
meta <- read.csv("data/merged_data/combined_metadata_filtered.csv")  # Meta


# Create a filter for samples that do NOT start with 'A' followed by digits (i.e. cc22)
keep <- !grepl("^A\\d+$", meta$sampleID)

# Subset meta
meta_cc22 <- meta[keep, ]

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




# Importance plot ----
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
      axis.text.y = element_text(size = 40, face = "bold", color = "black"),
      axis.text.x = element_text(size = 40, face = "bold", color = "black"),
      axis.title.y = element_text(size = 40, face = "bold", color = "black"),
      axis.title.x = element_text(size = 40, face = "bold", color = "black"),
      plot.title = element_text(size = 18, face = "bold", color = "black"),
      panel.background = element_rect(fill = "transparent"),
      plot.background = element_rect(fill = "transparent", colour = NA),
      plot.margin = grid::unit(c(20, 20, 20, 20), "pt")
    ))



# save
#ggsave(rf_importance,
#       filename = "code/code/figures/presentations/rf_importance.png",
#       bg = "transparent",
#       width = 30, height = 14, dpi = 300, limitsize = FALSE)




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
                       mainbar.y.label = "    Intersection Count",
                       sets.x.label = "Peak Occurrence Count",
                       text.scale = c(5, 5, 5, 5, 4, 5)))





# Toxicity ----
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
              width = 0.1, height = 0, size = 6, alpha = 0.9) +
  
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
    axis.title.x = element_text(colour = "black", size = 40),
    axis.title.y = element_text(colour = "black", size = 40),
    axis.text.x = element_text(colour = "black", size = 40, angle = 45, hjust = 1),
    axis.text.y = element_text(colour = "black", size = 40),
    axis.line = element_line(color = "black"),
    panel.background = element_rect(fill = "transparent"),
    plot.background = element_rect(fill = "transparent", colour = NA),
    plot.title = element_text(size = 5, face = "bold", color = "black"),
    plot.margin = unit(c(20, 15, 15, 15), "pt")
  )

print(tox_combo_plot)



# save
#ggsave(tox_combo_plot,
#       filename = "code/code/figures/presentations/tox_combo_plot.png",
#       bg = "transparent",
#       width = 30, height = 14, dpi = 300, limitsize = FALSE)





# Showing that all combinations are associated with high toxicity
df_combined$toxicity <- df_tox$toxicity[match(rownames(df_combined), rownames(df_tox))]

combo_high_percent <- df_combined %>%
  filter(combo_raw != "000") %>%
  mutate(toxicity_num = as.numeric(as.character(toxicity))) %>%  # Convert factor to numeric
  group_by(combo) %>%
  summarise(
    n = n(),
    n_high = sum(toxicity_num == 1, na.rm = TRUE),
    percent_high = round(100 * n_high / n, 1)
  ) %>%
  arrange(desc(percent_high))

# View result
print(combo_high_percent)





# Smoothed Spectra ----
# Load data
load("data/merged_data/scaled_spectra.RData")
load("data/merged_data/filtered_matrix.Rdata")

# Define samples and m/z range
samples <- c("ASARM102", "ASASM12", "ASASM181")
mz_range <- c(2900, 3100)
target_peaks <- c(2977, 3006, 3026)

# Custom colours
custom_colors <- c(
  "Peak 3026" = "#6E4318",   # brown
  "Peak 3006" = "darkblue",  # dark blue
  "Peak 2977" = "purple"     # purple
)

# Match closest actual m/z in the matrix to target peaks
matrix_mz <- as.numeric(colnames(filtered_matrix))
matched_mz <- sapply(target_peaks, function(x) matrix_mz[which.min(abs(matrix_mz - x))])
names(matched_mz) <- paste0("Peak ", target_peaks)

# Extract corresponding intensities from the matrix
peak_points <- lapply(samples, function(sample_id) {
  df <- data.frame(
    mass = as.numeric(matched_mz),
    intensity = as.numeric(filtered_matrix[sample_id, as.character(matched_mz)]),
    sample = sample_id,
    peak = names(matched_mz)  # "Peak 2977", etc.
  )
  df
}) %>% bind_rows() %>% na.omit()



# Prepare combined spectra for plotting
df_list <- lapply(samples, function(sample_id) {
  spec <- scaled_spectra[[sample_id]]
  data.frame(
    mass = spec$mass,
    intensity = spec$intensity,
    sample = sample_id
  ) %>%
    filter(mass >= mz_range[1], mass <= mz_range[2])
})
df_combined <- bind_rows(df_list)



# Plot
(tox_peak_plot <- ggplot(df_combined, aes(x = mass, y = intensity)) +
  geom_line(colour = "black") +
  geom_point(data = peak_points, 
             aes(x = mass, y = intensity, colour = peak),
             size = 10) +
  scale_colour_manual(values = custom_colors, name = "Peak") +
  facet_wrap(~ sample, ncol = 1, scales = "free_y") +
  scale_x_continuous(breaks = seq(2900, 3100, by = 25)) +
  labs(title = "Stacked Scaled Spectra (m/z 2900–3100)",
       x = "\nm/z", y = "Intensity\n") +
    theme(
      axis.title.x = element_text(colour = "black", size = 40),
      axis.title.y = element_text(colour = "black", size = 40),
      axis.text.x = element_text(colour = "black", size = 30, angle = 45, hjust = 1),
      axis.text.y = element_text(colour = "black", size = 40),
      axis.line = element_line(color = "black"),
      panel.background = element_rect(fill = "transparent"),
      plot.background = element_rect(fill = "transparent", colour = NA),
      legend.background = element_rect(fill = "transparent", colour = NA),
      plot.margin = unit(c(20, 25, 15, 15), "pt"),
      strip.text = element_blank(),
      strip.background = element_blank(),
      legend.title = element_text(size = 36),
      legend.text = element_text(size = 32),
      legend.key.size = unit(2, "cm"),
      panel.spacing = unit(3, "lines"),         # more space between stacked plots
    ))


# save
#ggsave(tox_peak_plot,
#       filename = "code/code/figures/presentations/tox_peak_plot.png",
#       bg = "transparent",
#       width = 30, height = 14, dpi = 300, limitsize = FALSE)





# Same plot, not transparent ----
# Plot
(tox_peak_plot_fig <- ggplot(df_combined, aes(x = mass, y = intensity)) +
   geom_line(colour = "black") +
   geom_point(data = peak_points, 
              aes(x = mass, y = intensity, colour = peak),
              size = 10) +
   scale_colour_manual(values = custom_colors, name = "Peak") +
   facet_wrap(~ sample, ncol = 1, scales = "free_y") +
   scale_x_continuous(breaks = seq(2900, 3100, by = 25)) +
   labs(x = "\nm/z", y = "Intensity\n") +
   theme(
     axis.title.x = element_text(colour = "black", size = 40),
     axis.title.y = element_text(colour = "black", size = 40),
     axis.text.x = element_text(colour = "black", size = 30, angle = 45, hjust = 1),
     axis.text.y = element_text(colour = "black", size = 40),
     axis.line = element_line(color = "black"),
     panel.background = element_rect(fill = "white"),  # changed from "transparent"
     plot.background = element_rect(fill = "white", colour = NA),  # changed from "transparent"
     legend.background = element_rect(fill = "white", colour = NA),  # changed from "transparent"
     plot.margin = unit(c(20, 25, 15, 15), "pt"),
     strip.text = element_blank(),
     strip.background = element_blank(),
     legend.title = element_text(size = 36),
     legend.text = element_text(size = 32),
     legend.key.size = unit(2, "cm"),
     panel.spacing = unit(3, "lines")
   ))

# Save (non-transparent background)
ggsave(tox_peak_plot_fig,
       filename = "code/code/figures/toxicity/tox_peak_plot_fig.png",
       bg = "white",  # changed from "transparent"
       width = 30, height = 14, dpi = 300, limitsize = FALSE)

# Raw spectra ----
# Load raw spectra
load("data/merged_data/combined_spectra_raw_filtered.RData")  # loads combined_spectra_raw_filtered
load("data/merged_data/filtered_matrix.Rdata")

# Define samples and m/z range
samples <- c("ASARM102", "ASASM12", "ASASM181")
mz_range <- c(2900, 3100)
target_peaks <- c(2977, 3006, 3026)

# Custom colours
custom_colors <- c(
  "Peak 3026" = "#6E4318",   # brown
  "Peak 3006" = "darkblue",  # dark blue
  "Peak 2977" = "purple"     # purple
)

# Match closest actual m/z in the matrix to target peaks
matrix_mz <- as.numeric(colnames(filtered_matrix))
matched_mz <- sapply(target_peaks, function(x) matrix_mz[which.min(abs(matrix_mz - x))])
names(matched_mz) <- paste0("Peak ", target_peaks)

# Extract corresponding intensities from the matrix
# Extract peak intensities directly from raw spectra
peak_points <- lapply(samples, function(sample_id) {
  spec <- combined_spectra_raw_filtered[[sample_id]]
  
  # For each matched m/z, find closest value in raw spectrum
  peak_df <- lapply(seq_along(matched_mz), function(i) {
    target_mz <- matched_mz[i]
    peak_label <- names(matched_mz)[i]
    
    # Find closest index in raw spectrum
    idx <- which.min(abs(spec$mass - target_mz))
    
    data.frame(
      mass = spec$mass[idx],
      intensity = spec$intensity[idx],
      sample = sample_id,
      peak = peak_label
    )
  }) %>% bind_rows()
  
  peak_df
}) %>% bind_rows()


# Build data from raw spectra
df_list_raw <- lapply(samples, function(sample_id) {
  spec <- combined_spectra_raw_filtered[[sample_id]]
  data.frame(
    mass = spec$mass,
    intensity = spec$intensity,
    sample = sample_id
  ) %>%
    filter(mass >= mz_range[1], mass <= mz_range[2])
})
df_combined_raw <- bind_rows(df_list_raw)

# Plot
(tox_peak_plot_raw <- ggplot(df_combined_raw, aes(x = mass, y = intensity)) +
    geom_line(colour = "black") +
    geom_point(data = peak_points, 
               aes(x = mass, y = intensity, colour = peak),
               size = 10) +
    scale_colour_manual(values = custom_colors, name = "Peak") +
    facet_wrap(~ sample, ncol = 1, scales = "free_y") +
    scale_x_continuous(breaks = seq(2900, 3100, by = 25)) +
    labs(title = "Stacked Raw Spectra (m/z 2900–3100)",
         x = "\nm/z", y = "Intensity\n") +
    theme(
      axis.title.x = element_text(colour = "black", size = 40),
      axis.title.y = element_text(colour = "black", size = 40),
      axis.text.x = element_text(colour = "black", size = 30, angle = 45, hjust = 1),
      axis.text.y = element_text(colour = "black", size = 40),
      axis.line = element_line(color = "black"),
      panel.background = element_rect(fill = "white"),  # changed from "transparent"
      plot.background = element_rect(fill = "white", colour = NA),  # changed from "transparent"
      legend.background = element_rect(fill = "white", colour = NA),  # changed from "transparent"
      plot.margin = unit(c(20, 25, 15, 15), "pt"),
      strip.text = element_blank(),
      strip.background = element_blank(),
      legend.title = element_text(size = 36),
      legend.text = element_text(size = 32),
      legend.key.size = unit(2, "cm"),
      panel.spacing = unit(3, "lines")
    ))


# Save
#ggsave(tox_peak_plot_raw,
#       filename = "code/code/figures/presentations/tox_peak_plot_raw.png",
#       bg = "white",
#       width = 30, height = 14, dpi = 300, limitsize = FALSE)








# ----
# ----
# ----
# Extra ----
# Raw spectra 2 ----
# Load raw spectra
load("data/merged_data/combined_spectra_raw_filtered.RData")  # loads combined_spectra_raw_filtered
load("data/merged_data/filtered_matrix.Rdata")

# Define samples and m/z range
samples <- c("ASARM70", "ASARM193", "ASASM392")
mz_range <- c(2900, 3100)
target_peaks <- c(2977, 3006, 3026)

# Custom colours
custom_colors <- c(
  "Peak 3026" = "#6E4318",   # brown
  "Peak 3006" = "darkblue",  # dark blue
  "Peak 2977" = "purple"     # purple
)

# Match closest actual m/z in the matrix to target peaks
matrix_mz <- as.numeric(colnames(filtered_matrix))
matched_mz <- sapply(target_peaks, function(x) matrix_mz[which.min(abs(matrix_mz - x))])
names(matched_mz) <- paste0("Peak ", target_peaks)

# Extract peak intensities directly from raw spectra
peak_points <- lapply(samples, function(sample_id) {
  spec <- combined_spectra_raw_filtered[[sample_id]]
  
  peak_df <- lapply(seq_along(matched_mz), function(i) {
    target_mz <- matched_mz[i]
    peak_label <- names(matched_mz)[i]
    
    idx <- which.min(abs(spec$mass - target_mz))
    
    data.frame(
      mass = spec$mass[idx],
      intensity = spec$intensity[idx],
      sample = sample_id,
      peak = peak_label
    )
  }) %>% bind_rows()
  
  peak_df
}) %>% bind_rows()

# Build data from raw spectra
df_list_raw <- lapply(samples, function(sample_id) {
  spec <- combined_spectra_raw_filtered[[sample_id]]
  data.frame(
    mass = spec$mass,
    intensity = spec$intensity,
    sample = sample_id
  ) %>%
    filter(mass >= mz_range[1], mass <= mz_range[2])
})
df_combined_raw <- bind_rows(df_list_raw)

# Plot
(tox_peak_plot_raw_2 <- ggplot(df_combined_raw, aes(x = mass, y = intensity)) +
    geom_line(colour = "black") +
    geom_point(data = peak_points %>% filter(sample != "ASASM392"),  # remove peaks from bottom facet
               aes(x = mass, y = intensity, colour = peak),
               size = 10) +
    scale_colour_manual(values = custom_colors, name = "Peak") +
    facet_wrap(~ sample, ncol = 1, scales = "free_y") +
    scale_x_continuous(breaks = seq(2900, 3100, by = 25)) +
    labs(title = "Stacked Raw Spectra (m/z 2900–3100)",
         x = "\nm/z", y = "Intensity\n") +
    theme(
      axis.title.x = element_text(colour = "black", size = 40),
      axis.title.y = element_text(colour = "black", size = 40),
      axis.text.x = element_text(colour = "black", size = 30, angle = 45, hjust = 1),
      axis.text.y = element_text(colour = "black", size = 40),
      axis.line = element_line(color = "black"),
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white", colour = NA),
      legend.background = element_rect(fill = "white", colour = NA),
      plot.margin = unit(c(20, 25, 15, 15), "pt"),
      strip.text = element_blank(),
      strip.background = element_blank(),
      legend.title = element_text(size = 36),
      legend.text = element_text(size = 32),
      legend.key.size = unit(2, "cm"),
      panel.spacing = unit(3, "lines")
    ))

# Save
#ggsave(tox_peak_plot_raw_2,
#       filename = "code/code/figures/presentations/tox_peak_plot_raw_2.png",
#       bg = "white",
#       width = 30, height = 14, dpi = 300, limitsize = FALSE)


# Smoothed_2 ----
# Load data
load("data/merged_data/scaled_spectra.RData")
load("data/merged_data/filtered_matrix.Rdata")

# Define samples and m/z range
samples <- c("ASARM70", "ASARM193", "ASASM392")mz_range <- c(2900, 3100)
target_peaks <- c(2977, 3006, 3026)

# Custom colours
custom_colors <- c(
  "Peak 3026" = "#6E4318",   # brown
  "Peak 3006" = "darkblue",  # dark blue
  "Peak 2977" = "purple"     # purple
)

# Match closest actual m/z in the matrix to target peaks
matrix_mz <- as.numeric(colnames(filtered_matrix))
matched_mz <- sapply(target_peaks, function(x) matrix_mz[which.min(abs(matrix_mz - x))])
names(matched_mz) <- paste0("Peak ", target_peaks)

# Extract corresponding intensities from the matrix
peak_points <- lapply(samples, function(sample_id) {
  df <- data.frame(
    mass = as.numeric(matched_mz),
    intensity = as.numeric(filtered_matrix[sample_id, as.character(matched_mz)]),
    sample = sample_id,
    peak = names(matched_mz)  # "Peak 2977", etc.
  )
  df
}) %>% bind_rows() %>% na.omit()



# Prepare combined spectra for plotting
df_list <- lapply(samples, function(sample_id) {
  spec <- scaled_spectra[[sample_id]]
  data.frame(
    mass = spec$mass,
    intensity = spec$intensity,
    sample = sample_id
  ) %>%
    filter(mass >= mz_range[1], mass <= mz_range[2])
})
df_combined <- bind_rows(df_list)



# Plot
(tox_peak_plot_2 <- ggplot(df_combined, aes(x = mass, y = intensity)) +
    geom_line(colour = "black") +
    geom_point(data = peak_points, 
               aes(x = mass, y = intensity, colour = peak),
               size = 10) +
    scale_colour_manual(values = custom_colors, name = "Peak") +
    facet_wrap(~ sample, ncol = 1, scales = "free_y") +
    scale_x_continuous(breaks = seq(2900, 3100, by = 25)) +
    labs(title = "Stacked Scaled Spectra (m/z 2900–3100)",
         x = "\nm/z", y = "Intensity\n") +
    theme(
      axis.title.x = element_text(colour = "black", size = 40),
      axis.title.y = element_text(colour = "black", size = 40),
      axis.text.x = element_text(colour = "black", size = 30, angle = 45, hjust = 1),
      axis.text.y = element_text(colour = "black", size = 40),
      axis.line = element_line(color = "black"),
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white", colour = NA),
      legend.background = element_rect(fill = "white", colour = NA),
      plot.margin = unit(c(20, 25, 15, 15), "pt"),
      strip.text = element_blank(),
      strip.background = element_blank(),
      legend.title = element_text(size = 36),
      legend.text = element_text(size = 32),
      legend.key.size = unit(2, "cm"),
      panel.spacing = unit(3, "lines"),         # more space between stacked plots
    ))


# save
#ggsave(tox_peak_plot_2,
#       filename = "code/code/figures/presentations/tox_peak_plot_2.png",
#       bg = "white",
#       width = 30, height = 14, dpi = 300, limitsize = FALSE)

