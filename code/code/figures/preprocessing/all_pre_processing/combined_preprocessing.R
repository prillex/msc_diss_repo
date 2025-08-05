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


# Define nice, consistent theme (no transparency)
pretty_theme <- theme_minimal() +
  theme(
    axis.title.x = element_text(size = 20, colour = "black"),
    axis.title.y = element_text(size = 35, colour = "black"),
    axis.text.x = element_text(size = 25, colour = "black"),
    axis.text.y = element_text(size = 16, colour = "black"),
    panel.grid = element_blank(),
    axis.line = element_line(colour = "black", linewidth = 1.2),
    plot.title = element_text(size = 22, hjust = 0.5),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA)
  )

# Set x-axis limits: from 2000 to max of raw spectrum
xlims <- c(2000, 10000)
ylims <- c(0,2000)

# Individual panels
p1 <- ggplot() +
  geom_line(aes(x = raw$mass, y = raw$intensity), colour = "black") +
  labs(title = "Raw", x = NULL, y = "Absolute Intensity\n") +
  xlim(xlims) + 
  ylim(ylims) +
  pretty_theme

p2 <- ggplot() +
  geom_line(aes(x = smoothed$mass, y = smoothed$intensity), colour = "black") +
  labs(title = "Smoothed", x = NULL, y = "Absolute Intensity\n") +
  xlim(xlims) +
  ylim(ylims) +
  pretty_theme 

p3 <- ggplot() +
  geom_line(aes(x = scaled$mass, y = scaled$intensity), colour = "black") +
  labs(title = "Baseline Corrected & Scaled", x = NULL, y = " Relative Intensity\n") +
  xlim(xlims) +
  pretty_theme 

p4 <- ggplot() +
  geom_segment(
    aes(x = peaks$mz, xend = peaks$mz, y = 0, yend = peaks$intensity),
    colour = "darkblue", alpha = 0.7
  ) +
  geom_point(
    aes(x = peaks$mz, y = peaks$intensity),
    colour = "darkblue", size = 1.5
  ) +
  labs(title = "Selected Peaks", x = NULL, y = "Relative Intensity\n") +
  xlim(xlims) +
  pretty_theme 

p5 <- ggplot(df_filtered, aes(x = mz, y = intensity)) +
  geom_segment(aes(x = mz, xend = mz, y = 0, yend = intensity), 
               colour = "#6E4318", alpha = 0.7) +
  geom_point(colour = "#6E4318", size = 1.5) +
  labs(title = "Filtered Peaks", x = "m/z", y = "Realtive Intensity\n") +
  xlim(xlims) +
  pretty_theme

# Combine vertically
final_plot <- (p1 / p2 / p3 / p4 / p5) +
  plot_annotation(title = "Preprocessing Pipeline: Isolate A060") &
  theme(plot.title = element_text(size = 26, hjust = 0.5))

# Show combined plot
final_plot



# Save individual plots
p1
ggsave(
  filename = "code/code/figures/preprocessing/all_pre_processing/raw.png",
  plot = p1,
  width = 16, height = 4, dpi = 300, limitsize = FALSE
)


p2
p2 <- p2 + theme(plot.margin = margin(t = 10, r = 10, b = 20, l = 10))

ggsave(
  filename = "code/code/figures/preprocessing/all_pre_processing/smoothed.png",
  plot = p2,
  width = 16, height = 4, dpi = 300, limitsize = FALSE
)


p3
ggsave(
  filename = "code/code/figures/preprocessing/all_pre_processing/corrected_scaled.png",
  plot = p3,
  width = 16, height = 4, dpi = 300, limitsize = FALSE
)



p4
ggsave(
  filename = "code/code/figures/preprocessing/all_pre_processing/selected_peaks.png",
  plot = p4,
  width = 16, height = 4, dpi = 300, limitsize = FALSE
)

p5
ggsave(
  filename = "code/code/figures/preprocessing/all_pre_processing/filtered_peaks.png",
  plot = p5,
  width = 16, height = 4, dpi = 300, limitsize = FALSE
)
