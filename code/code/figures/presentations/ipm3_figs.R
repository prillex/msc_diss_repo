# Plot for presentation

load("data/merged_data/combined_spectra_raw_filtered.RData")
load("data/merged_data/smoothed_spectra.RData")

isolate_id <- "A060"

raw <- combined_spectra_raw_filtered[[isolate_id]]
smoothed <- smoothed_spectra[[isolate_id]]

beforeafter_df <- data.frame(mz = raw$mass, intensity = raw$intensity)
after_df <- data.frame(mz = smoothed$mass, intensity = smoothed$intensity)

# Generate plot
(overlay_plot <- ggplot() +
    geom_line(data = beforeafter_df, aes(x = mz, y = intensity), colour = "black", linewidth = 0.5) +
    geom_line(data = after_df, aes(x = mz, y = intensity), colour = "red", linewidth = 1) +
    xlab("\nMass to charge ratio (m/z)") +
    ylab("Signal Intensity\n") +
    xlim(1950, 8000) +
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
      plot.margin = unit(c(1, 1, 1, 1), "cm")
    ))

ggsave(overlay_plot, filename = "code/code/figures/presentations/MS_smoothed.png", bg = "transparent", 
       width = 35, height = 14, dpi = 300)  # Try adjusting these values












# Savitkzy
# Load required libraries
library(ggplot2)
library(signal)  # for Savitzky-Golay filter

# Load data
load("data/merged_data/combined_spectra_raw_filtered.RData")
isolate_id <- "A060"

# Extract raw spectrum
raw <- combined_spectra_raw_filtered[[isolate_id]]

# Apply Savitzky-Golay filter
# Choose a window size (must be odd) and polynomial order
window_size <- 21  # adjust as needed; must be odd
poly_order <- 3     # order of polynomial for fitting

# Apply SG filter to intensity
sg_intensity <- sgolayfilt(raw$intensity, p = poly_order, n = window_size)

# Create data frames for ggplot
beforeafter_df <- data.frame(mz = raw$mass, intensity = raw$intensity)
after_df <- data.frame(mz = raw$mass, intensity = sg_intensity)

# Generate plot
(overlay_plot <- ggplot() +
    geom_line(data = beforeafter_df, aes(x = mz, y = intensity), colour = "black", linewidth = 0.5) +
    geom_line(data = after_df, aes(x = mz, y = intensity), colour = "red", linewidth = 1) +
    xlab("\nMass to charge ratio (m/z)") +
    ylab("Signal Intensity\n") +
    xlim(1950, 8000) +
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
      plot.margin = unit(c(1, 1, 1, 1), "cm")
    ))

# Save plot
ggsave(overlay_plot,
       filename = "code/code/figures/presentations/MS_smoothed_SG.png",
       bg = "transparent",
       width = 30, height = 14, dpi = 300, limitsize = FALSE)


