# Example mass spectra
# Transparent

load("data/merged_data/combined_spectra_filtered.RData")  # called combined_spectra_filtered

library(ggplot2)



# Example
# Convert the spectrum data to a data frame
spectrum_df <- data.frame(
  mz = combined_spectra_filtered[["A030"]]$mass,
  intensity = combined_spectra_filtered[["A030"]]$intensity
)

# Create the ggplot
(spectrum_plot <- ggplot(spectrum_df, aes(x = mz, y = intensity)) +
  geom_line(colour = "black", linewidth = 0.5) +
  xlab("\nMass to charge ratio (m/z)") +
  ylab("Signal Intensity\n") +
  xlim(1950, 10000) +
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

#ggsave(spectrum_plot, filename = "code/code/figures/presentations/spectrum_plot.png", bg = "transparent", 
#       width = 35, height = 14, dpi = 300)  # Try adjusting these values



# Baseline correction
poor_baseline <- data.frame(
  mz = combined_spectra_filtered[["A022"]]$mass,
  intensity = combined_spectra_filtered[["A022"]]$intensity
)

# Create the ggplot
(poor_baseline <- ggplot(poor_baseline, aes(x = mz, y = intensity)) +
    geom_line(colour = "black", linewidth = 0.5) +
    xlab("\nMass to charge ratio (m/z)") +
    ylab("Signal Intensity\n") +
    xlim(1950, 10000) +
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


#ggsave(poor_baseline, filename = "code/code/figures/presentations/poor_baseline.png", bg = "transparent", 
#       width = 35, height = 14, dpi = 300)  # Try adjusting these values



# Local Maxima

# Baseline correction
local_maxima <- data.frame(
  mz = combined_spectra_filtered[["A022"]]$mass,
  intensity = combined_spectra_filtered[["A022"]]$intensity
)

# Create the ggplot
(local_maxima <- ggplot(local_maxima, aes(x = mz, y = intensity)) +
    geom_line(colour = "black", linewidth = 0.5) +
    xlab("\nMass to charge ratio (m/z)") +
    ylab("Signal Intensity\n") +
    xlim(4000, 6000) +
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



# Wavelet transform
# Before
# Convert the spectrum data to a data frame
beforeafter_df <- data.frame(
  mz = combined_spectra_filtered[["ASASM430"]]$mass,
  intensity = combined_spectra_filtered[["ASASM430"]]$intensity
)

# Create the ggplot
(before_plot <- ggplot(beforeafter_df, aes(x = mz, y = intensity)) +
    geom_line(colour = "black", linewidth = 0.5) +
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



#ggsave(before_plot, filename = "code/code/figures/presentations/ASASM430_before.png", bg = "transparent", 
#       width = 35, height = 14, dpi = 300)  # Try adjusting these values

# After

# Install if not already installed
if (!requireNamespace("MassSpecWavelet", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("MassSpecWavelet")
}

library(MassSpecWavelet)
library(ggplot2)

# Ensure mz and intensity are correctly defined
mz <- beforeafter_df$mz
intensity <- beforeafter_df$intensity

cwt_result <- cwt(intensity, scales = 1:64, wavelet = "mexh")
dim(cwt_result)  # 21450 x 64

# Smooth using high scales
smoothed_signal <- rowSums(Re(cwt_result[, 50:64, drop = FALSE]))
smoothed_signal <- pmax(smoothed_signal, 0)
smoothed_signal_norm <- (smoothed_signal - min(smoothed_signal)) / (max(smoothed_signal) - min(smoothed_signal))

# Confirm lengths match
length(smoothed_signal) == length(beforeafter_df$mz)  # Should be TRUE

# Combine into a dataframe
# Create the normalized smoothed dataframe
after_df <- data.frame(
  mz = beforeafter_df$mz,
  intensity = smoothed_signal_norm
)

# After
(After_plot <- ggplot(after_df, aes(x = mz, y = intensity)) +
    geom_line(colour = "black", linewidth = 0.5) +
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


# First make sure smoothed_signal and mz have same length
length(smoothed_signal) == length(beforeafter_df$mz)  # Should return TRUE

# Create a new data frame for the smoothed line
after_df <- data.frame(
  mz = beforeafter_df$mz,
  intensity = smoothed_signal_norm
)

# Overlay the smoothed line on top of the raw spectrum
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

ggsave(overlay_plot, filename = "code/code/figures/presentations/ASASM430_after.png", bg = "transparent", 
       width = 35, height = 14, dpi = 300)  # Try adjusting these values
