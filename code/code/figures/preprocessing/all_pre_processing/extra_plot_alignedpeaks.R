# Load aligned peak matrix and bins
load("data/merged_data/aligned_peaks.RData")  # contains: mat, bins

# Extract aligned peaks for A060
aligned_peaks <- mat["A060", ]
df_aligned <- data.frame(
  mz = bins,
  intensity = as.numeric(aligned_peaks)
)
df_aligned <- na.omit(df_aligned)

# Define consistent theme
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

# Define axis limits
xlims <- c(2000, 10000)
ylims <- c(0, 2000)

# Panel 6: Aligned Peaks (Post-Binning) with formatting
p6 <- ggplot(df_aligned, aes(x = mz, y = intensity)) +
  geom_segment(aes(xend = mz, y = 0, yend = intensity), 
               colour = "blue", alpha = 0.7) +
  geom_point(colour = "blue", size = 1.5) +
  labs(title = "Aligned Peaks (Post-Binning)", x = "m/z", y = "Absolute Intensity\n") +
  xlim(xlims) +
  ylim(ylims) +
  pretty_theme


# Just remove ylim()
p6 <- ggplot(df_aligned, aes(x = mz, y = intensity)) +
  geom_segment(aes(xend = mz, y = 0, yend = intensity), 
               colour = "blue", alpha = 0.7) +
  geom_point(colour = "blue", size = 1.5) +
  labs(title = "Aligned Peaks (Post-Binning)", x = "m/z", y = "Absolute Intensity\n") +
  xlim(xlims) +
  pretty_theme



ggsave(
  filename = "code/code/figures/preprocessing/all_pre_processing/aligned_peaks.png",
  plot = p6,
  width = 16, height = 4, dpi = 300, limitsize = FALSE
)
