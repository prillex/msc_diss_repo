library(tidyverse)
library(reticulate)

# Setup Python
use_virtualenv("~/.virtualenvs/r-reticulate", required = TRUE)
py_run_string("globals().clear()")

# Extract spectra
specs <- list(
  ASARM110 = scaled_spectra[["ASARM110"]],
  A060     = scaled_spectra[["A060"]]
)

# Detect peaks and collect in a tibble
peak_df <- imap_dfr(specs, function(spec, name) {
  py$mz <- spec$mass
  py$intensity <- spec$intensity
  
  py_run_string("
import numpy as np
from scipy.signal import find_peaks

mz = np.array(mz)
intensity = np.array(intensity)

peaks, _ = find_peaks(intensity, prominence=0.02, height=0.02)
mz_fp = mz[peaks]
intensity_fp = intensity[peaks]
")
  
  tibble(
    iso       = name,
    mz        = py$mz_fp,
    intensity = py$intensity_fp
  )
})

# Also get full spectra for reference lines
spectra_df <- imap_dfr(specs, function(spec, name) {
  tibble(
    iso       = name,
    mz        = spec$mass,
    intensity = spec$intensity
  )
})

# Color mapping
peak_colors <- c(ASARM110 = "darkgreen", A060 = "darkorange")

# Plot using ggplot
peak_plot <- ggplot() +
  geom_line(data = spectra_df,
            aes(x = mz, y = intensity),
            color = "grey", linewidth = 0.5) +
  geom_point(data = peak_df,
             aes(x = mz, y = intensity, color = iso),
             shape = 17, size = 1.5) +
  facet_wrap(~iso, ncol = 1, scales = "free_y") +
  scale_color_manual(values = peak_colors) +
  labs(
    x = "m/z",
    y = "Relative Intensity\n"
  ) +
  theme_minimal(base_family = "Arial") +
  theme(
    strip.text       = element_text(size = 20, face = "bold", hjust = 0),
    axis.title       = element_text(size = 24),
    axis.text        = element_text(size = 20),
    panel.grid       = element_blank(),
    legend.position  = "none",
    panel.spacing.y  = unit(1.5, "lines"),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

# Save plot
ggsave("code/code/figures/preprocessing/peak_selection/peak_detection_plot.png",
       plot = peak_plot, width = 16, height = 8, dpi = 300, bg = "white")
