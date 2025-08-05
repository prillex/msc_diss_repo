# load the scaled data

library(reticulate)

# 1. Use your virtualenv
use_virtualenv("~/.virtualenvs/r-reticulate", required = TRUE)
py_config()

# 2. Clear the Python environment
py_run_string("globals().clear()")

# 3. Prepare spectra
spec1 <- scaled_spectra[["ASARM110"]]
spec2 <- scaled_spectra[["A060"]]

# 4. Set up a 2-row plot
par(mfrow = c(2, 1), mar = c(4, 4, 2, 1))  # adjust margins if needed

for (i in 1:2) {
  spec <- if (i == 1) spec1 else spec2
  name <- if (i == 1) "ASARM110" else "A060"
  
  mz <- spec$mass
  intensity <- spec$intensity
  
  # Push to Python
  py$mz <- mz
  py$intensity <- intensity
  
  # Run peak detection
  py_run_string("
import numpy as np
from scipy.signal import find_peaks

mz = np.array(mz)
intensity = np.array(intensity)

peaks, properties = find_peaks(intensity, prominence=0.02, height=0.02)

mz_fp = mz[peaks]
intensity_fp = intensity[peaks]
  ")
  
  # Plot
  plot(mz, intensity, type = "l", col = "grey40",
       main = paste0("Prominence Peaks (", name, ")"),
       xlab = "m/z", ylab = "Intensity", xlim = c(2900, 3100))
  
  points(py$mz_fp, py$intensity_fp, col = "darkred", pch = 17, cex = 0.7)
}


# Plot full spectrum ----
library(reticulate)

# 1. Activate virtual environment
use_virtualenv("~/.virtualenvs/r-reticulate", required = TRUE)
py_config()

# 2. Clear Python environment
py_run_string("globals().clear()")

# 3. Extract spectra
spec1 <- scaled_spectra[["ASARM110"]]
spec2 <- scaled_spectra[["A060"]]

# 4. Set up stacked plotting
par(mfrow = c(2, 1), mar = c(4, 4, 2, 1))  # rows = 2, columns = 1

for (i in 1:2) {
  spec <- if (i == 1) spec1 else spec2
  name <- if (i == 1) "ASARM110" else "A060"
  
  mz <- spec$mass
  intensity <- spec$intensity
  
  # Push vectors to Python
  py$mz <- mz
  py$intensity <- intensity
  
  # Run Python peak detection
  py_run_string("
import numpy as np
from scipy.signal import find_peaks

mz = np.array(mz)
intensity = np.array(intensity)

peaks, properties = find_peaks(intensity, prominence=0.02, height=0.02)

mz_fp = mz[peaks]
intensity_fp = intensity[peaks]
")
  
  # Plot full spectrum
  plot(mz, intensity, type = "l", col = "grey40",
       main = paste0("Prominence Peaks (", name, ")"),
       xlab = "m/z", ylab = "Intensity")
  
  # Overlay detected peaks
  points(py$mz_fp, py$intensity_fp, col = "darkred", pch = 17, cex = 0.7)
}


















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








library(tidyverse)
library(reticulate)

# Activate Python env
use_virtualenv("~/.virtualenvs/r-reticulate", required = TRUE)
py_run_string("globals().clear()")

# Extract A060 spectrum
spec <- scaled_spectra[["A060"]]
mz <- spec$mass
intensity <- spec$intensity

# Restrict to m/z 2900–3100
mz_trimmed <- mz[mz >= 2900 & mz <= 3100]
int_trimmed <- intensity[mz >= 2900 & mz <= 3100]

# Define prominence thresholds
prominences <- c(0.01, 0.02, 0.05, 0.1)
colors <- c("darkgreen", "darkorange", "darkblue", "purple")
names(colors) <- paste0("Prominence = ", prominences)

# Prepare peaks from Python
peak_df_list <- lapply(prominences, function(p) {
  py$mz <- mz_trimmed
  py$intensity <- int_trimmed
  py$p <- p
  
  py_run_string("
import numpy as np
from scipy.signal import find_peaks

mz = np.array(mz)
intensity = np.array(intensity)
peaks, _ = find_peaks(intensity, prominence=p)
mz_fp = mz[peaks]
intensity_fp = intensity[peaks]
")
  
  tibble(
    mz = py$mz_fp,
    intensity = py$intensity_fp,
    prominence = paste0("Prominence = ", p)
  )
})

# Combine peaks
peak_df <- bind_rows(peak_df_list)

# Full spectrum data (trimmed range)
zoom_spec <- tibble(mz = mz_trimmed, intensity = int_trimmed)

# Plot
(prominence_effect <- ggplot() +
  geom_line(data = zoom_spec, aes(x = mz, y = intensity), color = "grey40") +
  geom_point(data = peak_df,
             aes(x = mz, y = intensity, color = prominence),
             size = 3) +
  facet_wrap(~ prominence, ncol = 2) +
  scale_color_manual(values = colors) +
  labs(
    title = "A060: Peak Detection at Different Prominence Thresholds",
    subtitle = "Zoomed-in view from 2900 to 3100 m/z",
    x = "m/z",
    y = "Relative Intensity\n"
  ) +
  theme_minimal(base_family = "Arial") +
  theme(
    strip.text       = element_text(face = "bold", size = 16),
    plot.title       = element_text(size = 20, face = "bold", hjust = 0),
    plot.subtitle    = element_text(size = 20, hjust = 0),
    axis.title       = element_text(size = 24),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text        = element_text(size = 20),
    legend.position  = "none"
  ))

#ggsave("code/code/figures/preprocessing/peak_selection/prominence_effect.png",
#       plot = prominence_effect, width = 20, height = 10, dpi = 300, bg = "white")



