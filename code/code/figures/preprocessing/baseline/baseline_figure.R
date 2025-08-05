# Matt Prill
# MSc Applied Data Science
# Baseline correction graph


# libraries ----
library(tidyverse)

# Data ----
load("data/merged_data/combined_spectra_filtered.RData")
load("")

# Example spectrum
A060 <- smoothed_spectra[["A060"]]
mass <- A060$mass
intensity <- A060$intensity


# The two SNIP functions ----

# SNIP baseline estimation function (normal log)
custom_snip <- function(intensity, iterations = 100) {
  log_intensity <- log1p(intensity)  # log(1 + x) to avoid log(0)
  baseline <- log_intensity  # Starting point for estimated baseline (log intensity)
  
  n <- length(intensity)  # initialised
  for (k in 1:iterations) {
    for (i in (k + 1):(n - k)) {  # for each intensity value excl boundaries
      avg <- (baseline[i - k] + baseline[i + k]) / 2  # calculate avg intensities for 1 point either side of i
      if (baseline[i] > avg) {  # if intensity is higher than avg
        baseline[i] <- avg  # its replaced with the average
      }
    }
  }
  expm1(baseline)  # inverse log1p to return to intensity scale
}



#  SNIP baseline estimation function (LLS)
custom_snip_LLS <- function(intensity, iterations = 100) {
  log_intensity <- log1p(log1p(sqrt(intensity + 1)))  # using LLS transform
  baseline <- log_intensity  # Starting point for estimated baseline (log intensity)
  
  n <- length(intensity)  # initialised
  for (k in 1:iterations) {
    for (i in (k + 1):(n - k)) {  # for each intensity value excl boundaries
      avg <- (baseline[i - k] + baseline[i + k]) / 2  # calculate avg intensities for 1 point either side of i
      if (baseline[i] > avg) {  # if intensity is higher than avg
        baseline[i] <- avg  # its replaced with the average
      }
    }
  }
  #  inverse of LLS operator for backtransform
  corrected_baseline <- ((exp(exp(baseline) - 1) - 1)^2) - 1
  return(corrected_baseline)
}




# Applying the baseline corrections ----
# Apply custom SNIP (Log)
baseline_est <- custom_snip(intensity, iterations = 100)
corrected <- intensity - baseline_est

# Apply custom SNIP (LLS)
baseline_est_LLS <- custom_snip_LLS(intensity, iterations = 100)
corrected_LLS <- intensity - baseline_est_LLS

# Visualising the differences ----

# 3 stacked plots
par(mfrow = c(3, 1), mar = c(4, 4, 3, 2))  # 3 rows, shared margins

# A: Raw spectrum
plot(mass, intensity, type = "l", col = "black", main = "A) Raw Spectrum",
     xlab = "m/z", ylab = "Intensity", xlim = c(min(mass), 8000))

# B: Baseline-corrected (log1p SNIP)
plot(mass, corrected, type = "l", col = "darkgreen", main = "B) SNIP-Corrected (log1p)",
     xlab = "m/z", ylab = "Corrected Intensity", xlim = c(min(mass), 8000))

# C: Baseline-corrected (LLS SNIP)
plot(mass, corrected_LLS, type = "l", col = "darkorange", main = "C) SNIP-Corrected (LLS)",
     xlab = "m/z", ylab = "Corrected Intensity", xlim = c(min(mass), 8000))


# Final Plot ----
library(ggplot2)
library(patchwork)

# Create dataframes for each plot
df_raw <- data.frame(mass = mass, intensity = intensity)
df_snip <- data.frame(mass = mass, intensity = corrected)
df_lls <- data.frame(mass = mass, intensity = corrected_LLS)

# A) Raw spectrum
p_raw <- ggplot(df_raw, aes(x = mass, y = intensity)) +
  geom_line(colour = "black") +
  labs(title = "A) Uncorrected Spectrum", x = "\nm/z", y = "Intensity\n") +
  xlim(min(mass), 8000) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(colour = "black", size = 40),
    axis.title.y = element_text(colour = "black", size = 40),
    axis.text.x = element_text(colour = "black", size = 30),
    axis.text.y = element_text(colour = "black", size = 30),
    panel.grid = element_blank(),
    plot.title = element_text(size = 36, face = "bold", hjust = 0)
  )

# B) SNIP-corrected (log1p)
p_snip <- ggplot(df_snip, aes(x = mass, y = intensity)) +
  geom_line(colour = "darkgreen") +
  labs(title = "B) SNIP-Corrected (log1p)", x = "\nm/z", y = "Corrected Intensity\n") +
  xlim(min(mass), 8000) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(colour = "black", size = 40),
    axis.title.y = element_text(colour = "black", size = 40),
    axis.text.x = element_text(colour = "black", size = 30),
    axis.text.y = element_text(colour = "black", size = 30),
    panel.grid = element_blank(),
    plot.title = element_text(size = 36, face = "bold", hjust = 0)
  )

# C) SNIP-corrected (LLS)
p_lls <- ggplot(df_lls, aes(x = mass, y = intensity)) +
  geom_line(colour = "darkorange") +
  labs(title = "C) SNIP-Corrected (LLS)", x = "\nm/z", y = "Corrected Intensity\n") +
  xlim(min(mass), 8000) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(colour = "black", size = 40),
    axis.title.y = element_text(colour = "black", size = 40),
    axis.text.x = element_text(colour = "black", size = 30),
    axis.text.y = element_text(colour = "black", size = 30),
    panel.grid = element_blank(),
    plot.title = element_text(size = 36, face = "bold", hjust = 0)
  )

# Combine into one plot
baseline_plot_lss_snip <- (p_raw / p_snip / p_lls) + plot_layout(guides = "collect")

# Save
ggsave(
  baseline_plot_lss_snip,
  filename = "code/code/figures/preprocessing/baseline/baseline_plot_lss_snip.png",
  bg = "white",
  width = 30, height = 20, dpi = 300, limitsize = FALSE
)



# TopHat function  ----
# Example spectrum
A060 <- combined_spectra_filtered[["A060"]]
mass <- A060$mass
intensity <- A060$intensity


# Top-hat functions 
rolling_min <- function(x, window) {
  n <- length(x)                      # initialise
  half_win <- floor(window / 2)       # half window size
  result <- numeric(n)                # initialise
  for (i in 1:n) {                    # local window boundaries
    left <- max(1, i - half_win)      # left of center
    right <- min(n, i + half_win)     # right of center
    result[i] <- min(x[left:right])   # minimum applied over window 
  }
  return(result)
}

# same as above but max is made the center
rolling_max <- function(x, window) {
  n <- length(x)
  half_win <- floor(window / 2)
  result <- numeric(n)
  for (i in 1:n) {
    left <- max(1, i - half_win)
    right <- min(n, i + half_win)
    result[i] <- max(x[left:right])
  }
  return(result)
}

# Choose a structuring element/window size
window_size <- 21  # must be odd to change the middle intensity

# Erosion
eroded <- rolling_min(intensity, window_size)

# Dilation of eroded signal
opened <- rolling_max(eroded, window_size)

# TopHat correction
tophat_corrected <- intensity - opened






# Top hat vs LLS-SNIP ----
#SNIP
baseline_est_LLS <- custom_snip_LLS(intensity, iterations = 100)  # snip (see function earlier)
corrected_LLS <- intensity - baseline_est_LLS

#Top-hat
window_size <- 21 
eroded <- rolling_min(intensity, window_size)
opened <- rolling_max(eroded, window_size)
tophat_corrected <- intensity - opened




# Visualising the differences ----
# A: Raw spectrum
plot(mass, intensity, type = "l", col = "black", main = "A) Raw Spectrum",
     xlab = "m/z", ylab = "Intensity", xlim = c(2000, 8000))

# B: Baseline-corrected (LLS SNIP)
plot(mass, corrected_LLS, type = "l", col = "blue", main = "B) SNIP-Corrected (LLS)",
     xlab = "m/z", ylab = "Corrected Intensity", xlim = c(2000, 8000))

# C: Baseline-corrected (Top-hat)
plot(mass, tophat_corrected, type = "l", col = "blue", main = "C) TopHat Baseline Corrected", 
    xlab = "m/z", ylab = "Corrected Intensity", xlim = c(2000, 8000))



# Final Plot

# Create dataframes
df_raw     <- data.frame(mass = mass, intensity = intensity)
df_lls     <- data.frame(mass = mass, intensity = corrected_LLS)
df_tophat  <- data.frame(mass = mass, intensity = tophat_corrected)

# A) Raw
p_raw <- ggplot(df_raw, aes(x = mass, y = intensity)) +
  geom_line(colour = "black") +
  labs(title = "A) Uncorrected", x = "\nm/z", y = "Intensity\n") +
  xlim(2000, 8000) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(colour = "black", size = 40),
    axis.title.y = element_text(colour = "black", size = 40),
    axis.text.x = element_text(colour = "black", size = 30),
    axis.text.y = element_text(colour = "black", size = 30),
    panel.grid = element_blank(),
    plot.title = element_text(size = 36, face = "bold", hjust = 0)
  )

# B) LLS SNIP
p_lls <- ggplot(df_lls, aes(x = mass, y = intensity)) +
  geom_line(colour = "darkorange") +
  labs(title = "B) SNIP", x = "\nm/z", y = "Corrected Intensity\n") +
  xlim(2000, 8000) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(colour = "black", size = 40),
    axis.title.y = element_text(colour = "black", size = 40),
    axis.text.x = element_text(colour = "black", size = 30),
    axis.text.y = element_text(colour = "black", size = 30),
    panel.grid = element_blank(),
    plot.title = element_text(size = 36, face = "bold", hjust = 0)
  )

# C) Top-hat
p_tophat <- ggplot(df_tophat, aes(x = mass, y = intensity)) +
  geom_line(colour = "darkgreen") +
  labs(title = "C) TopHat", x = "\nm/z", y = "Corrected Intensity\n") +
  xlim(2000, 8000) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(colour = "black", size = 40),
    axis.title.y = element_text(colour = "black", size = 40),
    axis.text.x = element_text(colour = "black", size = 30),
    axis.text.y = element_text(colour = "black", size = 30),
    panel.grid = element_blank(),
    plot.title = element_text(size = 36, face = "bold", hjust = 0)
  )

# Combine
baseline_comparison_plot <- (p_raw / p_lls / p_tophat) + plot_layout(guides = "collect")

# Save
ggsave(
  baseline_comparison_plot,
  filename = "code/code/figures/preprocessing/baseline/baseline_comparison_plot.png",
  bg = "white",
  width = 30, height = 20, dpi = 300, limitsize = FALSE
)





