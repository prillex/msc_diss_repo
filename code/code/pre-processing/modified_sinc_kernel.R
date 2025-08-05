# Matt Prill
# MSc Applied Data Science
# Dissertation
# Comparing MS Kernel parameters
# Note -  this does not contain the final sinc kernel code.

# Seed ----
seed(1234)  # For reproducibility 
# Libraries ----
library(tidyverse)
library(reticulate)

# Data ----
load("data/merged_data/combined_spectra_raw_filtered.RData")  # called combined_spectra_filtered



# Interpolation ----
# Resampling onto a uniform m/z axis 

mass_ranges <- sapply(combined_spectra_raw_filtered, function(df) range(df$mass))
global_min <- max(mass_ranges[1, ])  # max of mins to avoid extrapolation
global_max <- min(mass_ranges[2, ])  # min of maxes

# Defining grid
n_points <- 21450  # or fewer if needed
mass_grid <- seq(from = global_min, to = global_max, length.out = n_points)
step_size <- (global_max - global_min) / (n_points - 1)
step_size  # ~0.84



# Interpolating
resample_spectrum <- function(df, mass_grid) {
  interpolated <- approx(x = df$mass, y = df$intensity, xout = mass_grid, rule = 2)$y  # only interpolate overlapping regions
  tibble(mass = mass_grid, intensity = interpolated)
}

# Apply to all spectra
resampled_spectra <- map(combined_spectra_raw_filtered, resample_spectrum, mass_grid = mass_grid)
names(resampled_spectra) <- names(combined_spectra_raw_filtered)  # preserve IDs


# MS Smoothing Kernel Functions ---- 
# --- MS Kernel Functions ---
generate_ms_kernel <- function(m, n = 6, alpha = 4) {  # # m = window size, n = how fast the sinc function oscillates, alpha = Gaussian-based weighting decay rate
  x <- seq(-m, m) / (m + 1)    # normalise (scale) window from -1 to 1 (incase i change m)
  w <- exp(-alpha * x^2) +     # Building the window - Gaussian centered at 0
    exp(-alpha * (x - 2)^2) +  # gaussian centered at +2
    exp(-alpha * (x + 2)^2) -  # gaussian centered at -2
    2 * exp(-alpha * 1^2)      # normalisation offset
  sinc <- function(x) ifelse(x == 0, 1, sin(pi * x) / (pi * x))  # sinc function
  kernel <- w * sinc(n * x)  # multiply window by function (apply sinc modulation)
  a <- 0.00172; b <- 0.02437; c <- 1.64375  # constants for correction (emprical values taken from Schmid 2022)
  correction <- a + b / (c - (m + 1))  # empirical sinusoidal correction (Schmid 2022)
  kernel <- kernel + correction * sin((2 + 1) * pi * x)  # add correction for ripple artifacts (cancels out frequency-domain ripples)
  kernel / sum(kernel)  # must normalise so a weight of >1 is not applied to the spectra
}

# To stretch out data at the edges
linear_extrapolate <- function(y, m) {
  left_lm <- lm(y[1:m] ~ seq_len(m))  # fit line to start of signal
  right_lm <- lm(y[(length(y) - m + 1):length(y)] ~ seq_len(m))  # fit line to end
  left_ext <- predict(left_lm, newdata = data.frame(seq_len = seq(-m, -1)))  # extrapolate left
  right_ext <- predict(right_lm, newdata = data.frame(seq_len = seq(length(y) + 1, length(y) + m)))  # extrapolate right
  c(left_ext, y, right_ext)  # combine padded signal
}

# Smoothing data with convolution
ms_smoother <- function(y, m = 30, n = 6) {  # Larger m = more smoothing, risk of blurring sharp peaks. Smaller m = keeps sharp peaks, but might not reduce noise well.
  kernel <- generate_ms_kernel(m, n)  # create kernel
  y_ext <- linear_extrapolate(y, m)  # Extend the input vector so convolution won’t break at the edges.
  conv <- stats::convolve(y_ext, rev(kernel), type = "filter")  # apply convolution
  conv[(m + 1):(length(conv) - m)]  # remove padding
}

# Re-pad end so that smooth spectra has same length as before
pad_ms <- function(smoothed, m, original_len) {
  pad_left <- rep(NA, m)  # pad start
  pad_right <- rep(NA, original_len - length(smoothed) - m)  # pad end
  c(pad_left, smoothed, pad_right)  # combine
}

# --- Extract and Smooth A060 ---
mz <- resampled_spectra[["A060"]]$mass
intensity <- resampled_spectra[["A060"]]$intensity
smoothed <- ms_smoother(intensity, m = 30, n = 6)

# Trim m/z to match smoothed length (convolution shortens by m at both ends)
m <- 30
mz_trimmed <- mz[(m + 1):(length(mz) - m)]

# Overlay plot: original and smoothed
par(mfrow = c(1, 1))
plot(mz, intensity, type = "l", col = "grey60", lwd = 1.2,
     main = "A060 – Raw vs Smoothed Spectrum", xlab = "m/z", ylab = "Intensity",
     xlim = c(2000,8000))
lines(mz_trimmed, smoothed, col = "blue", lwd = 2)
legend("topright", legend = c("Raw", "Smoothed (MS Kernel)"),
       col = c("grey60", "blue"), lwd = c(1.2, 2), bty = "n")


# Plotting the kernel itself ----
# use your existing kernel generator
generate_ms_kernel <- function(m, n = 6, alpha = 4) {
  x <- seq(-m, m) / (m + 1)
  w <- exp(-alpha * x^2) +
    exp(-alpha * (x - 2)^2) +
    exp(-alpha * (x + 2)^2) -
    2 * exp(-alpha * 1^2)
  sinc <- function(x) ifelse(x == 0, 1, sin(pi * x) / (pi * x))
  kernel <- w * sinc(n * x)
  a <- 0.00172; b <- 0.02437; c <- 1.64375
  correction <- a + b / (c - (m + 1))
  kernel <- kernel + correction * sin((2 + 1) * pi * x)
  kernel / sum(kernel)
}

# parameters
m <- 30
n <- 6
alpha <- 4

# generate x and kernel
x <- seq(-m, m) / (m + 1)
kernel <- generate_ms_kernel(m, n, alpha)

# plot
plot(x, kernel, type = "l", col = "blue", lwd = 2,
     main = "Modified Sinc Kernel (MS)",
     xlab = "Normalised window (−1 to 1)",
     ylab = "Kernel value")
abline(h = 0, lty = 2, col = "grey60")











# Multiple plots of diff kernels ----

library(ggplot2)
library(dplyr)

# kernel generator
generate_ms_kernel <- function(m, n = 6, alpha = 4) {
  x <- seq(-m, m) / (m + 1)
  w <- exp(-alpha * x^2) +
    exp(-alpha * (x - 2)^2) +
    exp(-alpha * (x + 2)^2) -
    2 * exp(-alpha * 1^2)
  sinc <- function(x) ifelse(x == 0, 1, sin(pi * x) / (pi * x))
  kernel <- w * sinc(n * x)
  a <- 0.00172; b <- 0.02437; c <- 1.64375
  correction <- a + b / (c - (m + 1))
  kernel <- kernel + correction * sin((2 + 1) * pi * x)
  kernel / sum(kernel)
}

# settings for 4 kernels
settings <- list(
  list(name = "A: Standard (α = 4, n = 6)", alpha = 4, n = 6),
  list(name = "B: High α (α = 10, n = 6)", alpha = 50, n = 6),
  list(name = "C: High n (α = 4, n = 12)", alpha = 4, n = 12),
  list(name = "D: Distorted (α = 0.5, n = 1)", alpha = 0.5, n = 1)
)

# generate kernel data
m <- 30
plot_data <- lapply(settings, function(s) {
  x_vals <- seq(-m, m) / (m + 1)
  y_vals <- generate_ms_kernel(m, n = s$n, alpha = s$alpha)
  data.frame(x = x_vals, y = y_vals, setting = s$name)
}) %>% bind_rows()

# plot with ggplot
ggplot(plot_data, aes(x, y)) +
  geom_line(colour = "blue", size = 1) +
  facet_wrap(~ setting, ncol = 2) +
  theme_minimal(base_size = 12) +
  labs(title = "Effect of alpha and frequency on MS Kernel Morphology",
       x = "Normalised window (−1 to 1)",
       y = "Kernel value") +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60")

