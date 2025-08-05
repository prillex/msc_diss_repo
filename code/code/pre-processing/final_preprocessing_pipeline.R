# Matt Prill
# MSc Applied Data Science
# Final Pre-processing pipeline (excluding merging data)
# WARNING -  The code may take ~30 minutes to run

# Seed ----
set.seed(1234)  # For reproducibility 

# Libraries ----
library(tidyverse)  # all sorts
library(reticulate)  # python
library(patchwork)  # plot
library(beepr)  # ding!


# Data ----
# Just raw spectra (RData with mass (to charge) and intensity)
load("data/merged_data/combined_spectra_raw_filtered.RData")  # called combined_spectra_raw_filtered

# Interpolation ----
# Resampling onto a uniform m/z axis 

mass_ranges <- sapply(combined_spectra_raw_filtered, function(df) range(df$mass))
global_min <- max(mass_ranges[1, ])  # max of mins to avoid extrapolation
global_max <- min(mass_ranges[2, ])  # min of maxes

# Defining grid
n_points <- 21450  # max no. data points per isolate
mass_grid <- seq(from = global_min, to = global_max, length.out = n_points)
step_size <- (global_max - global_min) / (n_points - 1)
step_size  # ~0.84


# Interpolating
resample_spectrum <- function(df, mass_grid) {
  interpolated <- approx(x = df$mass, y = df$intensity, xout = mass_grid, rule = 2)$y  # only interpolate overlapping regions
  tibble(mass = mass_grid, intensity = interpolated)
}  # Function

# Apply to all spectra
resampled_spectra <- map(combined_spectra_raw_filtered, resample_spectrum, mass_grid = mass_grid)
names(resampled_spectra) <- names(combined_spectra_raw_filtered)  # preserve isolate names (unnecessary for my spectra)



# Modified Sinc Kernel Smoothing ----
# Creating the functions
generate_ms_kernel <- function(m, n = 6, alpha = 4) {  # m = window size, n = how fast the sinc function oscillates, alpha = Gaussian-based weighting decay rate
  x <- seq(-m, m) / (m + 1)  # scale window from -1 to 1 (in case I change m)
  w <- exp(-alpha * x^2) +  # Building the window - Gaussian centered at 0
    exp(-alpha * (x - 2)^2) +  # gaussian centered at +2
    exp(-alpha * (x + 2)^2) -  # gaussian centered at -2
    2 * exp(-alpha * 1^2)  # normalise offset (to make it a band-pass filter)
  sinc <- function(x) ifelse(x == 0, 1, sin(pi * x) / (pi * x))  # sinc function
  kernel <- w * sinc(n * x)  # multiply window by function (apply sinc modulation)
  a <- 0.00172; b <- 0.02437; c <- 1.64375  # constants for correction (emprical values taken from Schmid et al., 2022)
  correction <- a + b / (c - (m + 1))  # empirical sinusoidal correction (Schmid 2022)
  kernel <- kernel + correction * sin((2 + 1) * pi * x)  # add correction for ripple artifacts (cancels out frequency-domain ripples)
  kernel / sum(kernel)  # must normalise so a weight of >1 is not applied to the spectra
}

# To stretch out data at the edges (to avoid edge distortion during smoothing)
linear_extrapolate <- function(y, m) {
  left_lm <- lm(y[1:m] ~ seq_len(m))  # fit line to start of signal
  right_lm <- lm(y[(length(y) - m + 1):length(y)] ~ seq_len(m))  # fit line to end
  left_ext <- predict(left_lm, newdata = data.frame(seq_len = seq(-m, -1)))  # extrapolate left
  right_ext <- predict(right_lm, newdata = data.frame(seq_len = seq(length(y) + 1, length(y) + m)))  # extrapolate right
  c(left_ext, y, right_ext)  # combine padded signal
}

# Smoothing data with convolution function (centred via "filter")
ms_smoother <- function(y, m = 30, n = 6) {
  kernel <- generate_ms_kernel(m, n)  # modified sinc kernel function
  y_ext <- linear_extrapolate(y, m)   # linear extrapolate function
  conv <- stats::convolve(y_ext, rev(kernel), type = "filter")  # apply convolution
  conv  # Return the smoothed signal - already same length as input
}

# Applying to the spectra
smoothed_spectra <- map(resampled_spectra, function(spec) {
  intensity <- spec$intensity
  mz <- spec$mass
  m <- 30
  smoothed <- ms_smoother(intensity, m = m, n = 6)
  list(mass = mz, intensity = smoothed)  # only keep smoothed
})

# Preserve isolate names
names(smoothed_spectra) <- names(combined_spectra_raw_filtered)  # preserve isolate names (uncessary for my spectra)
#save(smoothed_spectra, file = "data/merged_data/smoothed_spectra.RData")


# Baseline Correction ----
# SNIP baseline correction function 
custom_snip_LLS <- function(intensity, iterations = 100) {
  intensity[intensity < 0] <- 0  # make any negative values 0
  log_intensity <- log1p(log1p(sqrt(intensity + 1)))  # LLS transform
  baseline <- log_intensity        # Set initial baseline
  n <- length(intensity)
  for (k in 1:iterations) {        # SNIP for the given number of iterations
    for (i in (k + 1):(n - k)) {   # For each intensity, skipping edges to avoid out-of-bounds
      avg <- (baseline[i - k] + baseline[i + k]) / 2  # calculate average intensity of adjacent neighbours
      if (baseline[i] > avg) {     # If intensity is bigger than avg of two either side,
        baseline[i] <- avg         # make it the avg of those two (clipped)
      }
    }
  }
  corrected_baseline <- ((exp(exp(baseline) - 1) - 1)^2) - 1  # inverse LLS
  return(corrected_baseline)
}

# Apply SNIP baseline correction, preserving names
baseline_corrected_spectra <- map(smoothed_spectra, function(spectrum) {
  intensity <- spectrum$intensity  # already smoothed spectra
  intensity[is.na(intensity)] <- 0  # turn any NAs to 0
  mz <- spectrum$mass
  
  estimated_baseline <- custom_snip_LLS(intensity)   # apply previous SNIP function to estimate baseline
  corrected_intensity <- intensity - estimated_baseline  # subtract estimated baseline from intensity
  
  list(  # back to same format as before
    mass = mz,
    intensity = corrected_intensity
  )
})

# Warning, all have marginally negative values now
# Preserve isolate names
names(baseline_corrected_spectra) <- names(combined_spectra_raw_filtered)  # preserve isolate names (unnecessary for my spectra)
# save(baseline_corrected_spectra, file = "data/merged_data/baseline_corrected_spectra.RData")





# Normalisation ----
# Load and inspect
load("data/merged_data/baseline_corrected_spectra.RData")
str(baseline_corrected_spectra)

#rescale function
safe_rescale <- function(x) {
  r <- range(x, na.rm = TRUE, finite = TRUE)
  if (r[1] == r[2]) {
    warning("Flat spectrum detected — all intensities are equal.")  # to say if there's a dodgy spectrum (max and min same)
    return(rep(0, length(x)))
  } else {
    return((x - r[1]) / (r[2] - r[1]))
  }
}


scaled_spectra <- lapply(baseline_corrected_spectra, function(spec) {
  list(
    mass = spec$mass,
    intensity = safe_rescale(spec$intensity)
  )
})


# save(scaled_spectra, file = "data/merged_data/scaled_spectra.RData")







# Peak Selection ----
load("data/merged_data/scaled_spectra.RData")

# python set up
use_virtualenv("~/.virtualenvs/r-reticulate", required = TRUE)
py_config()
py_run_string("globals().clear()")

# python peak picking 
py_run_string("
import numpy as np
from scipy.signal import find_peaks

def pick_peaks(mz, intensity):
    mz = np.array(mz)
    intensity = np.array(intensity)
    peaks, props = find_peaks(intensity, prominence = 0.02, height = 0.02)
    return mz[peaks].tolist(), intensity[peaks].tolist()
")

# run for each spectrum
peak_list <- lapply(scaled_spectra, function(spec) {
  mz <- spec$mass
  intensity <- spec$intensity
  py$mz <- mz
  py$intensity <- intensity
  py_run_string("mz_out, int_out = pick_peaks(mz, intensity)")
  list(mz = py$mz_out, intensity = py$int_out)
})

# the range of peaks counts per isolate is 33 (A125) to 894 (A060)


# save(peak_list, file = "data/merged_data/peak_list.RData")












# Peak Alignment ----
load("data/merged_data/peak_list.RData")  # peak list 
load("data/merged_data/scaled_spectra.RData")


# alignment function
align_peaks_custom <- function(peak_list, tol) {  # defining function, input = list, each containing isoalte's peaks
  all_mz <- sort(unlist(lapply(peak_list, `[[`, "mz")))  # collect mz peaks from across isolates, sort by ascending
  bins <- c(); i <- 1  # where average of mz values will be stored
  while(i <= length(all_mz)) {  # for every peak, 
    start <- all_mz[i]; # starting from the first mz
    j <- i   # j = index of current peak of interest
    while(j+1 <= length(all_mz)  # while in bounds i.e. up to penultimate peak
          && all_mz[j+1] - start <= tol)   # and the next peak is within tolerance value of current peak of interest
      j <- j + 1  # then move forward to the next peak
    bins <- c(bins, mean(all_mz[i:j])); # otherwise calculate mean of bin's peaks
    i <- j + 1  # then move forward to the next ungrouped peak (stops overlapping bins)
  }
  iso <- names(peak_list)  # prep matrix - isolate names
  mat <- matrix(NA_real_, nrow = length(iso), ncol=length(bins),  # rows = isolate names, column = bins (the calculated mean mz)
                dimnames = list(iso, sprintf("%.4f", bins)))  # Give to 4 decimal places
  for(sp in iso) {  # for each  isolate
    mzv <- peak_list[[sp]]$mz;   # extract m/z
    iv <- peak_list[[sp]]$intensity  # extract intensity 
    for(k in seq_along(mzv)) {  # for each peak in the isolate
      idx <- which.min(abs(bins - mzv[k]))  # find the closest bin  (min absolute distance between a given peak and the bins). assign to idx
      if(abs(bins[idx] - mzv[k]) <= tol) mat[sp, idx] <- iv[k]  # if the distance to the bin is within tol, assign the intensity to corresponding matrix cell, otherwise ignore
    }
  }
  list(bins=bins, mat=mat)
}

# run alignment at 5 da
res <- align_peaks_custom(peak_list, tol = 5)  # tol is the mz within which peaks are aligned same

# summarise
cat("bins:", length(res$bins), "\n")             # no. distinct features at this point
cat("total peaks:", sum(!is.na(res$mat)), "\n")  # no. peaks assingned to the bins above
print(rowSums(!is.na(res$mat)))

# assign for saving
bins <- res$bins
mat  <- res$mat

# save(bins, mat, file="data/merged_data/aligned_peaks.RData")
beep(8)








# Labelling (Start of Pipeline 2) ----
load("data/merged_data/aligned_peaks.RData")  # called mat
meta <- read.csv("data/merged_data/combined_metadata_filtered.csv")  # metadata

# Checking no isolates have been lost
length(meta$sampleID)  # same no. samples
nrow(mat)   # same no. samples

all(rownames(mat) == meta$sampleID)  # Same samples but in a different order

meta <- meta[match(rownames(mat), meta$sampleID), ]  #  Matches them up again
all(rownames(mat) == meta$sampleID)  # proof it's fixed

labels <- meta$Outcome30days  # Extract response variable






# Visualising the Pipeline ----
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

# panels
p1 <- ggplot() +
  geom_line(aes(x = raw$mass, y = raw$intensity), colour = "black") +
  labs(title = "Raw", x = NULL, y = "Intensity") +
  theme_minimal()

p2 <- ggplot() +
  geom_line(aes(x = smoothed$mass, y = smoothed$intensity), colour = "black") +
  labs(title = "Smoothed", x = NULL, y = "Intensity") +
  theme_minimal()

p3 <- ggplot() +
  geom_line(aes(x = scaled$mass, y = scaled$intensity), colour = "black") +
  labs(title = "Baseline Corrected & Scaled", x = NULL, y = "Intensity") +
  theme_minimal()

p4 <- ggplot() +
  geom_segment(
    aes(x = peaks$mz, xend = peaks$mz, y = 0, yend = peaks$intensity),
    colour = "purple", alpha = 0.7
  ) +
  geom_point(
    aes(x = peaks$mz, y = peaks$intensity),
    colour = "purple", size = 1.5
  ) +
  labs(title = "Selected Peaks", x = "m/z", y = "Intensity") +
  theme_minimal()
  
  p5 <- ggplot(df_filtered, aes(x = mz, y = intensity)) +
  geom_segment(aes(x = mz, xend = mz, y = 0, yend = intensity), 
               colour = "darkorange", alpha = 0.7
               ) +
  geom_point(colour = "darkorange", size = 1.5) +
  labs(title = "Filtered Peaks (Post-Frequency Threshold)", x = "m/z", y = "Intensity") +
    xlim(range(raw$mass)) +
    theme_minimal()

  

# combine all vertically
(p1 / p2 / p3 / p4/ p5) +
  plot_annotation(title = "Preprocessing Pipeline: Isolate A060")

  

  
  
  # Isolate ASARM100 data
  raw <- combined_spectra_raw_filtered[["ASARM100"]]
  smoothed <- smoothed_spectra[["ASARM100"]]
  corrected <- baseline_corrected_spectra[["ASARM100"]]
  scaled <- scaled_spectra[["ASARM100"]]
  peaks <- peak_list[["ASARM100"]]
  filtered_peaks  <- filtered_matrix["ASARM100", ]
  df_filtered <- data.frame(
    mz = as.numeric(names(filtered_peaks)),
    intensity = as.numeric(filtered_peaks)
  )
  df_filtered <- na.omit(df_filtered)  # Optional: remove NAs if present
  
  # Panels
  p1 <- ggplot() +
    geom_line(aes(x = raw$mass, y = raw$intensity), colour = "black") +
    labs(title = "Raw", x = NULL, y = "Intensity") +
    theme_minimal()
  
  p2 <- ggplot() +
    geom_line(aes(x = smoothed$mass, y = smoothed$intensity), colour = "black") +
    labs(title = "Smoothed", x = NULL, y = "Intensity") +
    theme_minimal()
  
  p3 <- ggplot() +
    geom_line(aes(x = scaled$mass, y = scaled$intensity), colour = "black") +
    labs(title = "Baseline Corrected & Scaled", x = NULL, y = "Intensity") +
    theme_minimal()
  
  p4 <- ggplot() +
    geom_segment(
      aes(x = peaks$mz, xend = peaks$mz, y = 0, yend = peaks$intensity),
      colour = "purple", alpha = 0.7
    ) +
    geom_point(
      aes(x = peaks$mz, y = peaks$intensity),
      colour = "purple", size = 1.5
    ) +
    labs(title = "Selected Peaks", x = "m/z", y = "Intensity") +
    xlim(range(raw$mass)) +
    theme_minimal()
  
  p5 <- ggplot(df_filtered, aes(x = mz, y = intensity)) +
    geom_segment(aes(x = mz, xend = mz, y = 0, yend = intensity), 
                 colour = "darkorange", alpha = 0.7) +
    geom_point(colour = "darkorange", size = 1.5) +
    labs(title = "Filtered Peaks (Post-Frequency Threshold)", x = "m/z", y = "Intensity") +
    xlim(range(raw$mass)) +
    theme_minimal()
  
  # Combine all 5 plots vertically
  (p1 / p2 / p3 / p4 / p5) +
    plot_annotation(title = "Preprocessing Pipeline: Isolate ASARM100")
  
  
  
  
  
  
  
  
  
