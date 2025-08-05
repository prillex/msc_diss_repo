# Matt Prill
# MSc Dissertation
# Pre-processing mass spectra
# Up to now, the mass spectra have been merged and all redundant isolates removed

set.seed(1234)

# Libraries ----
library(reticulate)  # Python
library(MALDIquant)  # Baseline Correction etc.


# Python environment configuration ----
use_virtualenv("~/.virtualenvs/r-reticulate", required = TRUE)
py_config()  # Double checks pyhton path


# The data ----
# Cork/cc22 spectra
load("data/merged_data/combined_spectra_raw_filtered.RData")  # called combined_spectra_filtered
str(combined_spectra_raw_filtered)



# Apply the SNIP function to all specrtra
# ensure data is combined_spectra_raw_filtered








# (MALDIQUANT demo) & plot example ----
# Baseline correction (Testing and comparing different methods)
A001 <- combined_spectra_filtered[["A001"]]  # Practice
A001_spectrum <- createMassSpectrum(mass = A001$mass, intensity = A001$intensity)

## Test different baseline estimation methods
# SNIP
bSnip <- estimateBaseline(A001_spectrum, method = "SNIP", iterations = 100)  # ESTIMATED BASELINE
A001_corrected <- removeBaseline(A001_spectrum, method = "SNIP", iterations = 100)  # BASELINE REMOVED

# Plot
plot(A001_spectrum, main = "SNIP Baseline Correction", xlab = "m/z", ylab = "Intensity")
lines(bSnip, col = "red", lwd = 2)                      # The estimated baseline
lines(A001_corrected, col = "blue", lwd = 2)            # The corrected spectrum
legend("topright",
       legend = c("Original", "Estimated Baseline", "Corrected"),
       col = c("black", "red", "blue"),
       lty = 1, lwd = 2)



### TopHat ###
bTopHat <- estimateBaseline(A001_spectrum, method = "TopHat", halfWindowSize = 20)
A001_corrected_tophat <- removeBaseline(A001_spectrum, method = "TopHat", halfWindowSize = 20)

plot(A001_spectrum, main = "TopHat Baseline Correction", xlab = "m/z", ylab = "Intensity")
lines(bTopHat, col = "red", lwd = 2)
lines(A001_corrected_tophat, col = "blue", lwd = 2)
legend("topright", legend = c("Original", "Estimated Baseline", "Corrected"),
       col = c("black", "red", "blue"), lty = 1, lwd = 2)



### ConvexHull ###
bConvexHull <- estimateBaseline(A001_spectrum, method = "ConvexHull")
A001_corrected_convex <- removeBaseline(A001_spectrum, method = "ConvexHull")

plot(A001_spectrum, main = "ConvexHull Baseline Correction", xlab = "m/z", ylab = "Intensity")
lines(bConvexHull, col = "red", lwd = 2)
lines(A001_corrected_convex, col = "blue", lwd = 2)
legend("topright", legend = c("Original", "Estimated Baseline", "Corrected"),
       col = c("black", "red", "blue"), lty = 1, lwd = 2)



### Median ###
bMedian <- estimateBaseline(A001_spectrum, method = "median")
A001_corrected_median <- removeBaseline(A001_spectrum, method = "median")

plot(A001_spectrum, main = "Median Baseline Correction", xlab = "m/z", ylab = "Intensity")
lines(bMedian, col = "red", lwd = 2)
lines(A001_corrected_median, col = "blue", lwd = 2)
legend("topright", legend = c("Original", "Estimated Baseline", "Corrected"),
       col = c("black", "red", "blue"), lty = 1, lwd = 2)



# Comparison Plot
corrected_snip <- removeBaseline(A001_spectrum, method = "SNIP", iterations = 100)
corrected_tophat <- removeBaseline(A001_spectrum, method = "TopHat", halfWindowSize = 20)
corrected_convex <- removeBaseline(A001_spectrum, method = "ConvexHull")
corrected_median <- removeBaseline(A001_spectrum, method = "median")

# 2×2 grid (overlays)
par(mfrow = c(2, 2))  # 2 rows, 2 columns


plot(A001_spectrum, main = "SNIP", xlab = "m/z", ylab = "Intensity", xlim = c(2000, 10000))
lines(corrected_snip, col = "blue", lwd = 2)

plot(A001_spectrum, main = "TopHat", xlab = "m/z", ylab = "Intensity", xlim = c(2000, 10000))
lines(corrected_tophat, col = "blue", lwd = 2)

plot(A001_spectrum, main = "ConvexHull", xlab = "m/z", ylab = "Intensity", xlim = c(2000, 10000))
lines(corrected_convex, col = "blue", lwd = 2)

plot(A001_spectrum, main = "Median", xlab = "m/z", ylab = "Intensity", xlim = c(2000, 10000))
lines(corrected_median, col = "blue", lwd = 2)


# 2×2 grid (3 methods)
par(mfrow = c(2, 2))  # 2 rows, 2 columns

# Raw
plot(A001_spectrum, main = "Raw Spectrum", xlab = "m/z", ylab = "Intensity", col = "black", xlim = c(2000, 8000))

# 2. SNIP
plot(corrected_snip, main = "SNIP Baseline Corrected", xlab = "m/z", ylab = "Intensity", col = "blue", xlim = c(2000, 8000))

# 3. TopHat
plot(corrected_tophat, main = "TopHat Baseline Corrected", xlab = "m/z", ylab = "Intensity", col = "blue", xlim = c(2000, 8000))

# 4. ConvexHull
plot(corrected_convex, main = "ConvexHull Baseline Corrected", xlab = "m/z", ylab = "Intensity", col = "blue", xlim = c(2000, 8000))

# Reset plotting layout
par(mfrow = c(1, 1))






# The same but for a more nosisy spectra
A060 <- combined_spectra_filtered[["A060"]]  # Practice
A060_spectrum <- createMassSpectrum(mass = A060$mass, intensity = A060$intensity)

# Comparison Plot
  corrected_snip_A060 <- removeBaseline(A060_spectrum, method = "SNIP", iterations = 100)
corrected_tophat_A060 <- removeBaseline(A060_spectrum, method = "TopHat", halfWindowSize = 20)
corrected_convex_A060 <- removeBaseline(A060_spectrum, method = "ConvexHull")
corrected_median_A060 <- removeBaseline(A060_spectrum, method = "median")

# 2×2 grid (3 methods)
par(mfrow = c(2, 2))  # 2 rows, 2 columns

# Raw
plot(A060_spectrum, main = "Raw Spectrum", xlab = "m/z", ylab = "Intensity", col = "black", xlim = c(2000, 8000))

# 2. SNIP
plot(corrected_snip_A060, main = "SNIP Baseline Corrected", xlab = "m/z", ylab = "Intensity", col = "blue", xlim = c(2000, 8000))

# 3. TopHat
plot(corrected_tophat_A060, main = "TopHat Baseline Corrected", xlab = "m/z", ylab = "Intensity", col = "blue", xlim = c(2000, 8000))

# 4. ConvexHull
plot(corrected_convex_A060, main = "ConvexHull Baseline Corrected", xlab = "m/z", ylab = "Intensity", col = "blue", xlim = c(2000, 8000))

# Reset plotting layout
par(mfrow = c(1, 1))



# THINK - DO I NEED TO MAKE SURE I DONT DO THIS TO THE ALREADY CORRECTED SPECTRA?


# My own SNIP functions ----

# own snip baseline estimation function (normal log)
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
  return(expm1(baseline))  # inverse log1p to return to intensity scale
}

# Load a single spectrum for demonstration
A060 <- combined_spectra_filtered[["A060"]]
mass <- A060$mass
intensity <- A060$intensity

# Apply custom SNIP
baseline_est <- custom_snip(intensity, iterations = 100)
corrected <- intensity - baseline_est

# Plot original, baseline, and corrected
par(mfrow = c(1, 1))
plot(mass, intensity, type = "l", col = "black", main = "Custom SNIP Baseline Correction", xlab = "m/z", ylab = "Intensity")
lines(mass, baseline_est, col = "red", lwd = 2)
lines(mass, corrected, col = "blue", lwd = 2)
legend("topright", legend = c("Original", "Estimated Baseline", "Corrected"), col = c("black", "red", "blue"), lty = 1, lwd = 2)




# own snip baseline estimation function (LLS)
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


baseline_est_LLS <- custom_snip_LLS(intensity, iterations = 100)
corrected_LLS <- intensity - baseline_est_LLS

# Plot original, baseline, and corrected
par(mfrow = c(1, 1))
plot(mass, intensity, type = "l", col = "black", main = "Custom SNIP Baseline Correction", xlab = "m/z", ylab = "Intensity")
lines(mass, baseline_est_LLS, col = "red", lwd = 2)
lines(mass, corrected_LLS, col = "blue", lwd = 2)
legend("topright", legend = c("Original", "Estimated Baseline", "Corrected"), col = c("black", "red", "blue"), lty = 1, lwd = 2)







# TopHat  ----
# Load a single spectrum for demonstration
A060 <- combined_spectra_filtered[["A060"]]
mass <- A060$mass
intensity <- A060$intensity


# Functions 
rolling_min <- function(x, window) {
  n <- length(x)  # initialise
  half_win <- floor(window / 2)  # half window size
  result <- numeric(n)  # initialise
  for (i in 1:n) {  # # local window boundaries
    left <- max(1, i - half_win)  # left of center
    right <- min(n, i + half_win)  # right of center
    result[i] <- min(x[left:right])  # minimum applied over window 
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
window_size <- 21  # must be odd

# Step 1: Erosion
eroded <- rolling_min(intensity, window_size)

# Step 2: Dilation of eroded signal
opened <- rolling_max(eroded, window_size)

# Step 3: TopHat correction
tophat_corrected <- intensity - opened


plot(mass, tophat_corrected, type = "l", main = "TopHat Baseline Corrected", ylab = "Corrected Intensity", xlab = "m/z")














# DWT/UDWT smoothing ----
load("data/merged_data/combined_spectra_raw_filtered.RData")  # called combined_spectra_filtered

# Library
library(waveslim)

# Check Python and import libraries
py_config()

# function
py_run_string("
import numpy as np
import pywt

def wavelet_smooth(signal, wavelet='db4', level=4, threshold_scale=1.0, mode='soft'):
    level = int(level)
    signal = np.asarray(signal)
    coeffs = pywt.wavedec(signal, wavelet, level=level)
    
    # Estimate noise level
    sigma = np.median(np.abs(coeffs[-1])) / 0.6745
    uthresh = threshold_scale * sigma * np.sqrt(2 * np.log(len(signal)))
    
    # Thresholding
    for i in range(1, len(coeffs)):
        coeffs[i] = pywt.threshold(coeffs[i], uthresh, mode=mode)
    
    smoothed = pywt.waverec(coeffs, wavelet)
    
    # Trim to original length
    return smoothed[:len(signal)].tolist()
")



# Extract raw spectrum
dfraw <- combined_spectra_filtered[["A060"]]
mass <- dfraw$mass
intensity <- dfraw$intensity

# Convert intensity to Python and run smoothing
smoothed <- py$wavelet_smooth(intensity, wavelet = "db4", level = 4, threshold_scale = 1.5, mode = "soft")



# Plot both
par(mfrow = c(2, 1))
plot(mass, intensity, type = "l", col = "black", main = "A) Raw Spectrum",
     xlab = "m/z", ylab = "Intensity", xlim = c(min(mass), 10000))
plot(mass, smoothed, type = "l", col = "red", main = "B) Python Wavelet Smoothed",
     xlab = "m/z", ylab = "Smoothed Intensity", xlim = c(min(mass), 10000))
par(mfrow = c(1, 1))



# UDWT ----
load("data/merged_data/combined_spectra_raw_filtered.RData")  # called combined_spectra_filtered

py_run_string("
import numpy as np
import pywt

def wavelet_smooth_udwt(signal, wavelet='db4', level=4, threshold_scale=1.0, mode='soft'):
    signal = np.asarray(signal)
    original_len = len(signal)

    # Pad to even length
    if original_len % 2 != 0:
        signal = np.append(signal, 0)

    # Get maximum allowed level
    max_level = pywt.swt_max_level(len(signal))
    level = min(level, max_level)

    if level < 1:
        raise ValueError(f'Requested level too high. Max SWT level for signal is {max_level}.')

    # SWT (Undecimated DWT)
    coeffs = pywt.swt(signal, wavelet, level=level)

    # Estimate noise
    sigma = np.median(np.abs(coeffs[-1][1])) / 0.6745
    threshold = threshold_scale * sigma * np.sqrt(2 * np.log(original_len))

    # Threshold detail coefficients
    coeffs_thresh = []
    for cA, cD in coeffs:
        cD_thresh = pywt.threshold(cD, threshold, mode=mode)
        coeffs_thresh.append((cA, cD_thresh))

    # Inverse SWT
    smoothed = pywt.iswt(coeffs_thresh, wavelet)

    return smoothed[:original_len].tolist()
")

dfraw <- combined_spectra_raw_filtered[["A060"]]
mass <- dfraw$mass
intensity <- dfraw$intensity

smoothed <- py$wavelet_smooth_udwt(intensity,
                                   wavelet = "haar",  # good balance between smoothness and resolution
                                   level = 4,  # Level 4 decomposition  chosen to balance noise suppression with signal preservation 
                                   threshold_scale = 1000,  # supress high freq noise but preserve peaks
                                   mode = "soft")  # ensure a smooth signal profile


# Plot
# Choose a zoom window
mz_min <- 3000
mz_max <- 3050

# Get indices within this m/z range
idx <- which(mass >= mz_min & mass <= mz_max)

# Plot the zoomed raw and smoothed spectra
par(mfrow = c(2, 1))  # Stack vertically

plot(mass[idx], intensity[idx], type = "l", col = "black",
     main = paste("A) Raw Spectrum (", mz_min, "-", mz_max, " m/z)", sep = ""),
     xlab = "m/z", ylab = "Intensity")

plot(mass[idx], smoothed[idx], type = "l", col = "red",
     main = paste("B) UDWT Smoothed (Python) (", mz_min, "-", mz_max, " m/z)", sep = ""),
     xlab = "m/z", ylab = "Smoothed Intensity")

par(mfrow = c(1, 1))  # Reset layout




# UDWT 2 ----

# UDWT Smoothing ----
library(reticulate)

# defined python udwt smoothing function 
py_run_string("
import numpy as np
import pywt

def wavelet_smooth_udwt(signal, wavelet='db4', level=4, threshold_scale=1.0, mode='soft'):
    signal = np.asarray(signal)
    original_len = len(signal)

    # pad to even length (required for swt)
    if original_len % 2 != 0:
        signal = np.append(signal, 0)

    # get maximum allowed decomposition level for swt
    max_level = pywt.swt_max_level(len(signal))
    level = min(level, max_level)  # ensure level does not exceed max allowed

    if level < 1:
        raise ValueError(f'requested level too high. max swt level for signal is {max_level}.')

    # perform stationary wavelet transform (udwt)
    coeffs = pywt.swt(signal, wavelet, level=level)

    # estimate noise from the finest detail coefficients
    sigma = np.median(np.abs(coeffs[-1][1])) / 0.6745

    # calculate threshold for soft/hard thresholding
    threshold = threshold_scale * sigma * np.sqrt(2 * np.log(original_len))

    # threshold detail coefficients to suppress noise
    coeffs_thresh = []
    for cA, cD in coeffs:
        cD_thresh = pywt.threshold(cD, threshold, mode=mode)
        coeffs_thresh.append((cA, cD_thresh))

    # reconstruct smoothed signal via inverse swt
    smoothed = pywt.iswt(coeffs_thresh, wavelet)

    return smoothed[:original_len].tolist()
")

# create an empty list to store smoothed spectra for all samples
smoothed_list <- list()

# loop over each data frame in the original list
for (name in names(combined_spectra_raw_filtered)) {
  df <- combined_spectra_raw_filtered[[name]]
  
  # extract intensity vector to smooth
  intensity <- df$intensity
  
  # apply udwt smoothing using the python function
  smoothed_signal <- py$wavelet_smooth_udwt(
    intensity,
    wavelet = "db4",           # good balance between smoothness and resolution
    level = 4,                 # level 4 decomposition chosen to balance noise suppression with signal preservation
    threshold_scale = 2.0,     # suppress high freq noise but preserve peaks
    mode = "soft"              # ensure a smooth signal profile without artifacts
  )
  
  # store result as data frame with original mass and smoothed intensity
  smoothed_list[[name]] <- data.frame(
    mass = df$mass,
    intensity_smoothed = unlist(smoothed_signal)
  )
}

# smoothed_list now contains smoothed spectra for all samples in the same structure as input
str(smoothed_list)










# Savitzky–Golay with windowed weights (SWG) ----

library(dplyr)
library(purrr)

# Hann-square window
hann_square_weights <- function(m) {
  x <- -m:m
  cos(pi * x / (2 * m))^4
}

# SGW smoother (for one spectrum)
sgw_smoother <- function(y, m = 25, degree = 4, window_func = hann_square_weights) {
  x <- -m:m
  weights <- window_func(m)
  A <- outer(x, 0:degree, `^`)
  W <- diag(weights)
  
  # Solve weighted least squares for convolution kernel
  coeffs <- solve(t(A) %*% W %*% A) %*% t(A) %*% W
  kernel <- coeffs[1, ]  # value at center (x=0)
  
  # Apply convolution filter (use sides=2 for centered)
  stats::filter(y, rev(kernel), sides = 2)
}

# Apply to resampled spectra
resampled_spectra_smooth <- map(resampled_spectra, ~ {
  y_smooth <- sgw_smoother(.x$intensity, m = 25, degree = 4)
  .x %>% mutate(intensity_smooth = as.numeric(y_smooth))
})



#Extract and smooth A060
iso_A060 <- resampled_spectra[["A060"]]
intensity_smooth <- sgw_smoother(iso_A060$intensity, m = 25, degree = 4)

# Combine into a single tibble
plot_df <- iso_A060 %>%
  mutate(
    smoothed = as.numeric(intensity_smooth),
    raw = intensity
  )

# Plot raw and smoothed overlaid
ggplot(plot_df, aes(x = mass)) +
  geom_line(aes(y = raw), color = "grey60", alpha = 0.6, size = 0.5) +
  geom_line(aes(y = smoothed), color = "steelblue", size = 0.6) +
  labs(
    title = "SGW Smoothing on Isolate A060",
    x = "m/z",
    y = "Intensity",
    caption = "Raw (grey) vs SGW-smoothed (blue)"
  ) +
  theme_minimal()







# Modified Sinc Kernel (MS) ----

library(dplyr)
library(ggplot2)

# --- MS Smoother Components ---

# Modified sinc kernel generator
generate_ms_kernel <- function(m, n = 6, alpha = 4) {
  x <- seq(-m, m) / (m + 1)
  
  # Gaussian-like window with zero-edge enforcement
  w <- exp(-alpha * x^2) +
    exp(-alpha * (x - 2)^2) +
    exp(-alpha * (x + 2)^2) -
    2 * exp(-alpha * 1^2)
  
  # Sinc function
  sinc <- function(x) ifelse(x == 0, 1, sin(pi * x) / (pi * x))
  kernel <- w * sinc(n * x)
  
  # Optional correction term from paper (Table 1)
  a <- 0.00172
  b <- 0.02437
  c <- 1.64375
  correction <- a + b / (c - (m + 1))
  kernel <- kernel + correction * sin((2 + 1) * pi * x)
  
  # Normalize to preserve peak height
  kernel / sum(kernel)
}

# Linear extrapolation at edges
linear_extrapolate <- function(y, m) {
  left_lm <- lm(y[1:m] ~ seq_len(m))
  right_lm <- lm(y[(length(y) - m + 1):length(y)] ~ seq_len(m))
  left_ext <- predict(left_lm, newdata = data.frame(seq_len = seq(-m, -1)))
  right_ext <- predict(right_lm, newdata = data.frame(seq_len = seq(length(y) + 1, length(y) + m)))
  c(left_ext, y, right_ext)
}

# MS smoother
ms_smoother <- function(y, m = 30, n = 6) {
  kernel <- generate_ms_kernel(m, n)
  y_ext <- linear_extrapolate(y, m)
  conv <- stats::convolve(y_ext, rev(kernel), type = "filter")
  conv[(m + 1):(length(conv) - m)]  # Trim extrapolated edges
}

# Pad smoothed vector to match original length
pad_ms <- function(smoothed, m, original_len) {
  pad_left <- rep(NA, m)
  pad_right <- rep(NA, original_len - length(smoothed) - m)
  c(pad_left, smoothed, pad_right)
}

# --- Apply to A060 ---

# Get isolate A060
iso_A060 <- resampled_spectra[["A060"]]

# Smooth and pad
m_val <- 30
intensity_ms <- ms_smoother(iso_A060$intensity, m = m_val, n = 6)
intensity_ms_padded <- pad_ms(intensity_ms, m = m_val, original_len = nrow(iso_A060))

# Add smoothed column to tibble
plot_df <- iso_A060 %>%
  mutate(ms_smoothed = intensity_ms_padded)

# --- Plot ---
ggplot(plot_df, aes(x = mass)) +
  geom_line(aes(y = intensity), color = "grey60", alpha = 0.5, size = 0.4) +
  geom_line(aes(y = ms_smoothed), color = "darkorange", size = 0.6) +
  labs(
    title = "Modified Sinc Kernel Smoothing — Isolate A060",
    x = "m/z",
    y = "Intensity",
    caption = "Raw (grey) vs MS-smoothed (orange)"
  ) +
  theme_minimal()


























# Transforms (MALDI) ----
A001 <- combined_spectra_filtered[["A001"]]  # Practice

## sqrt transform (for variance stabilisation)
s2 <- transformIntensity(A001, method="sqrt")

# 21 point Savitzky-Golay-Filter for smoothing spectra
s3 <- smoothIntensity(s2, method="SavitzkyGolay", halfWindowSize=10)

# remove baseline
s4 <- removeBaseline(s3, method="SNIP", iterations=100)



# Peak detection
p <- detectPeaks(s4, method="MAD", halfWindowSize=20, SNR=2)
par(mfrow=c(2,3))

xlim <- range(mass(A001)) # use same xlim on all plots for better comparison
plot(A001, main="1: raw", sub="", xlim=xlim)

plot(s2, main="2: variance stabilization", sub="", xlim=xlim)

plot(s3, main="3: smoothing", sub="", xlim=xlim)

plot(s4, main="4: baseline correction", sub="", xlim=xlim)

plot(s4, main="5: peak detection", sub="", xlim=xlim)

points(p)


# label top 20 peaks
top20 <- intensity(p) %in% sort(intensity(p), decreasing=TRUE)[1:20]

labelPeaks(p, index=top20, underline=TRUE)

plot(p, main="6: peak plot", sub="", xlim=xlim)

labelPeaks(p, index=top20, underline=TRUE)

par(mfrow=c(1,1))


# Warping
## use only 4 spectra
spectra <- combined_spectra_filtered[seq(1, 16, by=4)]

## some preprocessing
## sqrt transform (for variance stabilization)
spectra <- transformIntensity(spectra, method="sqrt")

## 21 point Savitzky-Golay-Filter for smoothing spectra
spectra <- smoothIntensity(spectra, method="SavitzkyGolay", halfWindowSize=10)

## remove baseline
spectra <- removeBaseline(spectra, method="SNIP", iterations=100)

## calibrate intensity values by "Total Ion Current"
spectra <- calibrateIntensity(spectra, method="TIC")

## run peak detection
peaks <- detectPeaks(spectra, method="MAD", halfWindowSize=20, SNR=2)

## warping
par(mfrow=c(2, 2))

warpingFunctions <- determineWarpingFunctions(peaks, tolerance=0.001,
                                              plot=TRUE, plotInteractive=TRUE)

## warp peaks
warpedPeaks <- warpMassPeaks(peaks, warpingFunctions)

## compare some regions in a plot
par(mfrow=c(2, 2))

## helper function to avoid double coding
plotSpectra <- function(unwarped, warped, range) {
       plot(unwarped[[1]], main=paste0("unwarped spectra (mass ",
             paste0(range, collapse=":"), " Da)"),
                     xlim=range, ylim=c(0, 2e-3), type="n")
  color <- rainbow(length(unwarped))
  for (i in seq(along=unwarped)) {
              lines(unwarped[[i]], col=color[i])
           }

   plot(unwarped[[1]], main=paste0("warped spectra (mass ",
   paste0(range, collapse=":"), " Da)"),
  xlim=range, ylim=c(0, 2e-3), type="n")
  for (i in seq(along=warped)) {
   lines(warped[[i]], col=color[i])
  }
}

plotSpectra(spectra, warpedSpectra, c(4180, 4240))
plotSpectra(spectra, warpedSpectra, c(9200, 9400))
par(mfrow=c(1, 1))

# Calibrate
spectra <- calibrateIntensity(spectra, method="TIC")

# Alignment
spectra <- alignSpectra(spectra,    
                        halfWindowSize=20, SNR=2,
                        tolerance=0.002, warpingMethod="lowess")

# (After alignment peak positions (mass) are similar but not identical. Binning
# is needed to make similar peak mass values identical.)
peaks <- binPeaks(peaks, tolerance=0.002)

# prepare for statistical analysis (label mass spectra w infection outcome)
filenames <- sapply(peaks, function(x)metaData(x)$file[1])

cancer <- grepl(pattern="/tumor/", x=filenames)

classes <- factor(ifelse(cancer, "cancer", "control"),
                  levels=c("cancer", "control"))




# Normalisation ----
# median-based
# Normalisation ----
# Load and inspect
load("data/merged_data/baseline_corrected_spectra.RData")
str(baseline_corrected_spectra)


#: Median-based cross-normalisation BEFORE max-min scaling
cross_normalised_spectra <- lapply(baseline_corrected_spectra, function(spec) {
  med <- median(spec$intensity, na.rm = TRUE)
  if (!is.finite(med) || med == 0) med <- 1  # prevent division by zero
  list(
    mass = spec$mass,
    intensity = spec$intensity / med
  )
})

#  max-min re scaling normalisation

#rescale function
safe_rescale <- function(x) {
  r <- range(x, na.rm = TRUE, finite = TRUE)
  if (r[1] == r[2]) {
    warning("Flat spectrum detected — all intensities are equal.")  # to say if theres a dodgy spectrum (max and min same)
    return(rep(0, length(x)))
  } else {
    return((x - r[1]) / (r[2] - r[1]))
  }
}


scaled_spectra <- lapply(cross_normalised_spectra, function(spec) {
  list(
    mass = spec$mass,
    intensity = safe_rescale(spec$intensity)
  )
})

# Save scaled_spectra











# Load and inspect
load("data/merged_data/baseline_corrected_spectra.RData")
str(baseline_corrected_spectra)




# STEP 1: Median-based cross-normalisation BEFORE max-min scaling
cross_normalised_spectra <- lapply(baseline_corrected_spectra, function(spec) {
  med <- median(spec$intensity, na.rm = TRUE)
  if (!is.finite(med) || med == 0) med <- 1  # prevent division by zero
  list(
    mass = spec$mass,
    intensity = spec$intensity / med
  )
})




# STEP 2: Apply safe [0, 1] rescaling

# Safe rescale function
safe_rescale <- function(x) {
  r <- range(x, na.rm = TRUE, finite = TRUE)
  if (r[1] == r[2]) {
    warning("Flat spectrum detected — all intensities are equal.")
    return(rep(0, length(x)))
  } else {
    return((x - r[1]) / (r[2] - r[1]))
  }
}


scaled_spectra <- lapply(cross_normalised_spectra, function(spec) {
  list(
    mass = spec$mass,
    intensity = safe_rescale(spec$intensity)
  )
})












# Choose spectra to compare
spec_names <- c("ASARM100", "A115")  # add more if needed

# Define m/z window for zoom-in
mz_min <- 4000
mz_max <- 8000

# Loop through and plot
par(mfrow = c(length(spec_names), 1), mar = c(4, 4, 2, 1))  # stacked vertically

for (name in spec_names) {
  spec_raw <- baseline_corrected_spectra[[name]]
  spec_cross <- scaled_spectra[[name]]
  
  # Scale raw (direct)
  spec_raw_scaled <- list(
    mass = spec_raw$mass,
    intensity = safe_rescale(spec_raw$intensity)
  )
  
  # Subset for zoom window
  subset_raw <- with(spec_raw_scaled, mass >= mz_min & mass <= mz_max)
  subset_cross <- with(spec_cross, mass >= mz_min & mass <= mz_max)
  
  # Plot
  plot(spec_raw_scaled$mass[subset_raw], spec_raw_scaled$intensity[subset_raw],
       type = "l", col = "grey40", lwd = 1.2,
       xlab = "m/z", ylab = "Intensity (scaled)",
       main = paste(name, ": Raw vs Median-Normalised"))
  
  lines(spec_cross$mass[subset_cross], spec_cross$intensity[subset_cross], col = "blue", lwd = 1.2)
  legend("topright", legend = c("Raw scaled", "Median cross-normalised + scaled"),
         col = c("grey40", "blue"), lty = 1, lwd = 1.2, cex = 0.8)
}








# Peak selection ----
# Load data
load("data/merged_data/scaled_spectra.Rdata")
str(scaled_spectra)

# Example spectrum
spec <- scaled_spectra[["A060"]]
mz <- spec$mass
intensity <- spec$intensity



# 1. Load reticulate
library(reticulate)

# 2. Point to your virtual environment
use_virtualenv("~/.virtualenvs/r-reticulate", required = TRUE)

# 3. Check config (should match your Python 3.13)
py_config()

# 4. Clear Python environment (optional but good for clean state)
py_run_string("globals().clear()")

# 5. Prepare R data
spec <- scaled_spectra[["ASARM110"]]
mz <- spec$mass
intensity <- spec$intensity

# 6. Push data from R to Python
py$mz <- mz
py$intensity <- intensity

# 7. Run peak picking in Python
py_run_string("
import numpy as np
from scipy.signal import find_peaks

# Convert R vectors to numpy arrays
mz = np.array(mz)
intensity = np.array(intensity)

# Detect peaks
peaks, properties = find_peaks(intensity, prominence = 0.02, 
height = 0.02)

# Extract peak values
mz_fp = mz[peaks]
intensity_fp = intensity[peaks]
")

length(py$mz_fp)

# 0.02 for cc22

plot(mz, intensity, type = "l", col = "grey30", main = "Prominence-Based Peaks (A060)",
     xlab = "m/z", ylab = "Intensity")
points(reticulate::py$mz_fp, reticulate::py$intensity_fp, col = "purple", pch =19, cex = 0.8)







# Smaller scale
# Set zoom range
mz_min <- 2900
mz_max <- 3100

# Subset full spectrum
subset_idx <- which(mz >= mz_min & mz <= mz_max)
mz_zoom <- mz[subset_idx]
intensity_zoom <- intensity[subset_idx]

# Subset peak data to zoom range
peak_idx_in_zoom <- which(py$mz_fp >= mz_min & py$mz_fp <= mz_max)
mz_peaks_zoom <- py$mz_fp[peak_idx_in_zoom]
intensity_peaks_zoom <- py$intensity_fp[peak_idx_in_zoom]

# Plot
plot(mz_zoom, intensity_zoom, type = "l", col = "black",
     xlab = "m/z", ylab = "Intensity", main = paste("Peak Detection (", mz_min, "-", mz_max, "m/z)", sep = ""))
points(mz_peaks_zoom, intensity_peaks_zoom, col = "purple", pch = 19) # peaks



# Peak Alignment ----
#library(tibble)


library(ggplot2)
library(dplyr)
library(tibble)



# 1. extract the full scaled spectrum for A060
scaled_full <- scaled_spectra[["A060"]]
scaled_df   <- tibble(mass = scaled_full$mass, intensity = scaled_full$intensity)

# 2. your custom alignment function
align_peaks_custom <- function(peak_list, tol = 0.002) {
  all_mz <- sort(unlist(lapply(peak_list, `[[`, "mz")))
  bins <- c(); i <- 1
  while(i <= length(all_mz)) {
    start <- all_mz[i]; j <- i
    while(j+1 <= length(all_mz) && all_mz[j+1] - start <= tol) j <- j + 1
    bins <- c(bins, mean(all_mz[i:j])); i <- j + 1
  }
  iso_names <- names(peak_list)
  mat <- matrix(NA_real_, nrow = length(iso_names), ncol = length(bins),
                dimnames = list(iso_names, sprintf("%.4f", bins)))
  for(spec in iso_names) {
    mz_vec  <- peak_list[[spec]]$mz
    int_vec <- peak_list[[spec]]$intensity
    for(k in seq_along(mz_vec)) {
      idx <- which.min(abs(bins - mz_vec[k]))
      if (abs(bins[idx] - mz_vec[k]) <= tol)
        mat[spec, idx] <- int_vec[k]
    }
  }
  list(bins = bins, intensity_mat = mat)
}

# 3. pick only A060 peaks
single_peak_list <- list(A060 = peak_list[["A060"]])

# 4. tolerances to compare (ascending)
tols <- c(1, 5, 25, 50)

# 5. align and build data for each tol
aligned_dfs <- lapply(tols, function(tol) {
  res <- align_peaks_custom(single_peak_list, tol = tol)
  tibble(
    mass      = res$bins,
    intensity = as.numeric(res$intensity_mat["A060", ]),
    tol       = factor(paste0("tol=", tol), levels = paste0("tol=", tols))
  ) %>%
    filter(!is.na(intensity))
})

plot_data <- bind_rows(aligned_dfs)

# 6. final plot: full spectrum + peak overlays
ggplot() +
  # full scaled spectrum line
  geom_line(data = scaled_df,
            aes(x = mass, y = intensity),
            colour = "black", size = 0.5) +
  # aligned peaks
  geom_linerange(data = plot_data,
                 aes(x = mass, ymin = 0, ymax = intensity),
                 colour = "purple", alpha = 0.7) +
  geom_point(data = plot_data,
             aes(x = mass, y = intensity),
             colour = "purple", size = 1.5) +
  facet_wrap(~ tol, ncol = 1, scales = "free_y") +
  labs(
    title    = "Scaled Spectrum of A060 with Custom-Aligned Peaks",
    subtitle = "tolerance = 0.01, 0.05, 0.25, 0.50",
    x        = "m/z",
    y        = "Intensity"
  ) +
  theme_minimal(base_family = "Arial") +
  theme(strip.text = element_text(face = "bold"))











# Labelling ----
load("data/merged_data/aligned_peaks.RData")  # called mat
meta <- read.csv("data/merged_data/combined_metadata_filtered.csv")  # metadata

# Checking no isolates have been lost
length(meta$sampleID)  # same no. samples
nrow(mat)   # same no. samples

all(rownames(mat) == meta$sampleID)  # Same samples but in a different order

meta <- meta[match(rownames(mat), meta$sampleID), ]  #  Matches them up again
all(rownames(mat) == meta$sampleID)  # proof it's fixed

labels <- meta$Outcome30days  # Extract response variable

class(labels)


# Dimensionality Reduction ----
load("data/merged_data/aligned_peaks.RData")  # called mat

# total peaks before filtering
total_peaks <- ncol(mat)
cat("Total peaks before filtering:", total_peaks, "\n")

# minimum frequency threshold (reference peak present in at least 15% (42) of isolates )
threshold <- 0.15
min_freq <- threshold * nrow(mat) 

# no. isolates with a non-NA (i.e., present) value for each peak
peak_counts <- colSums(!is.na(mat))  

# keep only peaks above the threshold
filtered_matrix <- mat[, peak_counts >= min_freq]

# total peaks after filtering
filtered_peaks_no <- ncol(filtered_matrix)
cat("Peaks after frequency filtering (≥15% of isolates):", filtered_peaks_no, "\n")

#write.csv(filtered_matrix, "data/merged_data/filtered_matrix.csv", row.names = TRUE)


# UMAP
# install if needed
library(umap)
set.seed(1234)

# run UMAP
filtered_matrix_bin <- ifelse(is.na(filtered_matrix), 0, 1)  # remove na, and make binary

nrow(filtered_matrix_bin)

umap_res <- umap(filtered_matrix_bin)

plot(umap_res$layout)
# For mortality
plot(umap_res$layout, col = as.factor(meta$Outcome30days), pch = 19,
     main = "UMAP of Filtered Peaks (Mortality)", xlab = "UMAP1", ylab = "UMAP2")
legend("topright", legend = levels(as.factor(meta$Outcome30days)),
       col = 1:length(unique(meta$Outcome30days)), pch = 19)

# For MRSA/MSSA
plot(umap_res$layout, col = as.factor(meta$Collection), pch = 19,
     main = "UMAP of Filtered Peaks (MRSA/MSSA)", xlab = "UMAP1", ylab = "UMAP2")
legend("topright", legend = levels(as.factor(meta$Collection)),
       col = 1:length(unique(meta$Collection)), pch = 19)

# For Source
plot(umap_res$layout, col = as.factor(meta$source), pch = 19,
     main = "UMAP of Filtered Peaks (Source)", xlab = "UMAP1", ylab = "UMAP2")
legend("topright", legend = levels(as.factor(meta$source)),
       col = 1:length(unique(meta$source)), pch = 19)
  






# Spectral hump
load("data/merged_data/peak_list.RData")
load("data/merged_data/scaled_spectra.RData")
meta <- read.csv("data/merged_data/combined_metadata_filtered.csv")  # metadata

scaled <- scaled_spectra[["ASARM107"]]

ggplot() +
  geom_line(aes(x = scaled$mass, y = scaled$intensity), colour = "black") +
  labs(title = "Baseline Corrected & Scaled", x = NULL, y = "Intensity") +
  xlim(c(2000,8000)) +
  theme_minimal()


# Checking hump
# Find the column index corresponding to m/z 2186.0680
mz_value <- "2186.0680"
peak_col_index <- which(colnames(filtered_matrix) == mz_value)

# Check if it was found
if (length(peak_col_index) == 0) {
  stop("m/z 2186.0680 not found in column names.")
}

# Define a threshold for peak presence
threshold <- 0.000001  # Adjust as needed

# Count how many spectra have intensity above the threshold at that m/z
present_count <- sum(filtered_matrix[, peak_col_index] > threshold, na.rm = TRUE)

# Calculate percentage
total_spectra <- nrow(filtered_matrix)
percentage_present <- (present_count / total_spectra) * 100

# Output result
cat(sprintf("Peak at m/z 2186.0680 is present in %.1f%% of spectra.\n", percentage_present))






# Define the m/z value and threshold
mz_value <- "2186.0680"
threshold <- 0.05

# Get the column index for the peak
peak_col_index <- which(colnames(filtered_matrix) == mz_value)

# Identify sample names starting with "A" followed by digits (e.g., "ASARM123")
is_cc22 <- grepl("^A\\d+", rownames(filtered_matrix))

# Subset the matrix
cc22_matrix <- filtered_matrix[is_cc22, ]
cork_matrix <- filtered_matrix[!is_cc22, ]

# Replace NA with 0 in both subsets
cc22_peak <- cc22_matrix[, peak_col_index]
cc22_peak[is.na(cc22_peak)] <- 0

cork_peak <- cork_matrix[, peak_col_index]
cork_peak[is.na(cork_peak)] <- 0

# Calculate % of spectra with peak present
cc22_percent <- sum(cc22_peak > threshold) / length(cc22_peak) * 100
cork_percent <- sum(cork_peak > threshold) / length(cork_peak) * 100

# Output
cat(sprintf("Peak at m/z 2186.0680 is present in %.1f%% of cc22 spectra.\n", cc22_percent))
cat(sprintf("Peak at m/z 2186.0680 is present in %.1f%% of cork spectra.\n", cork_percent))





