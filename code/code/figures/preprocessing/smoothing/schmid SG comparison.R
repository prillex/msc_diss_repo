library(dplyr)
library(ggplot2)
library(signal)
library(patchwork)



# this requires the resampled spectra (run code in final preprocessing) ----

# Just raw spectra (RData with mass (to charge) and intensity)
load("data/merged_data/combined_spectra_raw_filtered.RData")  # called combined_spectra_raw_filtered


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


# SG Setup ----
generate_ms_kernel <- function(m, n = 6, alpha = 4) {  # m = window size, n = how fast the sinc function oscillates, alpha = Gaussian-based weighting 
  x <- seq(-m, m) / (m + 1)  # scaled window from -1 to 1
  w <- exp(-alpha * x^2) +   # gaussian centred at 0
    exp(-alpha * (x - 2)^2) +  # gaussian centred at +2
    exp(-alpha * (x + 2)^2) -  # gaussian centred at -2
    2 * exp(-alpha * 1^2)      # normalisation offset
  sinc <- function(x) ifelse(x == 0, 1, sin(pi * x) / (pi * x))  # sinc function
  kernel <- w * sinc(n * x)  # apply sinc modulation
  a <- 0.00172; b <- 0.02437; c <- 1.64375  # constants for correction
  correction <- a + b / (c - (m + 1))  # empirical sinusoidal correction
  kernel <- kernel + correction * sin((2 + 1) * pi * x)  # add correction
  kernel / sum(kernel)  # normalise
}

linear_extrapolate <- function(y, m) {
  left_lm <- lm(y[1:m] ~ seq_len(m))  # fit line to start of signal
  right_lm <- lm(y[(length(y) - m + 1):length(y)] ~ seq_len(m))  # fit line to end
  left_ext <- predict(left_lm, newdata = data.frame(seq_len = seq(-m, -1)))  # extrapolate left
  right_ext <- predict(right_lm, newdata = data.frame(seq_len = seq(length(y) + 1, length(y) + m)))  # extrapolate right
  c(left_ext, y, right_ext)  # combine padded signal
}

ms_smoother <- function(y, m = 30, n = 6) {  # Larger m = more smoothing, risk of blurring sharp peaks. Smaller m = keeps sharp peaks, but might not reduce noise well.
  kernel <- generate_ms_kernel(m, n)  # create kernel
  y_ext <- linear_extrapolate(y, m)  # pad signal
  conv <- stats::convolve(y_ext, rev(kernel), type = "filter")  # apply convolution
  conv[(m + 1):(length(conv) - m)]  # remove padding
}

pad_ms <- function(smoothed, m, original_len) {
  pad_left <- rep(NA, m)  # pad start
  pad_right <- rep(NA, original_len - length(smoothed) - m)  # pad end
  c(pad_left, smoothed, pad_right)  # combine
}



# Get data and apply smoothing ---
iso_A060 <- resampled_spectra[["A060"]]

# Raw
intensity_raw <- iso_A060$intensity

# SG (standard, from signal pkg)
intensity_sg <- sgolayfilt(intensity_raw, p = 4, n = 51)  # window size must be odd

# MS
intensity_ms <- ms_smoother(intensity_raw, m = 30, n = 6)
intensity_ms_padded <- pad_ms(intensity_ms, m = 30, original_len = nrow(iso_A060))

# Plotting ----

df_plot <- iso_A060 %>%
  mutate(
    SG = intensity_sg,
    MS = intensity_ms_padded
  )

p1 <- ggplot(df_plot, aes(x = mass)) +
  geom_line(aes(y = intensity), color = "grey60", alpha = 0.4) +
  geom_line(aes(y = SG), color = "darkgreen") +
  labs(title = "Standard Savitzky–Golay", x = "m/z", y = "Intensity") +
  xlim(c(2000,8000))+
  theme_minimal()


p3 <- ggplot(df_plot, aes(x = mass)) +
  geom_line(aes(y = intensity), color = "grey60", alpha = 0.4) +
  geom_line(aes(y = MS), color = "darkorange") +
  labs(title = "Modified Sinc Kernel (MS)", x = "m/z", y = "Intensity") +
  xlim(c(2000,8000))+
  theme_minimal()

# Stack the plots
(p1 /  p3) + plot_layout(guides = "collect")








# extra bits on ms kernel ----

# (Chat GPT) This shows the need for normalisation of the kernel weight
# Reuse your kernel generator but *without* normalization for raw version
generate_ms_kernel_raw <- function(m, n = 6, alpha = 4) {
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
  tibble(x = x, kernel = kernel)
}

# Normalised version
generate_ms_kernel_normalized <- function(m, n = 6, alpha = 4) {
  k <- generate_ms_kernel_raw(m, n, alpha)
  k %>% mutate(kernel = kernel / sum(kernel))
}

# Plot both
library(ggplot2)
library(dplyr)

m <- 30
raw_k <- generate_ms_kernel_raw(m)
norm_k <- generate_ms_kernel_normalized(m)

raw_k$version <- "Unnormalized"
norm_k$version <- "Normalized"

ggplot(bind_rows(raw_k, norm_k), aes(x = x, y = kernel, color = version)) +
  geom_line(size = 0.8) +
  labs(title = "Effect of Normalizing the Modified Sinc Kernel",
       x = "x (normalized position in kernel)",
       y = "Kernel weight",
       color = "Kernel Version") +
  theme_minimal()





# Different kernels

# Custom kernel generator
generate_custom_kernel <- function(m, n, alpha) {
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
  tibble(x = x, kernel = kernel)
}

# Create all kernels with corresponding labels
kernels <- bind_rows(
  generate_custom_kernel(30, 6, 4)   %>% mutate(Label = "A) alpha = 4, n = 6"),
  generate_custom_kernel(30, 6, 10)  %>% mutate(Label = "B) alpha = 10, n = 6"),
  generate_custom_kernel(30, 12, 4)  %>% mutate(Label = "C) alpha = 4, n = 12"),
  generate_custom_kernel(30, 1, 0.5) %>% mutate(Label = "D) alpha = 0.5, n = 1")
)

# Define colour mapping
colour_map <- c(
  "A) alpha = 4, n = 6"   = "darkblue",
  "B) alpha = 10, n = 6"  = "darkgreen",
  "C) alpha = 4, n = 12"  = "darkorange",
  "D) alpha = 0.5, n = 1" = "purple"
)

# Final plot
ggplot(kernels, aes(x = x, y = kernel, colour = Label)) +
  geom_line(size = 1) +
  facet_wrap(~ Label, ncol = 2) +
  scale_colour_manual(values = colour_map, guide = "none") +
  labs(
    x = "Normalised Position (x)",
    y = "Kernel Weight"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 18, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 12)
  )

# Save
ggsave(
  ms_kernel_facet_plot,
  filename = "code/code/figures/preprocessing/smoothing/ms_kernel_facet_plot.png",
  bg = "white",
  width = 28, height = 20, dpi = 300
)






# Final plots ----
# Recreate p1 (Savitzky–Golay)
p1_fmt <- ggplot(df_plot, aes(x = mass)) +
  geom_line(aes(y = intensity), colour = "grey60", alpha = 0.4) +
  geom_line(aes(y = SG), colour = "darkgreen") +
  labs(title = "Savitzky–Golay Smoothing", x = "\nm/z", y = "Intensity\n") +
  xlim(c(2000, 8000)) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(colour = "black", size = 40),
    axis.title.y = element_text(colour = "black", size = 40),
    axis.text.x = element_text(colour = "black", size = 30),
    axis.text.y = element_text(colour = "black", size = 30),
    plot.title = element_text(size = 36, hjust = 0)
  )

# Recreate p3 (MS)
p3_fmt <- ggplot(df_plot, aes(x = mass)) +
  geom_line(aes(y = intensity), colour = "grey60", alpha = 0.4) +
  geom_line(aes(y = MS), colour = "darkorange") +
  labs(title = "Modified Sinc Kernel (MS)", x = "\nm/z", y = "Intensity\n") +
  xlim(c(2000, 8000)) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(colour = "black", size = 40),
    axis.title.y = element_text(colour = "black", size = 40),
    axis.text.x = element_text(colour = "black", size = 30),
    axis.text.y = element_text(colour = "black", size = 30),
    plot.title = element_text(size = 36, hjust = 0)
  )

# Stack plots
smoothing_comparison_plot <- (p1_fmt / p3_fmt) + plot_layout(guides = "collect")


# Save the plot
ggsave(
  smoothing_comparison_plot,
  filename = "code/code/figures/preprocessing/smoothing/smoothing_comparison_plot.png",
  bg = "white",
  width = 30, height = 14, dpi = 300, limitsize = FALSE
)










