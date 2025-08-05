library(ggplot2)
library(dplyr)
library(tibble)

# Load required data
load("data/merged_data/scaled_spectra.RData")
load("data/merged_data/peak_list.RData")

# 1. Extract the full scaled spectrum for A060
scaled_full <- scaled_spectra[["A060"]]
scaled_df   <- tibble(mass = scaled_full$mass, intensity = scaled_full$intensity)

# 2. Custom alignment function
align_peaks_custom <- function(peak_list, tol = 0.002) {
  all_mz <- sort(unlist(lapply(peak_list, `[[`, "mz")))
  bins <- c(); i <- 1
  while(i <= length(all_mz)) {
    start <- all_mz[i]; j <- i
    while(j + 1 <= length(all_mz) && all_mz[j + 1] - start <= tol) j <- j + 1
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

# 3. Pick only A060 peaks
single_peak_list <- list(A060 = peak_list[["A060"]])

# 4. Tolerances to compare (ascending)
tols <- c(1, 5, 25, 50)

# 5. Align and build data for each tolerance
aligned_dfs <- lapply(tols, function(tol) {
  res <- align_peaks_custom(single_peak_list, tol = tol)
  (tibble(
    mass      = res$bins,
    intensity = as.numeric(res$intensity_mat["A060", ]),
    tol       = factor(paste0("tol=", tol), levels = paste0("tol=", tols))
  )) %>%
    filter(!is.na(intensity))
})

# Combine all into one data frame
plot_data <- bind_rows(aligned_dfs)

# 6. final plot: full spectrum + peak overlays
# Final plot: Full scaled spectrum + peak overlays
peak_alignment_improved <- ggplot() +
  # Full scaled spectrum as a grey line
  geom_line(data = scaled_df,
            aes(x = mass, y = intensity),
            colour = "grey60", size = 0.5) +
  
  # Aligned peaks: vertical lines + points
  geom_linerange(data = plot_data,
                 aes(x = mass, ymin = 0, ymax = intensity),
                 colour = "purple", alpha = 0.5, linewidth = 1) +
  geom_point(data = plot_data,
             aes(x = mass, y = intensity),
             colour = "purple", size = 2.2) +
  
  # Facet by tolerance
  facet_wrap(~ tol, ncol = 1, scales = "free_y",
             labeller = as_labeller(function(x) paste0("Tolerance: ", gsub("tol=", "", x)))) +
  
  # Labels
  labs(
    title    = "Scaled Spectrum of A060 with Custom-Aligned Peaks",
    x        = "\nm/z",
    y        = "Relative Intensity\n"
  ) +
  
  theme_minimal(base_family = "Arial") +
  theme(
    plot.title = element_text(size = 24, face = "bold", hjust = 0),
    axis.title = element_text(size = 25, colour = "black"),
    axis.text  = element_text(size = 20, colour = "black"),
    strip.text = element_text(face = "bold", size = 20, hjust = 0),
    panel.grid = element_blank(),
    plot.background = element_rect(fill = "white", colour = NA)  # <- this line
  )



ggsave("code/code/figures/preprocessing/peak_alignment/peak_alignment_improved.png",
       plot = peak_alignment_improved, width = 20, height = 12, dpi = 300, bg = "white")





# Smaller range ----
# 1. extract scaled spectrum and restrict to 2900–3100
scaled_full <- scaled_spectra[["A060"]]
scaled_df <- tibble(
  mass = scaled_full$mass,
  intensity = scaled_full$intensity
) %>%
  filter(mass >= 2900, mass <= 3100)

# 2. custom alignment function (as before)
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
    mz_vec <- peak_list[[spec]]$mz
    int_vec <- peak_list[[spec]]$intensity
    for(k in seq_along(mz_vec)) {
      idx <- which.min(abs(bins - mz_vec[k]))
      if (abs(bins[idx] - mz_vec[k]) <= tol) 
        mat[spec, idx] <- int_vec[k]
    }
  }
  list(bins = bins, intensity_mat = mat)
}

# 3. build peak alignments for tolerances
single_peak_list <- list(A060 = peak_list[["A060"]])
tols <- c(1, 5, 10, 25)

aligned_dfs <- lapply(tols, function(tol) {
  res <- align_peaks_custom(single_peak_list, tol = tol)
  tibble(
    mass      = res$bins,
    intensity = as.numeric(res$intensity_mat["A060", ]),
    tol       = tol
  ) %>% 
    filter(!is.na(intensity), mass >= 2900, mass <= 3100)
})

plot_data <- bind_rows(aligned_dfs) %>%
  mutate(tol = factor(tol, levels = tols))

# 4. plot: scaled spectrum + peaks
peak_alignment_zoomed <- ggplot() +
  # scaled spectrum line (same for all panels)
  geom_line(data = scaled_df, aes(x = mass, y = intensity), colour = "black", size = 0.5) +
  # peaks
  geom_linerange(data = plot_data, aes(x = mass, ymin = 0, ymax = intensity),
                 colour = "purple", alpha = 0.7) +
  geom_point(data = plot_data, aes(x = mass, y = intensity),
             colour = "purple", size = 1.5) +
facet_wrap(~ tol, ncol = 1, scales = "free_y",
           labeller = as_labeller(function(x) paste0("Tolerance: ", x))) +
  labs(
    title    = "Scaled Spectrum (A060) with Aligned Peaks at Various Tolerances",
    subtitle = "m/z range 2900–3100",
    x        = "\nm/z",
    y        = "Relative Intensity\n",
    caption  = "Tolerance values shown in ascending order"
  ) +
  theme_minimal(base_family = "Arial") +
  theme(
    plot.title    = element_text(size = 24, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = 20, hjust = 0),
    plot.caption  = element_text(size = 16, hjust = 0),
    axis.title    = element_text(size = 25, colour = "black"),
    axis.text     = element_text(size = 20, colour = "black"),
    strip.text    = element_text(face = "bold", size = 20, hjust = 0),
    panel.grid    = element_blank(),
    plot.background = element_rect(fill = "white", colour = NA)  # solid white background
  )

# Save
ggsave("code/code/figures/preprocessing/peak_alignment/peak_alignment_zoomed.png",
       plot = peak_alignment_zoomed, width = 20, height = 12, dpi = 300, bg = "white")











# Plot for mutiple spectra 2900 - 3100


library(tidyverse)
library(tibble)

# 1. extract and trim scaled spectra for A060 and A090
isos <- c("A060", "ASARM103")
scaled_dfs <- map_dfr(isos, function(id) {
  tibble(
    iso       = id,
    mass      = scaled_spectra[[id]]$mass,
    intensity = scaled_spectra[[id]]$intensity
  ) %>%
    filter(mass >= 2950, mass <= 3050)
})

# 2. custom alignment (same as before)
align_peaks_custom <- function(peak_list, tol) {
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
  for(sp in iso_names) {
    mzv <- peak_list[[sp]]$mz
    iv  <- peak_list[[sp]]$intensity
    for(k in seq_along(mzv)) {
      idx <- which.min(abs(bins - mzv[k]))
      if (abs(bins[idx] - mzv[k]) <= tol)
        mat[sp, idx] <- iv[k]
    }
  }
  list(bins = bins, mat = mat)
}

# 3. Align for each tolerance value
tols <- c(1, 5, 10, 25)
aligned_list <- lapply(tols, function(tol) {
  res <- align_peaks_custom(peak_list[isos], tol = tol)
  bins <- res$bins
  mat  <- res$mat
  as_tibble(mat, rownames = "iso") %>%
    pivot_longer(-iso, names_to = "mass", values_to = "intensity") %>%
    mutate(
      mass = as.numeric(mass),
      tol  = tol
    ) %>%
    filter(!is.na(intensity), mass >= 2950, mass <= 3050)
})

plot_data <- bind_rows(aligned_list) %>%
  mutate(
    tol = factor(tol, levels = tols),
    iso = factor(iso, levels = isos)
  )

# 4. Final plot
(zoomed_alignment <- ggplot() +
  geom_line(data = scaled_dfs,
            aes(x = mass, y = intensity),
            colour = "black", size = 0.4) +
  geom_linerange(data = plot_data,
                 aes(x = mass, ymin = 0, ymax = intensity),
                 colour = "purple", alpha = 0.6) +
  geom_point(data = plot_data,
             aes(x = mass, y = intensity),
             colour = "purple", size = 1.2) +
  facet_grid(iso ~ tol, scales = "free_y",
             labeller = as_labeller(function(x) paste0("Tolerance: ", x))) +
  labs(
    x        = "\nm/z",
    y        = "Relative Intensity\n"
  ) +
  theme_minimal(base_family = "Arial") +
  theme(
    strip.text       = element_text(face = "bold", size = 18, hjust = 0),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title       = element_text(size = 22, colour = "black"),
    axis.text        = element_text(size = 16, colour = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title       = element_text(size = 24, face = "bold", hjust = 0),
    plot.subtitle    = element_text(size = 18, hjust = 0),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.spacing.y = unit(1.5, "lines"),
    panel.spacing.x = unit(1, "lines")
  ))


# Save
ggsave("code/code/figures/preprocessing/peak_alignment/zoomed_alignment.png",
       plot = zoomed_alignment, width = 20, height = 10, dpi = 300, bg = "white")



#  3 Spectra
# 1. extract and trim scaled spectra for A060, A090 & ASARM103
isos <- c("A060", "A090", "ASARM103")
scaled_dfs <- map_dfr(isos, function(id) {
  tibble(
    iso       = id,
    mass      = scaled_spectra[[id]]$mass,
    intensity = scaled_spectra[[id]]$intensity
  ) %>%
    filter(mass >= 2950, mass <= 3050)
})

# 2. custom alignment function (unchanged)
align_peaks_custom <- function(peak_list, tol) {
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
  for(sp in iso_names) {
    mzv <- peak_list[[sp]]$mz
    iv  <- peak_list[[sp]]$intensity
    for(k in seq_along(mzv)) {
      idx <- which.min(abs(bins - mzv[k]))
      if (abs(bins[idx] - mzv[k]) <= tol)
        mat[sp, idx] <- iv[k]
    }
  }
  list(bins = bins, mat = mat)
}

# 3. align for each tol and gather for plotting
tols <- c(1, 5, 10, 25)
aligned_list <- lapply(tols, function(tol) {
  res <- align_peaks_custom(peak_list[isos], tol = tol)
  as_tibble(res$mat, rownames = "iso") %>%
    pivot_longer(-iso, names_to = "mass", values_to = "intensity") %>%
    mutate(
      mass      = as.numeric(mass),
      tol       = paste0("tol=", tol)
    ) %>%
    filter(!is.na(intensity), mass >= 2950, mass <= 3050)
})
plot_data <- bind_rows(aligned_list) %>%
  mutate(
    tol = factor(tol, levels = paste0("tol=", tols)),
    iso = factor(iso, levels = isos)
  )

# 4. plot trimmed spectra + aligned peaks for 3 isolates
ggplot() +
  geom_line(data = scaled_dfs,
            aes(x = mass, y = intensity),
            colour = "black", size = 0.4) +
  geom_linerange(data = plot_data,
                 aes(x = mass, ymin = 0, ymax = intensity),
                 colour = "purple", alpha = 0.6) +
  geom_point(data = plot_data,
             aes(x = mass, y = intensity),
             colour = "purple", size = 1.2) +
  facet_grid(iso ~ tol, scales = "free_y") +
  labs(
    title    = "Trimmed (2950–3050 m/z) Scaled Spectra with Custom-Aligned Peaks",
    subtitle = "Isolates A060, A090 & ASARM103 across multiple tolerances",
    x        = "m/z",
    y        = "Intensity"
  ) +
  theme_minimal(base_family = "Arial") +
  theme(
    strip.text       = element_text(face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )



# Final alignment ----
# 0. prerequisites --------------------------------------------------------

# install.packages(c("tibble", "dplyr", "tidyr", "ggplot2"))  # if needed
library(tibble)
library(dplyr)
library(tidyr)
library(ggplot2)

# load your data
load("data/merged_data/scaled_spectra.RData")   # scaled_spectra: list of spectra
load("data/merged_data/peak_list.RData")        # peak_list: list of peak m/z + intensities

# 1. define custom alignment function -----------------------------------

align_peaks_custom <- function(peak_list, tol) {
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
  for(sp in iso_names) {
    mzv <- peak_list[[sp]]$mz
    iv  <- peak_list[[sp]]$intensity
    for(k in seq_along(mzv)) {
      idx <- which.min(abs(bins - mzv[k]))
      if (abs(bins[idx] - mzv[k]) <= tol)
        mat[sp, idx] <- iv[k]
    }
  }
  list(bins = bins, mat = mat)
}

# 2. run alignment at 5 Da threshold ------------------------------------

tol_value <- 5    # 5 m/z tolerance
res <- align_peaks_custom(peak_list, tol = tol_value)

# 3. checks and summaries -----------------------------------------------

# number of aligned feature bins
n_bins <- length(res$bins)
cat("number of aligned m/z bins:", n_bins, "\n")

# total non-NA peak assignments (i.e. total peaks across all isolates)
total_peaks <- sum(!is.na(res$mat))
cat("total peaks assigned across all spectra:", total_peaks, "\n")

# peaks per isolate
peaks_per_iso <- rowSums(!is.na(res$mat))
print(peaks_per_iso)
range(peaks_per_iso)

# 4. optional: quick visual check for one isolate (A060) ----------------

iso_id <- "A060"
df_iso <- tibble(
  mass      = res$bins,
  intensity = res$mat[iso_id, ]
) %>% filter(!is.na(intensity))

# overlay on the full scaled spectrum
scaled_iso <- tibble(
  mass      = scaled_spectra[[iso_id]]$mass,
  intensity = scaled_spectra[[iso_id]]$intensity
)

ggplot() +
  geom_line(data = scaled_iso, aes(x = mass, y = intensity),
            colour = "black", size = 0.5) +
  geom_linerange(data = df_iso,
                 aes(x = mass, ymin = 0, ymax = intensity),
                 colour = "purple", alpha = 0.7) +
  geom_point(data = df_iso,
             aes(x = mass, y = intensity),
             colour = "purple", size = 1.5) +
  labs(
    title = paste0("A060 scaled spectrum + aligned peaks (tol=", tol_value, " Da)"),
    x = "m/z", y = "Intensity"
  ) +
  theme_minimal(base_family = "Arial")


