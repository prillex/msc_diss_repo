# Matt Prill
# MSc Applied Data Science
# Dimensionality Reduction Pipeline Figure


# Load packages
library(tidyverse)
library(patchwork)

# === 1. Load all necessary data ===

# Load filtered peak matrix
load("data/merged_data/filtered_matrix.RData")

# Load uncorrelated matrix
load("data/merged_data/uncorrelated_matrix.Rdata")

# Load peak list
load("data/merged_data/peak_list.RData")

# Load RFECV-selected features (from previously written CSV)
df_rfecv <- read.csv("data/merged_data/nine_peaks.csv")

# --- Common Setup ---
library(tidyverse)

# A060 peaks
a060_peaks <- peak_list[["A060"]]
mz_orig <- a060_peaks$mz
int_orig <- a060_peaks$intensity

# Matching function (nearest m/z within tolerance)
get_matched_peaks <- function(target_mz, mz_source, int_source, tol = 1.0) {
  matched <- sapply(target_mz, function(mz_val) {
    diffs <- abs(mz_source - mz_val)
    if (min(diffs) <= tol) {
      return(int_source[which.min(diffs)])
    } else {
      return(NA)
    }
  })
  df <- tibble(mz = target_mz, intensity = matched) %>% drop_na()
  message(glue::glue("Matched {nrow(df)} of {length(target_mz)} peaks (tol = {tol})"))
  return(df)
}


# --- 1. Selected Peaks (original from peak_list) ---
df_peaks <- tibble(mz = mz_orig, intensity = int_orig)

# --- 2. Frequency-filtered (filtered_matrix) ---
mz_filtered <- as.numeric(colnames(filtered_matrix))
int_filtered <- as.numeric(filtered_matrix["A060", ])
df_filtered <- tibble(mz = mz_filtered, intensity = int_filtered) %>%
  drop_na()

# --- 3. Correlation-filtered (uncorrelated_matrix) ---
mz_uncor <- as.numeric(colnames(uncorrelated_matrix))
df_uncor <- get_matched_peaks(mz_uncor, mz_orig, int_orig)

# --- 4. RFECV (only selected peaks for A060) ---
mz_rfecv <- df_rfecv[1, -1]  # exclude ID column
names(mz_rfecv) <- gsub("^X", "", names(mz_rfecv))
mz_vals <- as.numeric(names(mz_rfecv))[mz_rfecv == 1]

df_rfecv_peaks <- get_matched_peaks(mz_vals, mz_orig, int_orig)


# --- Define Plot Theme ---
theme_spec <- theme(
  axis.title.x = element_text(colour = "black", size = 70),
  axis.title.y = element_text(colour = "black", size = 70),
  axis.text.x = element_text(colour = "black", size = 50),
  axis.text.y = element_text(colour = "black", size = 50),
  axis.line = element_line(colour = "black", linewidth = 3),
  panel.background = element_rect(fill = "white"),
  plot.background = element_rect(fill = "white", colour = NA),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  plot.margin = unit(c(1, 1, 1, 1), "cm"),
  legend.key.size = unit(3, 'cm'), #change legend key size
  legend.key.height = unit(2, 'cm'), #change legend key height
  legend.key.width = unit(5, 'cm'), #change legend key width
  legend.title = element_text(size = 20), #change legend title font size
  legend.text = element_text(size = 20),
  plot.title = element_text(size = 70, hjust = 0.5)
)

ylim_shared <- range(
  df_peaks$intensity,
  df_filtered$intensity,
  df_uncor$intensity,
  df_rfecv_peaks$intensity,
  na.rm = TRUE
)

# --- Build Plots ---
p1 <- ggplot(df_peaks, aes(x = mz, y = intensity)) +
  geom_segment(aes(xend = mz, y = 0, yend = intensity), color = "darkblue", alpha = 0.7) +
  geom_point(color = "darkblue", size = 3) +
  xlim(c(2000, 10000)) +
  coord_cartesian(ylim = ylim_shared) +
  labs(title = "Selected Peaks", x = NULL, y = "Relative Intensity\n") +
  theme_spec

p2 <- ggplot(df_filtered, aes(x = mz, y = intensity)) +
  geom_segment(aes(xend = mz, y = 0, yend = intensity), color = "darkorange", alpha = 0.7) +
  geom_point(color = "darkorange", size = 3) +
  xlim(c(2000, 10000)) +
  coord_cartesian(ylim = ylim_shared) +
  labs(title = "Frequency Filtered", x = NULL, y = "Relative Intensity\n") +
  theme_spec

p3 <- ggplot(df_uncor, aes(x = mz, y = intensity)) +
  geom_segment(aes(xend = mz, y = 0, yend = intensity), color = "purple", alpha = 0.7) +
  geom_point(color = "purple", size = 3) +
  xlim(c(2000, 10000)) +
  coord_cartesian(ylim = ylim_shared) +
  labs(title = "Correlation Filtered", x = NULL, y = "Relative Intensity\n") +
  theme_spec

p4 <- ggplot(df_rfecv_peaks, aes(x = mz, y = intensity)) +
  geom_segment(aes(xend = mz, y = 0, yend = intensity), color = "darkgreen", alpha = 0.7) +
  geom_point(color = "darkgreen", size = 3) +
  xlim(c(2000, 10000)) +
  coord_cartesian(ylim = ylim_shared) +
  labs(title = "Post-RFECV", x = "m/z", y = "Relative Intensity\n") +
  theme_spec

# --- Combine and Save ---
# Must save plots separately for formatting's sake


# Selected
ggsave(
  filename = "code/code/figures/preprocessing/peak_selection/p1.png",
  plot = p1,
  bg = "white",
  width = 36, height = 12, dpi = 300, limitsize = FALSE)

# Freq filtered
ggsave(
  filename = "code/code/figures/preprocessing/peak_selection/p2.png",
  plot = p2,
  bg = "white",
  width = 36, height = 12, dpi = 300, limitsize = FALSE)

# Correlation filtered
ggsave(
  filename = "code/code/figures/preprocessing/peak_selection/p3.png",
  plot = p3,
  bg = "white",
  width = 36, height = 12, dpi = 300, limitsize = FALSE)

# RFECV filtered
ggsave(
  filename = "code/code/figures/preprocessing/peak_selection/p4.png",
  plot = p4,
  bg = "white",
  width = 36, height = 12, dpi = 300, limitsize = FALSE)

