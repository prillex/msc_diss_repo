# Matt Prill
# MSc Applied Data Science Dissertation
# Bayesian Models

# Seed ----
set.seed(1234)

# Libraries ----
library(brms)
library(dplyr)
library(caret)
library(pROC)
library(tibble)
library(beepr)


# Bayes Models (Cross Validation) ----
# Data Prep ----
modelling_data <- read.csv("data/merged_data/modelling_data.csv", check.names = FALSE)
modelling_data$OutcomeBinary <- factor(modelling_data$Outcome30days == "Dead", levels = c(FALSE, TRUE), labels = c("Alive", "Dead"))
modelling_data <- modelling_data %>% mutate(Gender_Male = ifelse(Gender == "Male", 1, 0))

# Prepare metadata-only dataframe
meta_df <- modelling_data %>%
  select(OutcomeBinary, Age, Gender_Male) %>%
  na.omit()

# Prepare metadata + peaks dataframe
peak_cols <- grep("^\\d+$", colnames(modelling_data), value = TRUE)
new_peak_names <- paste0("peak_", peak_cols)
colnames(modelling_data)[colnames(modelling_data) %in% peak_cols] <- new_peak_names
peak_cols <- new_peak_names

bayes_df <- modelling_data %>%
  select(all_of(c("OutcomeBinary", "Age", "Gender_Male", peak_cols))) %>%
  na.omit()


# Cross Validation Folds ----
folds_meta_df <- caret::createFolds(meta_df$OutcomeBinary, k = 5, returnTrain = TRUE)
folds_bayes_df <- caret::createFolds(bayes_df$OutcomeBinary, k = 5, returnTrain = TRUE)

# ----
# ----
# ----
# Priors ----
# Uninformative Priors 
uninformative_priors <- c(
  prior(normal(0, 10), class = "b"),
  prior(normal(0, 10), class = "Intercept")
)

informative_priors <- c(
  prior(normal(0.03, 0.01), class = "b", coef = "Age"),
  prior(normal(-1.35, 0.5), class = "b", coef = "Gender_Male"),
  prior(normal(-1.51, 0.061), class = "Intercept")
)

# Function Bayes Models ----
run_cv_bayes <- function(folds, data, priors, label) {
  results <- lapply(folds, function(train_idx) {
    train_data <- data[train_idx, ]
    test_data  <- data[-train_idx, ]
    
    model <- brm(
      formula = OutcomeBinary ~ .,
      data = train_data,
      family = bernoulli(link = "logit"),
      prior = priors,
      chains = 3, cores = 3,
      iter = 2000, warmup = 200,
      seed = 123
    )
    
    probs <- fitted(model, newdata = test_data)[, "Estimate"]
    preds <- ifelse(probs > 0.5, "Dead", "Alive")
    
    list(
      preds = preds,
      probs = probs,
      truth = test_data$OutcomeBinary
    )
  })
  
  list(
    preds = unlist(lapply(results, \(x) x$preds)),
    probs = unlist(lapply(results, \(x) x$probs)),
    truth = unlist(lapply(results, \(x) as.character(x$truth)))
  )
}


# Bayesian Models ----
cv_meta_uninf     <- run_cv_bayes(folds_meta_df, meta_df, uninformative_priors, "Meta: Uninformative")
cv_meta_inf       <- run_cv_bayes(folds_meta_df, meta_df, informative_priors, "Meta: Informative")
cv_peaks_uninf    <- run_cv_bayes(folds_bayes_df, bayes_df, uninformative_priors, "Peaks: Uninformative")
cv_peaks_inf      <- run_cv_bayes(folds_bayes_df, bayes_df, informative_priors, "Peaks: Informative")

# ----
# ----
# ----
# Function Evaluation ----
evaluate_cv_model <- function(preds, probs, truth) {
  cm <- caret::confusionMatrix(
    factor(preds, levels = c("Alive", "Dead")),
    factor(truth, levels = c("Alive", "Dead")),
    positive = "Dead"
  )
  
  truth_bin <- ifelse(truth == "Dead", 1, 0)
  roc_obj <- pROC::roc(truth_bin, probs)
  auc_val <- as.numeric(pROC::auc(roc_obj))
  
  tibble(
    Accuracy = cm$overall["Accuracy"],
    F1_Dead  = cm$byClass["F1"],
    AUC      = auc_val
  )
}

# Evaluate Bayesian Models ----
results_meta_uninf     <- evaluate_cv_model(cv_meta_uninf$preds, cv_meta_uninf$probs, cv_meta_uninf$truth)
results_meta_inf       <- evaluate_cv_model(cv_meta_inf$preds, cv_meta_inf$probs, cv_meta_inf$truth)
results_peaks_uninf    <- evaluate_cv_model(cv_peaks_uninf$preds, cv_peaks_uninf$probs, cv_peaks_uninf$truth)
results_peaks_inf      <- evaluate_cv_model(cv_peaks_inf$preds, cv_peaks_inf$probs, cv_peaks_inf$truth)

# Combine into summary table
cv_results_summary <- bind_rows(
  results_meta_uninf,
  results_meta_inf,
  results_peaks_uninf,
  results_peaks_inf
) %>%
  mutate(Model = c("Meta: Uninformative", "Meta: Informative", "Peaks: Uninformative", "Peaks: Informative")) %>%
  select(Model, Accuracy, F1_Dead, AUC)

# Print
print(cv_results_summary)
# Convergence Checks ----
# Posterior Predictive Checks
pp_check(cv_meta_uninf)         # Meta: Uninformative
pp_check(cv_meta_inf)           # Meta: Informative
pp_check(cv_peaks_uninf)        # Peaks: Uninformative
pp_check(cv_peaks_inf)          # Peaks: Informative

# Caterpillar Plots
plot(cv_meta_uninf)
plot(cv_meta_inf)
plot(cv_peaks_uninf)
plot(cv_peaks_inf)






# ----
# ----
# ----
# Credibility Interval Plot (Model 4) ----
# Load required libraries
library(brms)
library(ggplot2)
library(dplyr)
library(tidyr)
library(posterior)
library(stringr)

# full model
bayes_peaks <- readRDS("code/code/analysis/models/bayes_peaks.rds")


# Extract posterior draws
bayes_peaks_ci <- as_draws_df(bayes_peaks) %>%   # UPDATE to your final model name
  select(starts_with("b_")) %>%
  pivot_longer(everything(), names_to = "term", values_to = "value") %>%
  group_by(term) %>%
  summarise(
    median = median(value),
    lower  = quantile(value, 0.025),
    upper  = quantile(value, 0.975),
    .groups = "drop"
  )

# Clean variable names
bayes_peaks_ci <- bayes_peaks_ci %>%
  mutate(
    term_clean = term %>%
      str_replace("^b_", "") %>%
      str_replace_all("_", " ") %>%
      str_replace("^peak", "Peak") %>%
      str_replace("^Gender male$", "Gender male") %>%
      str_replace("^Age$", "Age (Years)") %>%
      str_replace("^Intercept$", "Intercept")
  )

# Order terms by effect size
bayes_peaks_ci$term_clean <- factor(bayes_peaks_ci$term_clean, levels = bayes_peaks_ci$term_clean[order(bayes_peaks_ci$median)])

# Whisker width
cap_width <- 0.15

# Plot
CI_plot <- ggplot(bayes_peaks_ci, aes(x = term_clean, y = median)) +
  geom_segment(aes(x = term_clean, xend = term_clean, y = lower, yend = upper),
               linewidth = 1, colour = "darkblue") +
  geom_segment(aes(x = as.numeric(term_clean) - cap_width,
                   xend = as.numeric(term_clean) + cap_width,
                   y = lower, yend = lower),
               colour = "darkblue", linewidth = 1) +
  geom_segment(aes(x = as.numeric(term_clean) - cap_width,
                   xend = as.numeric(term_clean) + cap_width,
                   y = upper, yend = upper),
               colour = "darkblue", linewidth = 1) +
  geom_point(size = 3, colour = "darkblue") +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "black", linewidth = 0.9) +
  coord_flip() +
  labs(
    title = "95% Credible Intervals for Model Coefficients",
    x = "Predictor\n",
    y = "\nLog-Odds Estimate"
  ) +
  theme_minimal(base_family = "Arial") +
  theme(
    axis.title.x = element_text(colour = "black", size = 25),
    axis.title.y = element_text(colour = "black", size = 25),
    axis.text.x  = element_text(colour = "black", size = 20),
    axis.text.y  = element_text(colour = "black", size = 20),
    axis.line    = element_line(colour = "black", linewidth = 1.5),
    axis.ticks   = element_line(colour = "black", linewidth = 1.5),
    axis.ticks.length = unit(0.25, "cm"),
    plot.title   = element_text(size = 18, face = "bold", colour = "black"),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# Display
print(CI_plot)

# Save
ggsave("code/code/figures/outcome/CI_plot.png",
       plot = CI_plot, width = 14, height = 12, dpi = 300)


# Save using ggsave with wrap_plots()
#ggsave("code/code/figures/outcome/CI_int.png",
#       plot = CI_int,
#       width = 14, height = 12, dpi = 300)
